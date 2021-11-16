# -*- coding: utf-8 -*-

"""Processor for the INDRA database."""

import gzip
import json
import logging
import pickle
import click 
from more_click import verbose_option
from pathlib import Path
from tqdm import tqdm
from typing import Tuple, Union, Optional

import humanize
import pandas as pd
import pystow

from indra.ontology.bio import bio_ontology
from indra.databases.identifiers import ensure_prefix_if_needed
from indra_cogex.representation import Node, Relation
from indra_cogex.sources.processor import Processor

logger = logging.getLogger(__name__)
tqdm.pandas()


# If you don't have the data, get it from:
# 's3://bigmech/indra-db/dumps/2021-01-26/sif.pkl'


class DbProcessor(Processor):
    """Processor for the INDRA database."""

    name = "database"
    descriptions = {
        "indra_rel": "The source and target come from a variety of INDRA statements."
    }
    df: pd.DataFrame

    def __init__(self, dir_path: Union[None, str, Path] = None,
                 add_jsons: Optional[bool] = False):
        """Initialize the INDRA database processor.

        Parameters
        ----------
        dir_path :
            The path to the directory containing INDRA database SIF dump pickle
            and batches of statements (stored in batch*.json.gz files, only if
            add_jsons=True). If none given, will look in the default location.
        add_jsons :
            Whether to include statements JSONs in relation data. Default: False.
        """
        if dir_path is None:
            dir_path = pystow.join("indra", "db")
        elif isinstance(dir_path, str):
            dir_path = Path(dir_path)
        if add_jsons:
            self.stmt_fnames = dir_path.glob("batch*.json.gz")
            logger.info("Creating DB with Statement JSONs")
        else:
            self.stmt_fnames = None
            logger.info("Creating DB without Statement JSONs")
        sif_path = dir_path.joinpath("sif.pkl")
        with open(sif_path, "rb") as fh:
            logger.info("Loading %s" % sif_path)
            df = pickle.load(fh)
        logger.info("Loaded %s rows from %s", humanize.intword(len(df)), sif_path)
        self.df = df
        logger.info("Fixing ID and naming issues...")
        for side in "AB":
            # A lot of the names in the SIF dump are all over
            self.df[f"ag{side}_name"] = [
                bio_ontology.get_name(prefix, identifier)
                for prefix, identifier in self.df[
                    [f"ag{side}_ns", f"ag{side}_id"]
                ].values
            ]
            self.df[f"ag{side}_ns"], self.df[f"ag{side}_id"] = list(
                zip(
                    *[
                        fix_id(db_ns, db_id)
                        for db_ns, db_id in tqdm(
                            zip(list(df[f"ag{side}_ns"]), list(df[f"ag{side}_id"])),
                            total=len(df),
                            desc="Fixing IDs",
                        )
                    ]
                )
            )
        self.df["source_counts"] = self.df["source_counts"].apply(json.dumps)
        self.df = self.df.dropna(subset=['belief'])

    def get_nodes(self):  # noqa:D102
        df = pd.concat(
            [
                self.df[["agA_ns", "agA_id", "agA_name"]].rename(
                    columns={"agA_ns": "ns", "agA_id": "id", "agA_name": "name"}
                ),
                self.df[["agB_ns", "agB_id", "agB_name"]].rename(
                    columns={"agB_ns": "ns", "agB_id": "id", "agB_name": "name"}
                ),
            ],
            ignore_index=True,
        ).drop_duplicates()
        for db_ns, db_id, name in df.values:
            yield Node(db_ns, db_id, ["BioEntity"], dict(name=name))

    def get_relations(self):  # noqa:D102
        rel_type = "indra_rel"
        columns = [
            "agA_ns",
            "agA_id",
            "agB_ns",
            "agB_id",
            "stmt_type",
            "source_counts",
            "evidence_count",
            "belief",
            "stmt_hash",
        ]
        total_count = 0
        # If we want to add statement JSONs, process the statement batches and
        # map to records in SIF dataframe
        if self.stmt_fnames:
            # Remove duplicate hashes (e.g. reverse edges for Complexes)
            df = self.df.drop_duplicates(subset='stmt_hash', keep='first')
            # Convert to dict with hashes as keys
            df = df.set_index('stmt_hash')
            df_dict = df.to_dict(orient='index')
            for fname in tqdm(self.stmt_fnames):
                count = 0
                with gzip.open(fname, 'r') as fh:
                    # For each statement find corresponding row in df
                    for i, line in enumerate(fh.readlines()):
                        stmt = json.loads(line)
                        stmt_hash = int(stmt['matches_hash'])
                        try:
                            values = df_dict[stmt_hash]
                            data = {
                                "stmt_hash:long": stmt_hash,
                                "source_counts:string": values['source_counts'],
                                "evidence_count:int": values['evidence_count'],
                                "stmt_type:string": values['stmt_type'],
                                "belief:float": values['belief'],
                                "stmt_json:string": json.dumps(stmt),
                            }
                            count += 1
                            yield Relation(
                                values['agA_ns'],
                                values['agA_id'],
                                values['agB_ns'],
                                values['agB_id'],
                                rel_type,
                                data,
                            )
                        # This statement is not in df
                        except KeyError:
                            continue
                total_count += count
                logger.info(f'Got {count} relations from {i} records in {fname}')
        # Otherwise only process the SIF dataframe
        else:
            for (
                source_ns,
                source_id,
                target_ns,
                target_id,
                stmt_type,
                source_counts,
                evidence_count,
                belief,
                stmt_hash,
            ) in (
                self.df[columns].drop_duplicates().values
            ):
                data = {
                    "stmt_hash:long": stmt_hash,
                    "source_counts:string": source_counts,
                    "evidence_count:int": evidence_count,
                    "stmt_type:string": stmt_type,
                    "belief:float": belief,
                }
                total_count += 1
                yield Relation(
                    source_ns,
                    source_id,
                    target_ns,
                    target_id,
                    rel_type,
                    data,
                )
        logger.info(f'Got {total_count} total relations')         

    @classmethod
    def get_cli(cls) -> click.Command:
        """Get the CLI for this processor."""

        # Add custom option not available in other processors
        @click.command()
        @click.option("--add_jsons", is_flag=True)
        @verbose_option
        def _main(add_jsons: bool):
            click.secho(f"Building {cls.name}", fg="green", bold=True)
            processor = cls(add_jsons=add_jsons)
            processor.dump()

        return _main

def fix_id(db_ns: str, db_id: str) -> Tuple[str, str]:
    """Fix ID issues specific to the SIF dump."""
    if db_ns == "GO":
        if db_id.isnumeric():
            db_id = "0" * (7 - len(db_id)) + db_id
    if db_ns == "EFO" and db_id.startswith("EFO:"):
        db_id = db_id[4:]
    if db_ns == "UP" and db_id.startswith("SL"):
        db_ns = "UPLOC"
    if db_ns == "UP" and "-" in db_id and not db_id.startswith("SL-"):
        db_id = db_id.split("-")[0]
    if db_ns == "FPLX" and db_id == "TCF-LEF":
        db_id = "TCF_LEF"
    db_id = ensure_prefix_if_needed(db_ns, db_id)
    return db_ns, db_id
