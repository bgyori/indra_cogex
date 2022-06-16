# -*- coding: utf-8 -*-

"""Processor for the INDRA database."""

import codecs
import csv
import gzip
import json
import logging
import os
import pickle
import textwrap
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

import click
import humanize
import pandas as pd
import pystow
from indra.databases.identifiers import ensure_prefix_if_needed
from indra.ontology.bio import bio_ontology
from indra.util import batch_iter
from indra.util.statement_presentation import (
    db_sources,
    internal_source_mappings,
    reader_sources,
)
from more_click import verbose_option
from tqdm import tqdm

from indra_cogex.representation import Node, Relation
from indra_cogex.sources.processor import Processor

logger = logging.getLogger(__name__)
tqdm.pandas()


# If you don't have the data, get it from:
# 's3://bigmech/indra-db/dumps/2021-01-26/sif.pkl'


class DbProcessor(Processor):
    """Processor for the INDRA database."""

    name = "database"
    df: pd.DataFrame
    node_types = ["BioEntity"]

    def __init__(
        self, dir_path: Union[None, str, Path] = None, add_jsons: Optional[bool] = False
    ):
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
            self.stmts_fname = dir_path / "statements_with_evidences.tsv.gz"
            self.text_refs_fname = dir_path / "text_refs_for_reading.tsv.gz"
            logger.info(
                "Creating DB with Statement JSONs. Note that this "
                "requires at least 128GB of RAM"
            )
        else:
            self.stmts_fname = None
            self.text_refs_fname = None
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
        self.df = self.df.dropna(subset=["belief"])

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
        if self.stmts_fname:
            ensure_statements_with_evidences(self.stmts_fname.as_posix())
            ensure_text_refs_for_reading(self.text_refs_fname.as_posix())
            # Remove duplicate hashes (e.g. reverse edges for Complexes)
            df = self.df.drop_duplicates(subset="stmt_hash", keep="first")
            # Convert to dict with hashes as keys
            df = df.set_index("stmt_hash")
            df_dict = df.to_dict(orient="index")
            hashes_yielded = set()
            logger.info("Getting text refs from text refs file")
            text_refs = load_text_refs_for_reading_dict(self.text_refs_fname)
            with gzip.open(self.stmts_fname, "rt", encoding="utf-8") as fh:
                # For each statement find corresponding row in df
                reader = csv.reader(fh, delimiter="\t")
                for (
                    raw_stmt_id,
                    reading_id,
                    stmt_hash,
                    raw_json_str,
                    pa_json_str,
                ) in reader:
                    stmt_hash = int(stmt_hash)
                    # If we already yielded this statement, we can skip it
                    if stmt_hash in hashes_yielded:
                        continue
                    stmt_json = load_statement_json(pa_json_str)
                    try:
                        values = df_dict[stmt_hash]
                        sources = set(json.loads(values["source_counts"]))
                        # For statements with only evidence from medscan,
                        # we don't add an evidence and yield the statement
                        medscan_only = sources == {"medscan"}
                        if medscan_only:
                            stmt_json["evidence"] = []
                        # Otherwise, we know that eventually we will bump into
                        # an evidence we can use and so we skip any medscan
                        # ones without yielding the statement
                        else:
                            raw_json = load_statement_json(raw_json_str)
                            raw_json_ev = raw_json["evidence"][0]
                            if raw_json_ev["source_api"] == "medscan":
                                continue
                            elif reading_id != "\\N":
                                tr = text_refs[reading_id]
                                raw_json_ev["text_refs"] = tr
                                if "PMID" in raw_json_ev["text_refs"]:
                                    raw_json_ev["pmid"] = raw_json_ev["text_refs"][
                                        "PMID"
                                    ]
                            stmt_json["evidence"] = raw_json["evidence"]
                        data = {
                            "stmt_hash:long": stmt_hash,
                            "source_counts:string": values["source_counts"],
                            "evidence_count:int": values["evidence_count"],
                            "stmt_type:string": values["stmt_type"],
                            "belief:float": values["belief"],
                            "stmt_json:string": json.dumps(stmt_json),
                            "has_database_evidence:bool": any(
                                source in db_sources for source in sources
                            ),
                            "has_reader_evidence:bool": any(
                                source in reader_sources for source in sources
                            ),
                            "medscan_only:bool": medscan_only,
                            "sparser_only:bool": sources == {"sparser"},
                        }
                        total_count += 1
                        hashes_yielded.add(stmt_hash)
                        yield Relation(
                            values["agA_ns"],
                            values["agA_id"],
                            values["agB_ns"],
                            values["agB_id"],
                            rel_type,
                            data,
                        )
                    # This statement is not in df
                    except KeyError:
                        continue
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
                sources = set(json.loads(source_counts))
                data = {
                    "stmt_hash:long": stmt_hash,
                    "source_counts:string": source_counts,
                    "evidence_count:int": evidence_count,
                    "stmt_type:string": stmt_type,
                    "belief:float": belief,
                    "has_database_evidence:bool": any(
                        source in db_sources for source in sources
                    ),
                    "has_reader_evidence:bool": any(
                        source in reader_sources for source in sources
                    ),
                    "medscan_only:bool": sources == {"medscan"},
                    "sparser_only:bool": sources == {"sparser"},
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
        logger.info(f"Got {total_count} total relations")

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


class EvidenceProcessor(Processor):
    name = "indra_db_evidence"
    node_types = ["Evidence", "Publication"]

    def __init__(self):
        base_path = pystow.module("indra", "db")
        self.statements_path = base_path.join(name="statements_with_evidences.tsv.gz")
        self.text_refs_path = base_path.join(name="text_refs_for_reading.tsv.gz")
        self.sif_path = base_path.join(name="sif.pkl")
        self._stmt_id_pmid_links = {}
        # Check if files exist without loading them
        for path in [self.statements_path, self.text_refs_path, self.sif_path]:
            if not path.exists():
                raise FileNotFoundError(f"No such file: {path}")

    def get_nodes(self) -> Iterable[Node]:
        """Get nodes from the SIF file."""
        # Load the text ref lookup so that we can set text refs in
        # evidences
        logger.info("Getting text refs from text refs file")
        text_refs = load_text_refs_for_reading_dict(self.text_refs_path.as_posix())
        # Get a list of hashes from the SIF file so that we only
        # add nodes/relations for statements that are in the SIF file
        logger.info("Getting hashes from SIF file")
        with open(self.sif_path, "rb") as fh:
            sif = pickle.load(fh)
        sif_hashes = set(sif["stmt_hash"])
        logger.info("Getting statements from statements file")
        with gzip.open(self.statements_path, "rt", encoding="utf-8") as fh:
            # TODO test whether this is a reasonable size
            batch_size = 100000
            # TODO get number of batches from the total number of statements
            # rather than hardcoding
            total = 352
            reader = csv.reader(fh, delimiter="\t")
            for batch in tqdm(
                batch_iter(reader, batch_size=batch_size, return_func=list),
                total=total,
            ):
                node_batch = []
                for raw_stmt_id, reading_id, stmt_hash, raw_json_str, _ in batch:
                    stmt_hash = int(stmt_hash)
                    if stmt_hash not in sif_hashes:
                        continue
                    try:
                        raw_json = load_statement_json(raw_json_str)
                    except StatementJSONDecodeError as e:
                        logger.warning(e)
                    evidence = raw_json["evidence"][0]
                    # Set text refs and get Publication node
                    pubmed_node = None
                    if reading_id != "\\N":
                        tr = text_refs[reading_id]
                        evidence["text_refs"] = tr
                        if "PMID" in evidence["text_refs"]:
                            evidence["pmid"] = evidence["text_refs"]["PMID"]
                            self._stmt_id_pmid_links[raw_stmt_id] = evidence["pmid"]
                            pubmed_node = Node(
                                "PUBMED",
                                evidence["pmid"],
                                ["Publication"],
                                {
                                    "trid": tr.get("TRID"),
                                    "pmcid": tr.get("PMCID"),
                                    "doi": tr.get("DOI"),
                                    "pii": tr.get("PII"),
                                    "url": tr.get("URL"),
                                    "manuscript_id": tr.get("ManuscriptID"),
                                },
                            )
                        else:
                            evidence["pmid"] = None
                    else:
                        if evidence.get("pmid"):
                            self._stmt_id_pmid_links[raw_stmt_id] = evidence["pmid"]
                            pubmed_node = Node(
                                "PUBMED",
                                evidence["pmid"],
                                ["Publication"],
                            )

                    source_api = evidence["source_api"]
                    source_api = internal_source_mappings.get(source_api, source_api)
                    node_batch.append(
                        Node(
                            "indra_evidence",
                            raw_stmt_id,
                            ["Evidence"],
                            {
                                "evidence:string": json.dumps(evidence),
                                "stmt_hash:long": stmt_hash,
                                "source_api:str": source_api,
                            },
                        )
                    )
                    if pubmed_node:
                        node_batch.append(pubmed_node)
                yield node_batch

    def get_relations(self):
        for stmt_id, pmid in self._stmt_id_pmid_links.items():
            yield Relation("indra_evidence", stmt_id, "PUBMED", pmid, "has_citation")

    def _dump_nodes(self) -> Path:
        # This overrides the default implementation in Processor because
        # we want to process Evidence nodes in batches
        paths_by_type = {}
        nodes_by_type = defaultdict(list)
        # Process Evidence and Publication nodes differently
        evid_node_type = "Evidence"
        pmid_node_type = "Publication"
        nodes_path, nodes_indra_path, sample_path = self._get_node_paths(evid_node_type)
        paths_by_type[evid_node_type] = nodes_path
        # From each batch get the nodes by type but only process Evidence nodes at the moment
        for bidx, nodes in enumerate(self.get_nodes()):
            logger.info(f"Processing batch {bidx}")
            for node in nodes:
                nodes_by_type[node.labels[0]].append(node)
            # We'll append all batches to a single tsv file
            write_mode = "wt"
            if bidx > 0:
                sample_path = None
                write_mode = "at"
            nodes = sorted(
                nodes_by_type[evid_node_type], key=lambda x: (x.db_ns, x.db_id)
            )
            self._dump_nodes_to_path(nodes, nodes_path, sample_path, write_mode)
            # Remove Evidence nodes batch because we don't need to keep them in memory,
            # keep the Publication nodes since we haven't processed them yet
            nodes_by_type[evid_node_type] = []
        # Now process the Publication nodes
        nodes_path, nodes_indra_path, sample_path = self._get_node_paths(pmid_node_type)
        paths_by_type[pmid_node_type] = nodes_path
        with open(nodes_indra_path, "wb") as fh:
            pickle.dump(nodes, fh)
        nodes = sorted(nodes_by_type[pmid_node_type], key=lambda x: (x.db_ns, x.db_id))
        nodes_by_type[pmid_node_type] = nodes
        self._dump_nodes_to_path(nodes, nodes_path, sample_path)
        return paths_by_type, dict(nodes_by_type)

    @classmethod
    def _get_node_paths(cls, node_type: str) -> Path:
        if node_type == "Publication":
            return super()._get_node_paths(node_type)
        elif node_type == "Evidence":
            return (
                cls.module.join(name=f"nodes_{node_type}.tsv.gz"),
                None,
                cls.module.join(name=f"nodes_{node_type}_sample.tsv"),
            )


class StatementJSONDecodeError(Exception):
    pass


def load_statement_json(json_str: str, attempt: int = 1, max_attempts: int = 5) -> json:
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        if attempt < max_attempts:
            json_str = codecs.escape_decode(json_str)[0].decode()
            return load_statement_json(
                json_str, attempt=attempt + 1, max_attempts=max_attempts
            )
    raise StatementJSONDecodeError(
        f"Could not decode statement JSON after " f"{attempt} attempts: {json_str}"
    )


def load_text_refs_for_reading_dict(fname: str):
    text_refs = {}
    for line in tqdm(
        gzip.open(fname, "rt", encoding="utf-8"),
        desc="Processing text refs for readings into a lookup dictionary",
    ):
        ids = line.strip().split("\t")
        id_names = ["TRID", "PMID", "PMCID", "DOI", "PII", "URL", "MANUSCRIPT_ID"]
        d = {}
        rid = ids[0]
        for id_name, id_val in zip(id_names, ids[1:]):
            if id_val != "\\N":
                d[id_name] = id_val
        text_refs[rid] = d
    return text_refs


def ensure_statements_with_evidences(fname):
    if os.path.exists(fname):
        logger.info(f"Found existing statements with evidences in {fname}")
        return
    from indra_db import get_ro

    db = get_ro("primary")
    os.environ["PGPASSWORD"] = db.url.password
    logger.info(f"Dumping statements with evidences into {fname}")
    command = textwrap.dedent(
        f"""
        psql -d {db.url.database} -h {db.url.host} -U {db.url.username}
        -c "COPY (SELECT id, reading_id, mk_hash, encode(raw_json::bytea, 'escape'),
        encode(pa_json::bytea, 'escape') FROM readonly.fast_raw_pa_link) TO STDOUT"
        | gzip > {fname}
    """
    ).replace("\n", " ")
    os.system(command)


def ensure_text_refs_for_reading(fname):
    if os.path.exists(fname):
        logger.info(f"Found existing text refs for reading in {fname}")
        return
    from indra_db import get_ro

    db = get_ro("primary")
    os.environ["PGPASSWORD"] = db.url.password
    logger.info(f"Dumping text refs for reading into {fname}")
    command = textwrap.dedent(
        f"""
        psql -d {db.url.database} -h {db.url.host} -U {db.url.username}
        -c "COPY (SELECT rid, trid, pmid, pmcid, doi, pii, url, manuscript_id 
        FROM readonly.reading_ref_link) TO STDOUT"
        | gzip > {fname}
    """
    ).replace("\n", " ")
    os.system(command)
