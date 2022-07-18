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
from itertools import combinations
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union, List, Dict

import click
import humanize
import pandas as pd
import pystow
from indra.databases.identifiers import ensure_prefix_if_needed
from indra.ontology.bio import bio_ontology
from indra.statements import Agent
from indra.util import batch_iter
from indra.util.statement_presentation import db_sources, reader_sources
from more_click import verbose_option
from tqdm import tqdm

from indra_cogex.representation import Node, Relation
from indra_cogex.sources.processor import Processor

from indra_cogex.sources.indra_db.assembly import belief_scores_pkl_fname
from indra_cogex.sources.indra_db.raw_export import (
    source_counts_fname,
    unique_stmts_fname,
    stmts_from_json,
)

logger = logging.getLogger(__name__)
tqdm.pandas()


# If you don't have the data, run the script in raw_export.py and then in
# assembly.py (both in this directory) to get it.


class DbProcessor(Processor):
    """Processor for the INDRA database."""

    name = "database"
    node_types = ["BioEntity"]

    def __init__(self, dir_path: Union[None, str, Path] = None):
        """Initialize the INDRA database processor.

        Parameters
        ----------
        dir_path :
            The path to the directory containing unique and grounded
            statements as a *.tsv.gz file, source counts as a pickle file and
            belief scores as a pickle file.
        """
        if dir_path is None:
            dir_path = unique_stmts_fname.parent
        elif isinstance(dir_path, str):
            dir_path = Path(dir_path)
        self.stmts_fname = dir_path / unique_stmts_fname.name
        self.source_counts_fname = dir_path / source_counts_fname.name
        self.belief_scores_fname = dir_path / belief_scores_pkl_fname.name

    def get_nodes(self):  # noqa:D102
        # Read the unique statements from the file and yield unique agents
        # The file contains statements that have already been filtered for
        # ungrounded statements, so we can just use the agent list.
        batch_size = 100000
        with gzip.open(self.stmts_fname, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            seen_agents = set()  # Store ns:id pairs of seen agents

            for batch in batch_iter(reader, batch_size=batch_size, return_func=list):
                sj_list = [load_statement_json(sjs) for _, sjs in batch]
                stmts = stmts_from_json(sj_list)
                for stmt in stmts:
                    for agent in stmt.real_agent_list():
                        db_ns, db_id = agent.get_grounding()
                        if (db_ns, db_id) not in seen_agents:
                            # Todo: do we need to use bio_ontology.get_name()?
                            yield Node(
                                db_ns, db_id, ["BioEntity"], dict(name=agent.name)
                            )
                            seen_agents.add((db_ns, db_id))

    def get_relations(self, max_complex_members: int = 3):  # noqa:D102
        # todo: Should this method call the source scripts for statements,
        #  source counts and belief calculations? They will take 12+ hours
        #  altogether.
        rel_type = "indra_rel"
        total_count = 0

        # Load the source counts and belief scores into dictionaries
        logger.info("Loading source counts per hash")
        with self.source_counts_fname.open("rb") as f:
            source_counts = pickle.load(f)
        logger.info("Loading belief scores per hash")
        with self.belief_scores_fname.open("rb") as f:
            belief_scores = pickle.load(f)

        hashes_yielded = set()
        with gzip.open(self.stmts_fname, "rt") as fh:
            reader = csv.reader(fh, delimiter="\t")
            for sh_str, stmt_json_str in tqdm(reader, desc="Reading statements"):
                stmt_hash = int(sh_str)
                try:
                    source_count = source_counts[stmt_hash]
                    belief = belief_scores[stmt_hash]
                except KeyError:
                    # NOTE: this should not happen if files are generated
                    # properly and are up to date.
                    logger.warning(
                        f"Could not find source count or belief score for "
                        f"statement hash {stmt_hash}. Are the source files updated?"
                    )
                    continue
                # If we already yielded this statement, we can skip it
                # NOTE: statements should be unique.
                if stmt_hash in hashes_yielded:
                    logger.warning(
                        f"Found duplicate hash {stmt_hash} among unique statements"
                    )
                    continue
                stmt_json = load_statement_json(stmt_json_str)
                ev_list = stmt_json["evidence"]
                stripped_json = [ev_list[0]]
                data = {
                    "stmt_hash:long": stmt_hash,
                    "source_counts:string": json.dumps(source_count),
                    "evidence_count:int": len(ev_list),
                    "stmt_type:string": stmt_json["stmt_type"],
                    "belief:float": belief,
                    "stmt_json:string": json.dumps(stripped_json),
                    "has_database_evidence:bool": any(
                        source in db_sources for source in source_counts
                    ),
                    "has_reader_evidence:bool": any(
                        source in reader_sources for source in source_counts
                    ),
                    "medscan_only:bool": set(source_counts) == {"medscan"},
                    "sparser_only:bool": set(source_counts) == {"sparser"},
                }

                # Get the agents from the statement
                agents: List[Agent] = [
                    ag for ag in stmt_json["agent_list"] if ag is not None
                ]
                if len(agents) < 2:
                    logger.warning(
                        f"Statement {stmt_hash} has less than two agents. Skipping."
                    )
                    continue
                elif len(agents) == 2:
                    ag_a, ag_b = agents
                    ns_a, id_a = ag_a.get_grounding()
                    ns_b, id_b = ag_b.get_grounding()
                    yield Relation(
                        ns_a,
                        id_a,
                        ns_b,
                        id_b,
                        rel_type,
                        data,
                    )
                    total_count += 1
                elif 2 < len(agents) <= max_complex_members:
                    # This is a complex, so we need to yield a relation for
                    # each pair of agents of the complex
                    for ag_a, ag_b in combinations(agents, 2):
                        ns_a, id_a = ag_a.get_grounding()
                        ns_b, id_b = ag_b.get_grounding()
                        yield Relation(
                            ns_a,
                            id_a,
                            ns_b,
                            id_b,
                            rel_type,
                            data,
                        )
                        total_count += 1
                else:
                    logger.warning(
                        f"Statement {stmt_hash} has more than {max_complex_members} agents. Skipping."
                    )
                    continue

                hashes_yielded.add(stmt_hash)

        logger.info(
            f"Got {total_count} total relations from {len(hashes_yielded)} unique statements"
        )

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
                    node_batch.append(
                        Node(
                            "indra_evidence",
                            raw_stmt_id,
                            ["Evidence"],
                            {
                                "evidence:string": json.dumps(evidence),
                                "stmt_hash:long": stmt_hash,
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


def load_raw_stmt_id_for_reading_id_dict(fname: str) -> Dict[str, str]:
    reading_id_for_raw_stmt_id = {}
    with gzip.open(fname, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for raw_stmt_id, db_info_id, reading_id, stmt_json_raw in reader:

            reading_id_for_raw_stmt_id[reading_id] = raw_stmt_id
    return reading_id_for_raw_stmt_id


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
