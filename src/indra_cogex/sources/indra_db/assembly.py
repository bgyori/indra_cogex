import csv
import gzip
import math
import json
from typing import List, Set, Tuple

import networkx as nx
import tqdm
import codecs
import pystow
import sqlite3
import subprocess
from collections import defaultdict
from indra.ontology.bio import bio_ontology
from indra.statements import Statement
from indra.statements.io import stmt_from_json
from indra.preassembler import Preassembler

from .raw_export import unique_stmts_fname

StmtList = List[Statement]

base_folder = pystow.module("indra", "db")
refinements_fname = base_folder.join(name="refinements.tsv.gz")


class StatementJSONDecodeError(Exception):
    pass


def load_statement_json(
    json_str: str,
    attempt: int = 1,
    max_attempts: int = 5,
    remove_evidence: bool = False,
):
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        if attempt < max_attempts:
            json_str = codecs.escape_decode(json_str)[0].decode()
            sj = load_statement_json(
                json_str, attempt=attempt + 1, max_attempts=max_attempts
            )
            if remove_evidence:
                sj["evidence"] = []
            return sj
    raise StatementJSONDecodeError(
        f"Could not decode statement JSON after " f"{attempt} attempts: {json_str}"
    )


def get_stmts(db, limit, offset):
    cur = db.execute("select * from processed limit %s offset %s" % (limit, offset))
    stmts = [
        stmt_from_json(load_statement_json(sjs, remove_evidence=True))
        for _, sjs in tqdm.tqdm(cur.fetchall(), total=limit, desc="Loading statements")
    ]
    return stmts


def get_related(stmts: StmtList):
    stmts_by_type = defaultdict(list)
    for stmt in stmts:
        stmts_by_type[stmt.__class__.__name__].append(stmt)
    refinements = set()
    for _, stmts_this_type in stmts_by_type.items():
        refinements |= pa._generate_relation_tuples(stmts_this_type)
    return refinements


def get_related_split(stmts1: StmtList, stmts2: StmtList):
    stmts_by_type1 = defaultdict(list)
    stmts_by_type2 = defaultdict(list)
    for stmt in stmts1:
        stmts_by_type1[stmt.__class__.__name__].append(stmt)
    for stmt in stmts2:
        stmts_by_type2[stmt.__class__.__name__].append(stmt)
    refinements = set()
    for stmt_type, stmts_this_type1 in stmts_by_type1.items():
        stmts_this_type2 = stmts_by_type2.get(stmt_type)
        if not stmts_this_type2:
            continue
        refinements |= pa._generate_relation_tuples(
            stmts_this_type1 + stmts_this_type2, split_idx=len(stmts_this_type1) - 1
        )
    return refinements


def sqlite_approach():
    """
    Assembly notes:

    Step 1: Create a SQLITE DB

    sqlite3 -batch statements.db "create table processed (hash integer, stmt text);"
    zcat < unique_statements.tsv.gz | sqlite3 -cmd ".mode tabs" -batch statements.db ".import '|cat -' processed"
    sqlite3 -batch statements.db "create index processed_idx on processed (hash);"
    """
    db = sqlite3.connect(base_folder.join(name="statements.db"))

    cur = db.execute("select count(1) from processed")
    num_rows = cur.fetchone()[0]

    offset0 = 0
    num_batches = math.ceil(num_rows / batch_size)
    refinements = set()
    for i in tqdm.tqdm(range(num_batches)):
        offset1 = i * batch_size
        stmts1 = get_stmts(db, batch_size, offset1)
        refinements |= get_related(stmts1)
        for j in tqdm.tqdm(range(i + 1, num_batches)):
            offset2 = j * batch_size
            stmts2 = get_stmts(db, batch_size, offset2)
            refinements |= get_related_split(stmts1, stmts2)


def belief_calc(refinements_set: Set[Tuple[int, int]]):
    # Belief Calculations
    """
    Belief calculation idea:
    - Load refinements from dump -> we should probably make a networkx graph
      - that will allow us to call e.g., "descendants" or "ancestors" to get
        hashes of statements that refine a given statement
    - Load source counts from dump
    - Iterate over each unique statement
      - Find all the statements that refine it
      - Add up all the source counts for the statement itself and the statements
        that refine it
      - Pass mock evidences based on source counts (or is there an API directly
        for source counts??) to INDRA's Belief Engine to get a belief score
    - Dump dict of stmt_hash: belief score values
    """

    # Make a networkx graph of the refinements
    dg = nx.DiGraph()

    # Edges go from refiner to refined, i.e.
    dg.add_edges_from(refinements_set)

    #


if __name__ == "__main__":
    # Global variables
    bio_ontology.initialize()
    bio_ontology._build_transitive_closure()
    pa = Preassembler(bio_ontology)
    batch_size = int(1e6)

    """
    Step 1, alternative:

    Open two CSV readers for the unique_statements.tsv.gz and then move them
    forward in batches to cover all combinations of Statement batches.

    - First batch of Stmts, internal refinement finding
    - First batch (first reader) x Second batch (second reader)
    - First batch (first reader) x Third batch (second reader)
    - ...
    - One before last batch (first reader) x Last batch (second reader)
    ---> Giant list of refinement relation pairs (hash1, hash2)
    """

    # Count lines in the file
    num_rows = int(subprocess.check_output(
        ["zcat", unique_stmts_fname.as_posix(), "|", "wc", "-l"]
    ).split()[0])
    num_batches = math.ceil(num_rows / batch_size)
    refinements = set()

    # Loop statements: the outer index runs all batches while the inner index
    # runs outer index < inner index <= num_batches. This way the outer
    # index runs the "diagonal" of the combinations while the inner index runs
    # the "off-diagonal" of the combinations.

    # Open two csv readers to the same file
    if not refinements_fname.exists():
        with gzip.open(unique_stmts_fname, "rt") as fh1:
            reader1 = csv.reader(fh1, delimiter="\t")
            for outer_batch_ix in tqdm.tqdm(range(num_batches), total=num_batches,
                                            desc="Processing refinements"):

                # read in a batch from the first reader
                stmts1 = []
                for _ in range(batch_size):
                    try:
                        _, sjs = next(reader1)
                        stmt = stmt_from_json(
                            load_statement_json(sjs, remove_evidence=True)
                        )
                        stmts1.append(stmt)
                    except StopIteration:
                        break

                refinements |= get_related(stmts1)

                # Loop batches from second reader, starting at the
                with gzip.open(unique_stmts_fname, "rt") as fh2:
                    reader2 = csv.reader(fh2, delimiter="\t")
                    for inner_batch_ix in range(outer_batch_ix + 1, num_batches):
                        stmts2 = []
                        for _ in range(batch_size):
                            try:
                                _, sjs = next(reader2)
                                stmt = stmt_from_json(
                                    load_statement_json(sjs, remove_evidence=True)
                                )
                                stmts2.append(stmt)
                            except StopIteration:
                                break

                        refinements |= get_related_split(stmts1, stmts2)

        # Write out the refinements as a gzipped TSV file
        with gzip.open(refinements_fname.as_posix(), "wt") as f:
            tsv_writer = csv.writer(f, delimiter="\t")
            tsv_writer.writerows(refinements)
    else:
        with gzip.open(refinements_fname.as_posix(), "rt") as f:
            tsv_reader = csv.reader(f, delimiter="\t")
            refinements = set(tsv_reader)

    # Step 2: Calculate belief scores
    belief_calc(refinements)
