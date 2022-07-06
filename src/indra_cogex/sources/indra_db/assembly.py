import math
import json
import tqdm
import codecs
import pystow
import sqlite3
from collections import defaultdict
from indra.ontology.bio import bio_ontology
from indra.statements.io import stmt_from_json
from indra.preassembler import Preassembler

base_folder = pystow.module("indra", "db")


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


def get_related(stmts):
    stmts_by_type = defaultdict(list)
    for stmt in stmts:
        stmts_by_type[stmt.__class__.__name__].append(stmt)
    refinements = set()
    for _, stmts_this_type in stmts_by_type.items():
        refinements |= pa._generate_relation_tuples(stmts_this_type)
    return refinements


def get_related_split(stmts1, stmts2):
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


if __name__ == "__main__":
    db = sqlite3.connect(base_folder.join(name="statements.db"))
    bio_ontology.initialize()
    bio_ontology._build_transitive_closure()
    pa = Preassembler(bio_ontology)

    cur = db.execute("select count(1) from processed")
    num_rows = cur.fetchone()[0]
    batch_size = 1e6

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


"""
Assembly notes:

Step 1: Create a SQLITE DB

sqlite3 -batch statements.db "create table processed (hash integer, stmt text);"
zcat < unique_statements.tsv.gz | sqlite3 -cmd ".mode tabs" -batch statements.db ".import '|cat -' processed"
sqlite3 -batch statements.db "create index processed_idx on processed (hash);"
"""

"""
Other possible approaches

Open two CSV readers for the unique_statements.tsv.gz and then move them
forward in batches to cover all combinations of Statement batches.

- First batch of Stmts, internal refinement finding
- First batch (first reader) x Second batch (second reader)
- First batch (first reader) x Third batch (second reader)
- ...
- One before last batch (first reader) x Last batch (second reader)
---> Giant list of refinement relation pairs (hash1, hash2)
"""

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
