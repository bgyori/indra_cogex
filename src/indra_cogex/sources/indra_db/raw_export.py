import logging
import os
from collections import defaultdict
import csv
import gzip
import json
import tqdm
import codecs
import pandas
import pickle
import pystow
from indra.util import batch_iter
from indra.statements import stmts_from_json, stmts_to_json
from indra.tools import assemble_corpus as ac

base_folder = pystow.module("indra", "db")
reading_text_content_fname = base_folder.join(name="reading_text_content_meta.tsv.gz")
text_refs_fname = base_folder.join(name="text_refs_principal.tsv.gz")
raw_stmts_fname = base_folder.join(name="raw_statements.tsv.gz")
drop_readings_fname = base_folder.join(name="drop_readings.pkl")
reading_to_text_ref_map = base_folder.join(name="reading_to_text_ref_map.pkl")

processed_stmts_fname = base_folder.join(name="processed_statements.tsv.gz")
grounded_stmts_fname = base_folder.join(name="grounded_statements.tsv.gz")
unique_stmts_fname = base_folder.join(name="unique_statements.tsv.gz")
source_counts_fname = base_folder.join(name="source_counts.pkl")


logger = logging.getLogger(__name__)


class StatementJSONDecodeError(Exception):
    pass


def load_statement_json(json_str: str, attempt: int = 1, max_attempts: int = 5):
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


def reader_prioritize(reader_contents):
    drop = set()
    # We first organize the contents by source/text type
    versions_per_type = defaultdict(list)
    for reading_id, reader_version, source, text_type in reader_contents:
        versions_per_type[(source, text_type)].append((reader_version, reading_id))
    # Then, within each source/text_type key, we sort according to reader
    # version to be able to select the newest reader version for each key
    reading_id_per_type = {}
    for (source, text_type), versions in versions_per_type.items():
        if len(versions) > 1:
            sorted_versions = sorted(
                versions,
                key=lambda x: reader_versions[version_to_reader[x[0]]].index(x[0]),
                reverse=True,
            )
            drop |= {x[1] for x in sorted_versions[1:]}
            reading_id_per_type[(source, text_type)] = sorted_versions[0][1]
        else:
            reading_id_per_type[(source, text_type)] = versions[0][1]
    fulltexts = [
        content
        for content in reader_contents
        if content[3] == "fulltext" and content[0] not in drop
    ]
    not_fulltexts = [
        content
        for content in reader_contents
        if content[3] != "fulltext" and content[0] not in drop
    ]
    # There are 3 types of non-fulltext content: CORD-19 abstract, PubMed abstract
    # and PubMed title. If we have CORD-19, we prioritize it because it includes
    # the title so we drop any PubMed readings. Otherwise we don't drop anything
    # because we want to keep both the title and the abstract (which doesn't include
    # the title).
    if not fulltexts:
        if ("cord19_abstract", "abstract") in reading_id_per_type:
            if ("pubmed", "abstract") in reading_id_per_type:
                drop.add(reading_id_per_type[("pubmed", "abstract")])
            if ("pubmed", "title") in reading_id_per_type:
                drop.add(reading_id_per_type[("pubmed", "title")])
    # In case of fulltext, we can drop all non-fulltexts, and then drop
    # everything that is lower on the fulltext priority order
    else:
        priority = [
            "xdd-pubmed",
            "xdd",
            "xdd-biorxiv",
            "cord19_pdf",
            "elsevier",
            "cord19_pmc_xml",
            "manuscripts",
            "pmc_oa",
        ]
        drop |= {c[0] for c in not_fulltexts}
        sorted_fulltexts = sorted(
            fulltexts, key=lambda x: priority.index(x[2]), reverse=True
        )
        drop |= {c[0] for c in sorted_fulltexts[1:]}
    return drop


def load_text_refs_by_trid(fname: str):
    text_refs = {}
    for line in tqdm.tqdm(
        gzip.open(fname, "rt", encoding="utf-8"),
        desc="Processing text refs into a lookup dictionary",
    ):
        ids = line.strip().split("\t")
        id_names = ["TRID", "PMID", "PMCID", "DOI", "PII", "URL", "MANUSCRIPT_ID"]
        d = {}
        for id_name, id_val in zip(id_names, ids):
            if id_val != "\\N":
                if id_name == "TRID":
                    id_val = int(id_val)
                d[id_name] = id_val
        text_refs[int(ids[0])] = d
    return text_refs


def get_update(start_date):
    # start date e.g., '20220112 00:00:00'
    query = """
        SELECT tc.insert_date, tr.pmid, tr.pmcid
        FROM public.text_content as tc, public.text_ref as tr
        WHERE tc.text_ref_id = tr.id
        AND tc.insert_date >= start_date
    """


if __name__ == "__main__":
    reader_versions = {
        "sparser": [
            "sept14-linux\n",
            "sept14-linux",
            "June2018-linux",
            "October2018-linux",
            "February2020-linux",
            "April2020-linux",
        ],
        "reach": [
            "61059a-biores-e9ee36",
            "1.3.3-61059a-biores-",
            "1.6.1",
            "1.6.3-e48717",
        ],
        "trips": ["STATIC", "2019Nov14", "2021Jan26"],
        "isi": ["20180503"],
        "eidos": ["0.2.3-SNAPSHOT"],
        "mti": ["1.0"],
    }

    version_to_reader = {}
    for reader, versions in reader_versions.items():
        for version in versions:
            version_to_reader[version] = reader

    # STAGE 0: Dump tsv.gz files directly from the principal database.
    # These can be done from the command line, in the folder
    # that corresponds to pystow.join('indra', 'db')

    command_line = """
    Text refs

      psql -d indradb_test -h indradb-refresh.cwcetxbvbgrf.us-east-1.rds.amazonaws.com
      -U tester -c "COPY (SELECT id, pmid, pmcid, doi, pii, url, manuscript_id
      FROM public.text_ref) TO STDOUT" | gzip > text_refs_principal.tsv.gz

    Time estimate: ~2.5 mins

    Text content joined with reading

      psql -d indradb_test -h indradb-refresh.cwcetxbvbgrf.us-east-1.rds.amazonaws.com
      -U tester -c "COPY (SELECT rd.id, rd.reader_version, tc.id, tc.text_ref_id, 
      tc.source, tc.text_type FROM public.text_content as tc, public.reading as rd
      WHERE tc.id = rd.text_content_id) TO STDOUT" | gzip > reading_text_content_meta.tsv.gz

    Time estimate: ~15 mins

    Raw statements

      psql -d indradb_test -h indradb-refresh.cwcetxbvbgrf.us-east-1.rds.amazonaws.com 
      -U tester -c "COPY (SELECT id, db_info_id, reading_id, encode(json::bytea, 'escape') FROM public.raw_statements) 
      TO STDOUT" | gzip > raw_statements.tsv.gz

    Time estimate: ~30-40 mins
    """

    needed_files = [reading_text_content_fname, text_refs_fname, raw_stmts_fname]
    if any(not f.exists() for f in needed_files):
        missing = [f.as_posix() for f in needed_files if not f.exists()]
        print(command_line)
        raise FileNotFoundError(f"{', '.join(missing)} missing, please run "
                                f"the command(s) above to get them.")

    if not os.environ.get('INDRA_DB_LITE_LOCATION'):
        raise ValueError("Environment variable 'INDRA_DB_LITE_LOCATION' not set")

    # STAGE 1: We need to run statement distillation to figure out which
    # raw statements we should ignore based on the text content and
    # reader version used to produce it.
    if not drop_readings_fname.exists():
        df = pandas.read_csv(
            reading_text_content_fname,
            header=None,
            sep="\t",
            names=[
                "reading_id",
                "reader_version",
                "text_content_id",
                "text_ref_id",
                "text_content_source",
                "text_content_type",
            ],
        )
        df.sort_values("text_ref_id", inplace=True)

        drop_readings = set()
        trid = df["text_ref_id"].iloc[0]
        contents = defaultdict(list)

        # This takes around 1.5 hours
        for row in tqdm.tqdm(df.itertuples(), total=len(df)):
            if row.text_ref_id != trid:
                for reader, reader_contents in contents.items():
                    if len(reader_contents) < 2:
                        continue
                    drop_new = reader_prioritize(reader_contents)
                    # A sanity check to make sure we don't drop all
                    # the readings
                    assert len(drop_new) < len(reader_contents)
                    drop_readings |= drop_new
                contents = defaultdict(list)
            contents[version_to_reader[row.reader_version]].append(
                (
                    row.reading_id,
                    row.reader_version,
                    row.text_content_source,
                    row.text_content_type,
                )
            )
            trid = row.text_ref_id

        with open(drop_readings_fname, "wb") as fh:
            pickle.dump(drop_readings, fh)

        # Dump mapping of reading_id to text_ref_id
        reading_id_to_text_ref_id = dict(zip(df.reading_id, df.text_ref_id))
        with reading_to_text_ref_map.open("wb") as fh:
            pickle.dump(reading_id_to_text_ref_id, fh)

    else:
        with open(drop_readings_fname, "rb") as fh:
            drop_readings = pickle.load(fh)
        # Get mapping of reading_id to text_ref_id
        with reading_to_text_ref_map.open("rb") as fh:
            reading_id_to_text_ref_id = pickle.load(fh)

    text_refs = load_text_refs_by_trid(text_refs_fname.as_posix())

    # STAGE 2: We now need to iterate over raw statements and do preassembly
    if not processed_stmts_fname.exists() or not source_counts_fname.exists():
        logger.info("Preassembling statements and collecting source counts")
        source_counts = defaultdict(lambda: defaultdict(int))
        with gzip.open(raw_stmts_fname, "rt") as fh, gzip.open(
            processed_stmts_fname, "wt"
        ) as fh_out:
            reader = csv.reader(fh, delimiter="\t")
            writer = csv.writer(fh_out, delimiter="\t")
            for lines in tqdm.tqdm(batch_iter(reader, 10000), total=7000):
                stmts_jsons = []
                for raw_stmt_id, db_info_id, reading_id, stmt_json_raw in lines:
                    refs = None
                    if reading_id != "\\N":
                        # Skip if this is for a dropped reading
                        if int(reading_id) in drop_readings:
                            continue
                        text_ref_id = reading_id_to_text_ref_id.get(int(reading_id))
                        if text_ref_id:
                            refs = text_refs.get(text_ref_id)
                    stmt_json = load_statement_json(stmt_json_raw)
                    if refs:
                        stmt_json["evidence"][0]["text_refs"] = refs
                        if refs.get("PMID"):
                            stmt_json["evidence"][0]["pmid"] = refs["PMID"]
                    stmts_jsons.append(stmt_json)
                stmts = stmts_from_json(stmts_jsons)
                stmts = ac.fix_invalidities(stmts, in_place=True)
                stmts = ac.map_grounding(stmts)
                stmts = ac.map_sequence(stmts)
                for stmt in stmts:
                    stmt_hash = stmt.get_hash(refresh=True)
                    source_counts[stmt_hash][stmt.evidence[0].source_api] += 1
                rows = [(stmt.get_hash(), json.dumps(stmt.to_json())) for stmt in stmts]
                writer.writerows(rows)

        # Cast defaultdict to dict and pickle it
        logger.info("Dumping source counts")
        source_counts = dict(source_counts)
        with open(source_counts_fname.as_posix(), "wb") as fh:
            pickle.dump(source_counts, fh)
        # We can now remove reference to source counts
        del source_counts

    else:
        logger.info(
            f"Statements already preassembled at "
            f"{processed_stmts_fname.as_posix()}, source counts "
            f"already dumped to {source_counts_fname.as_posix()}, skipping..."
        )

    # STAGE 3: create grounded and unique dumps
    with gzip.open(processed_stmts_fname, "rt") as fh, gzip.open(
        grounded_stmts_fname, "wt"
    ) as fh_out_gr, gzip.open(unique_stmts_fname, "wt") as fh_out_uniq:
        seen_hashes = set()
        reader = csv.reader(fh, delimiter="\t")
        writer_gr = csv.writer(fh_out_gr, delimiter="\t")
        writer_uniq = csv.writer(fh_out_uniq, delimiter="\t")
        for sh, stmt_json_str in tqdm.tqdm(reader, total=65102088):
            stmt = stmts_from_json([load_statement_json(stmt_json_str)])[0]
            stmt = load_statement_json(stmt_json_str)
            if len(stmt.real_agent_list()) < 2:
                continue
            if all(
                (set(agent.db_refs) - {"TEXT", "TEXT_NORM"})
                for agent in stmt.real_agent_list()
            ):
                writer_gr.writerow((sh, stmt_json_str))
                if sh not in seen_hashes:
                    writer_uniq.writerow((sh, stmt_json_str))
            seen_hashes.add(sh)


"""Notes
Dumped files into Ubuntu ~/.data/indra/db/

I forgot to include the reader in the readings metadata but that can be
reversed out using the inverse indra_db.databases.reader_versions dict so
we don't need to fix that necessarily.

Summary of what is in the text content joined with reading file

('1.3.3-61059a-biores-', 58572394),
 ('October2018-linux', 51637830),
 ('sept14-linux', 19142899),
 ('June2018-linux', 10469377),
 ('April2020-linux', 7693672),
 ('0.2.3-SNAPSHOT', 6661714),
 ('20180503', 4647511),
 ('1.6.3-e48717', 3809726),
 ('2019Nov14', 1667674),
 ('1.0', 1017773),
 ('2021Jan26', 773324),
 ('1.6.1', 238176),
 ('February2020-linux', 61676),
 ('STATIC', 48053)

(('pubmed', 'abstract'), 75625251),
 (('pubmed', 'title'), 74614345),
 (('pmc_oa', 'fulltext'), 10286339),
 (('manuscripts', 'fulltext'), 2128472),
 (('elsevier', 'fulltext'), 1396977),
 (('cord19_abstract', 'abstract'), 1000864),
 (('cord19_pdf', 'fulltext'), 649855),
 (('cord19_pmc_xml', 'fulltext'), 551702),
 (('xdd-pubmed', 'fulltext'), 174212),
 (('xdd', 'fulltext'), 7313),
 (('xdd-biorxiv', 'fulltext'), 6469)

"""
