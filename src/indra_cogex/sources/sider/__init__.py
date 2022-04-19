# -*- coding: utf-8 -*-

"""Processor for SIDER."""

import re
from collections import Counter
from typing import Iterable

import gilda
import gilda.grounder
import pandas as pd
import pystow
from indra.databases import biolookup_client
from indra.databases.identifiers import get_ns_id_from_identifiers
from tabulate import tabulate
from tqdm import tqdm

from indra_cogex.representation import Node, Relation
from indra_cogex.sources import Processor
from indra_cogex.sources.utils import UmlsMapper

VERSION = "4.1"
SUBMODULE = pystow.module("indra", "cogex", "sider", VERSION)
URL = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"
SIDE_EFFECTS_HEADER = [
    # 'STITCH_FLAT_ID',
    "STITCH_STEREO_ID",
    # 'UMLS CUI from Label',
    "MedDRA Concept Type",
    "UMLS CUI from MedDRA",
    "MedDRA Concept name",
]


cid_to_pubchem_pattern = re.compile(r"^CID(?:0)+(\d+)$")


def stitch_stereo_to_pubchem(cid: str) -> str:
    assert cid.startswith("CID")
    return re.sub(cid_to_pubchem_pattern, "\\1", cid)


class SIDERSideEffectProcessor(Processor):
    """A processor for SIDER side effects."""

    name = "sider_side_effects"
    node_types = ["BioEntity"]

    def __init__(self):
        self.df = SUBMODULE.ensure_csv(
            url=URL,
            read_csv_kwargs=dict(
                usecols=[1, 3, 4, 5],
                names=SIDE_EFFECTS_HEADER,
                dtype=str,
            ),
        )
        self.df = self.df[self.df["STITCH_STEREO_ID"].notna()]
        self.df = self.df[self.df["UMLS CUI from MedDRA"].notna()]

        # Prepare chemicals
        self.df["pubchem_id"] = self.df["STITCH_STEREO_ID"].map(
            stitch_stereo_to_pubchem
        )
        del self.df["STITCH_STEREO_ID"]
        self.chemicals = {
            pubchem_id: Node.standardized(
                db_ns="PUBCHEM",
                db_id=pubchem_id,
                labels=["BioEntity"],
            )
            for pubchem_id in tqdm(
                self.df["pubchem_id"], unit_scale=True, desc="Caching chemicals"
            )
        }
        for node in tqdm(self.chemicals.values(), desc="Finding chemical names"):
            if node.data["name"] is None:
                node.data["name"] = biolookup_client.get_name(node.db_ns, node.db_id)

        umls_mapper = UmlsMapper()
        self.side_effects = {}
        for umls_id in self.df["UMLS CUI from MedDRA"].unique():
            # TODO replace with "standardize"
            prefix, identifier, name = umls_mapper.lookup(umls_id)
            db_ns, db_id = get_ns_id_from_identifiers(prefix, identifier)
            if db_ns is None:
                db_ns, db_id = prefix, identifier
            self.side_effects[umls_id] = Node.standardized(
                db_ns=db_ns, db_id=db_id, name=name, labels=["BioEntity"]
            )

    def get_nodes(self) -> Iterable[Node]:
        """Iterate over SIDER chemicals and side effects."""
        yield from self.chemicals.values()
        yield from self.side_effects.values()

    def get_relations(self) -> Iterable[Relation]:
        """Iterate over SIDER side effect annotations."""
        for pubchem_id, umls_id in self.df[
            ["pubchem_id", "UMLS CUI from MedDRA"]
        ].values:
            chemical = self.chemicals[pubchem_id]
            indication = self.side_effects[umls_id]
            yield Relation(
                chemical.db_ns,
                chemical.db_id,
                indication.db_ns,
                indication.db_id,
                "has_side_effect",
                dict(
                    source=self.name,
                    version=VERSION,
                ),
            )


def generate_curation_sheet():
    """This function does some statistics over the UMLS mappings out of SIDER and generates
    a curation sheet based on gilda when possible. Still, about half of the terms are only
    available in UMLS at the moment.
    """
    umls_mapper = UmlsMapper()

    df = SUBMODULE.ensure_csv(
        url=URL,
        read_csv_kwargs=dict(
            usecols=[1, 3, 4, 5],
            names=SIDE_EFFECTS_HEADER,
            dtype=str,
        ),
    )
    df["pubchem_id"] = df["STITCH_STEREO_ID"].map(stitch_stereo_to_pubchem)
    del df["STITCH_STEREO_ID"]

    print(tabulate(df.head(), headers=df.columns))
    umls_ids = df["UMLS CUI from MedDRA"].unique()
    unique_umls_id_count = len(umls_ids)
    print("Unique UMLS identifiers:", unique_umls_id_count)
    c = Counter()
    unnormalized = set()
    for umls_id in tqdm(umls_ids):
        db_ns, db_id, name = umls_mapper.lookup(umls_id)
        if db_ns is not None:
            c[db_ns] += 1
        else:
            unnormalized.add(umls_id)

    print("Most common prefixes mapped from UMLS to")
    print(tabulate(c.most_common()))

    number_normalized = sum(c.values())
    print(
        f"normalized  : {number_normalized}/{unique_umls_id_count} ({number_normalized / unique_umls_id_count:.2%})"
    )
    print(
        f"unnormalized: {len(unnormalized)}/{unique_umls_id_count} ({len(unnormalized) / unique_umls_id_count:.2%})"
    )

    print("mapping all rows")
    df["prefix"], df["identifier"], df["name"] = zip(
        *df["UMLS CUI from MedDRA"].map(umls_mapper.lookup)
    )
    df.to_csv(SUBMODULE.join(name="norm_df.tsv"), sep="\t", index=False)

    df_unnormalized = df.loc[df["prefix"].isna()]
    unique_unnormalized = {
        (umls_id, name)
        for umls_id, name in df_unnormalized[
            ["UMLS CUI from MedDRA", "MedDRA Concept name"]
        ].values
        if pd.notna(name) and pd.notna(umls_id)
    }
    rows = []
    for umls_id, name in sorted(unique_unnormalized):
        for result in gilda.ground(name):
            result: gilda.grounder.ScoredMatch
            term: gilda.grounder.Term = result.term
            rows.append(
                (
                    umls_id,
                    name,
                    term.db,
                    term.id,
                    term.entry_name,
                    result.score,
                )
            )
    curation_df = pd.DataFrame(
        rows, columns=["umls_id", "umls_name", "prefix", "identifier", "name", "score"]
    )
    curation_df.to_csv(SUBMODULE.join(name="curation_df.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    SIDERSideEffectProcessor.cli()
