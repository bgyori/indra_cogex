# -*- coding: utf-8 -*-

"""This script gets DisGeNet, converts to BEL, exports in a few formats."""

import logging

import pandas as pd
import pystow
from indra.databases import hgnc_client

from indra_cogex.representation import Node, Relation
from indra_cogex.sources import Processor
from indra_cogex.sources.utils import UmlsMapper

logger = logging.getLogger(__name__)

SUBMODULE = pystow.module("indra", "cogex", "disgenet")

DOWNLOAD_BASE = "https://www.disgenet.org/static/disgenet_ap1/files/downloads"
CURATED_DISEASE_GENES_ASSOCIATIONS_URL = (
    f"{DOWNLOAD_BASE}/curated_gene_disease_associations.tsv.gz"
)
# CURATED_DISEASE_VARIANT_ASSOCIATIONS_URL = f"{DOWNLOAD_BASE}/curated_variant_disease_associations.tsv.gz"


TARGET_KEYS = {
    "NofSnps": int,
    "DSI": str,
    "DPI": str,
    "diseaseType": str,
    "diseaseSemanticType": str,
    "score": str,
    "EI": str,
    "YearInitial": int,
    "YearFinal": int,
    "NofPmids": str,
    "source": lambda s: set(s.split(";")) if s and pd.notna(s) else s,
    "GD_dis1GeneIdSet": int,
    "GD_dis2GeneIdSet": int,
    "GD_commonGeneIdSet": int,
    "GD_jaccard": float,
    "GD_pvalue": float,
    "VD_variant1GeneIdSet": int,
    "VD_variant2GeneIdSet": int,
    "VD_commonVariantIdSet": int,
    "VD_jaccard": float,
    "VD_pvalue": float,
}


class DisgenetProcessor(Processor):
    """Processor for the DisGeNet database."""

    name = "disgenet"
    df: pd.DataFrame
    node_types = ["BioEntity"]
    relation = "gene_disease_association"

    def __init__(self):
        """Initialize the DisGeNet processor."""
        self.df = load_disgenet(CURATED_DISEASE_GENES_ASSOCIATIONS_URL)

    def get_nodes(self):  # noqa:D102
        diseases = {
            tuple(row)
            for row in self.df[["disease_prefix", "disease_id", "disease_name"]].values
        }
        for prefix, identifier, name in diseases:
            yield Node.standardized(
                db_ns=prefix, db_id=identifier, name=name, labels=["BioEntity"]
            )
        for hgnc_id in self.df["hgnc_id"].unique():
            yield Node("HGNC", hgnc_id, ["BioEntity"])

    def get_relations(self):  # noqa:D102
        columns = [
            "hgnc_id",
            "disease_prefix",
            "disease_id",
            "NofSnps",
            "NofPmids",
        ]
        for hgnc_id, disease_prefix, disease_id, snps, papers in self.df[
            columns
        ].values:
            data = {"snps:int": snps, "source": self.name, "papers:int": papers}
            yield Relation(
                "HGNC", hgnc_id, disease_prefix, disease_id, self.relation, data
            )


def load_disgenet(url, force: bool = False) -> pd.DataFrame:
    """Export disease-gene association file."""
    df = SUBMODULE.ensure_csv(
        url=url,
        read_csv_kwargs=dict(dtype={"geneId": str}),
        force=force,
    )

    df["hgnc_id"] = df["geneId"].map(
        lambda s: hgnc_client.get_hgnc_from_entrez(s.strip())
    )
    df = df[df["hgnc_id"].notna()]

    umls_mapper = UmlsMapper()
    (
        df["disease_prefix"],
        df["disease_id"],
        df["disease_name"],
    ) = zip(*df["diseaseId"].map(umls_mapper.standardize))
    df = df[df["disease_prefix"].notna()]
    return df
