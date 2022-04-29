# -*- coding: utf-8 -*-

"""Processor for the Human Phenotype Ontology Associations (HPOA) database."""

import logging
from typing import Optional

import bioregistry
import pandas as pd
import pyobo
from indra.databases import hgnc_client
from indra.ontology.standardize import standardize_name_db_refs
from indra.statements.agent import get_grounding

from ..processor import Processor
from ...representation import Node, Relation

__all__ = [
    "HpDiseasePhenotypeProcessor",
    "HpPhenotypeGeneProcessor",
]

logger = logging.getLogger(__name__)


# DISEASE_GENE_URL = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"


class HpDiseasePhenotypeProcessor(Processor):
    """Processor for the Human Phenotype Ontology Annotations (HPOA) database."""

    name = "hp_disease_phenotype"
    df: pd.DataFrame
    node_types = ["BioEntity"]
    rel_type = "has_phenotype"

    def __init__(
        self, df: Optional[pd.DataFrame] = None, version: Optional[str] = None
    ):
        """Initialize the HPOA processor."""
        self.df = get_hpoa_df(version=version) if df is None else df
        self.version = version

    def get_nodes(self):  # noqa:D102
        for prefix, identifier, name in df_unique(
            self.df[["disease_prefix", "disease_id", "disease_name"]]
        ):
            yield Node.standardized(
                db_ns=prefix, db_id=identifier, name=name, labels=["BioEntity"]
            )
        for prefix, identifier, name in df_unique(
            self.df[["phenotype_prefix", "phenotype_id", "phenotype_name"]]
        ):
            yield Node.standardized(
                db_ns=prefix, db_id=identifier, name=name, labels=["BioEntity"]
            )

    def get_relations(self):  # noqa:D102
        for (
            disease_prefix,
            disease_id,
            phenotype_prefix,
            phenotype_id,
        ), evidence in self.df.groupby(
            ["disease_prefix", "disease_id", "phenotype_prefix", "phenotype_id"]
        )[
            "evidence"
        ]:
            all_ecs = ",".join(sorted(set(evidence)))
            # Possible properties could be e.g., evidence codes
            data = {"evidence_codes:string": all_ecs, "source": self.name}
            if self.version:
                data["version"] = self.version
            yield Relation(
                disease_prefix,
                disease_id,
                phenotype_prefix,
                phenotype_id,
                self.rel_type,
                data,
            )


def df_unique(df: pd.DataFrame):
    """Get unique rows in the dataframe as sorted tuples."""
    return sorted({tuple(row) for row in df.values})


def get_hpoa_df(version: Optional[str] = None) -> pd.DataFrame:
    """Get the HPOA database as a preprocessed dataframe.

    Parameters
    ----------
    version :
        An optional version in form of a string YYYY-MM-DD. If not given,
        gets the latest data.
    Returns
    -------
    :
        A parsed and preprocessed HPO annotation dataframe
    """
    if version is None:
        url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
    else:
        url = (
            f"http://purl.obolibrary.org/obo/hp/releases/hpoa/{version}/phenotype.hpoa"
        )
    df = pd.read_csv(
        url,
        sep="\t",
        skiprows=5,
        usecols=[0, 3, 4, 5],
        names=["disease", "phenotype", "reference", "evidence"],
    )
    df.drop_duplicates(inplace=True)
    return process_hpoa_df(df)


def process_hpoa_df(df: pd.DataFrame) -> pd.DataFrame:
    """Process the HPOA dataframe, in place."""
    # 1. Get several mappings not currently available in INDRA's bioontology
    omim_to_mondo = pyobo.get_filtered_xrefs("mondo", "omim", flip=True)
    mondo_to_doid = pyobo.get_filtered_xrefs("mondo", "doid")
    orphanet_to_mondo = pyobo.get_filtered_xrefs("mondo", "orphanet", flip=True)
    omim_to_doid = pyobo.get_filtered_xrefs("doid", "omim", flip=True)

    # Prepare disease mappings from CURIE (either OMIM, Orphanet,
    # or DECIPHER) to DOID and MONDO
    disease_standards = {}
    for disease in df.disease:
        prefix, lui = bioregistry.parse_curie(disease)
        db_xrefs = {
            prefix: lui,
        }
        mondo_id = None
        if prefix == "omim":
            do_id = omim_to_doid.get(lui)
            if do_id:
                db_xrefs["DOID"] = f"DOID:{do_id}"
            mondo_id = omim_to_mondo.get(lui)
            if mondo_id:
                db_xrefs["MONDO"] = mondo_id
        if prefix == "orphanet":
            mondo_id = orphanet_to_mondo.get(lui)
            if mondo_id:
                db_xrefs["MONDO"] = mondo_id
        if mondo_id:
            do_id = mondo_to_doid.get(mondo_id)
            if do_id:
                db_xrefs["DOID"] = f"DOID:{do_id}"

        new_name, new_db_xrefs = standardize_name_db_refs(db_xrefs)
        standard_db, standard_id = get_grounding(new_db_xrefs)
        if not new_name:
            continue
        disease_standards[disease] = standard_db, standard_id, new_name

    (
        df["disease_prefix"],
        df["disease_id"],
        df["disease_name"],
    ) = zip(*df["disease"].map(lambda s: disease_standards.get(s, (None, None, None))))
    del df["disease"]
    # Remove unmappable entries (e.g., DECIPHER entries)
    df = df[df.disease_prefix.notna()]

    df = process_phenotypes(df)
    return df


def process_phenotypes(df: pd.DataFrame, column: str = "phenotype") -> pd.DataFrame:
    """Process the phenotype-gene dataframe, in place."""
    phenotype_standards = {}
    for hp_id in df.phenotype.unique():
        new_name, new_db_xrefs = standardize_name_db_refs({"HP": hp_id})
        standard_db, standard_id = get_grounding(new_db_xrefs)
        if not new_name:
            continue
        phenotype_standards[hp_id] = standard_db, standard_id, new_name

    (
        df["phenotype_prefix"],
        df["phenotype_id"],
        df["phenotype_name"],
    ) = zip(*df[column].map(lambda s: phenotype_standards.get(s, (None, None, None))))
    del df[column]
    df = df[df.phenotype_prefix.notna()]
    return df


class HpPhenotypeGeneProcessor(Processor):
    """"""

    name = "hp_phenotype_gene"
    df: pd.DataFrame
    node_types = ["BioEntity"]
    rel_type = "phenotype_has_gene"

    def __init__(
        self, df: Optional[pd.DataFrame] = None, version: Optional[str] = None
    ):
        """Initialize the HPOA processor."""
        self.df = get_phenotype_gene_df(version=version) if df is None else df
        self.version = version

    def get_nodes(self):  # noqa:D102
        for prefix, identifier, name in df_unique(
            self.df[["phenotype_prefix", "phenotype_id", "phenotype_name"]]
        ):
            yield Node.standardized(
                db_ns=prefix, db_id=identifier, name=name, labels=["BioEntity"]
            )
        for hgnc_id in self.df.hgnc_id.unique():
            yield Node.standardized(db_ns="hgnc", db_id=hgnc_id, labels=["BioEntity"])

    def get_relations(self):  # noqa:D102
        for (phenotype_prefix, phenotype_id, hgnc_id,), evidences in self.df.groupby(
            ["phenotype_prefix", "phenotype_id", "hgnc_id"]
        )["source"]:
            evidence = ",".join(sorted(set(evidences)))
            # Possible properties could be e.g., evidence codes
            data = {"evidence_codes:string": evidence, "source": self.name}
            if self.version:
                data["version"] = self.version
            yield Relation(
                phenotype_prefix,
                phenotype_id,
                "HGNC",
                hgnc_id,
                self.rel_type,
                data,
            )


def get_phenotype_gene_df(version: Optional[str] = None) -> pd.DataFrame:
    """Get the phenotype-gene annotations as a preprocessed dataframe.

    Parameters
    ----------
    version :
        An optional version in form of a string YYYY-MM-DD. If not given,
        gets the latest data.
    Returns
    -------
    :
        A parsed and preprocessed phenotype-gene annotation dataframe
    """
    if version is None:
        url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype_to_genes.txt"
    else:
        url = f"http://purl.obolibrary.org/obo/hp/releases/hpoa/{version}/phenotype_to_genes.txt"
    df = pd.read_csv(
        url,
        sep="\t",
        skiprows=2,
        names=["phenotype", "ncbigene", "source"],
        usecols=[0, 2, 5],
        dtype=str,
    )
    df.drop_duplicates(inplace=True)
    return process_phenotype_gene(df)


def process_phenotype_gene(df) -> pd.DataFrame:
    """"""
    df = process_phenotypes(df)
    df["hgnc_id"] = df["ncbigene"].map(hgnc_client.get_hgnc_from_entrez)
    df = df[df.hgnc_id.notna()]
    return df
