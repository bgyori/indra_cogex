# -*- coding: utf-8 -*-

"""A collection of analyses possible on gene lists (of HGNC identifiers) with scores.

For example, this could be applied to the log_2 fold scores from differential gene
expression experiments.

.. warning::

    This module requires the optional dependency ``gseapy``. Install with
    ``pip install gseapy``.
"""

from pathlib import Path
from typing import Dict, Optional, Set, Tuple, Union

import gseapy
import pandas as pd
import pyobo
from indra.databases import hgnc_client

from indra_cogex.client.enrichment.utils import (
    get_entity_to_regulators,
    get_entity_to_targets,
    get_go,
    get_phenotype_gene_sets,
    get_reactome,
    get_wikipathways,
)
from indra_cogex.client.neo4j_client import Neo4jClient, autoclient

__all__ = [
    "get_rat_scores",
    "get_mouse_scores",
    "get_human_scores",
    "go_gsea",
    "wikipathways_gsea",
    "reactome_gsea",
    "phenotype_gsea",
    "indra_upstream_gsea",
    "indra_downstream_gsea",
    "gsea",
]

GENE_SYMBOL_COLUMN_GUESSES = [
    "gene_name",
]
SCORE_COLUMN_GUESSES = [
    "log2FoldChange",
]


def _guess_symbol_col(df: pd.DataFrame) -> str:
    for guess in GENE_SYMBOL_COLUMN_GUESSES:
        if guess in df.columns:
            return guess
    raise ValueError(f"could not guess gene symbol column name from: {df.columns}")


def _guess_score_col(df: pd.DataFrame) -> str:
    for guess in SCORE_COLUMN_GUESSES:
        if guess in df.columns:
            return guess
    raise ValueError(f"could not guess score column name from: {df.columns}")


def get_rat_scores(
    path: Union[Path, str, pd.DataFrame],
    read_csv_kwargs: Optional[Dict[str, any]] = None,
    gene_symbol_column_name: Optional[str] = None,
    score_column_name: Optional[str] = None,
) -> Dict[str, float]:
    """Load a differential gene expression file with rat measurements.

    This function extracts the RGD gene symbols, maps them
    to RGD identifiers, uses PyOBO to map orthologs to HGNC,
    then returns the HGNC gene and scores as a dictionary.

    Parameters
    ----------
    path :
        Path to the file to read with :func:`pandas.read_csv`.
    read_csv_kwargs :
        Keyword arguments to pass to :func:`pandas.read_csv`
    gene_symbol_column_name :
        The name of the column with gene symbols. If none,
        will try and guess.
    score_column_name :
        The name of the column with scores. If none, will try
        and guess.

    Returns
    -------
    :
        A dictionary of mapped orthologus human gene HGNC IDs to
        scores.
    """
    return _get_species_scores(
        prefix="rgd",
        func=hgnc_client.get_hgnc_from_rat,
        path=path,
        read_csv_kwargs=read_csv_kwargs,
        gene_symbol_column_name=gene_symbol_column_name,
        score_column_name=score_column_name,
    )


def get_mouse_scores(
    path: Union[Path, str, pd.DataFrame],
    read_csv_kwargs: Optional[Dict[str, any]] = None,
    gene_symbol_column_name: Optional[str] = None,
    score_column_name: Optional[str] = None,
) -> Dict[str, float]:
    """Load a differential gene expression file with mouse measurements.

    This function extracts the MGI gene symbols, maps them
    to MGI identifiers, uses PyOBO to map orthologs to HGNC,
    then returns the HGNC gene and scores as a dictionary.

    Parameters
    ----------
    path :
        Path to the file to read with :func:`pandas.read_csv`.
    read_csv_kwargs :
        Keyword arguments to pass to :func:`pandas.read_csv`
    gene_symbol_column_name :
        The name of the column with gene symbols. If none,
        will try and guess.
    score_column_name :
        The name of the column with scores. If none, will try
        and guess.

    Returns
    -------
    :
        A dictionary of mapped orthologus human gene HGNC IDs to
        scores.
    """
    return _get_species_scores(
        prefix="mgi",
        func=hgnc_client.get_hgnc_from_mouse,
        path=path,
        read_csv_kwargs=read_csv_kwargs,
        gene_symbol_column_name=gene_symbol_column_name,
        score_column_name=score_column_name,
    )


def get_human_scores(
    path: Union[Path, str, pd.DataFrame],
    read_csv_kwargs: Optional[Dict[str, any]] = None,
    gene_symbol_column_name: Optional[str] = None,
    score_column_name: Optional[str] = None,
) -> Dict[str, float]:
    """Load a differential gene expression file with human measurements.

    Parameters
    ----------
    path :
        Path to the file to read with :func:`pandas.read_csv`.
    read_csv_kwargs :
        Keyword arguments to pass to :func:`pandas.read_csv`
    gene_symbol_column_name :
        The name of the column with gene symbols. If none,
        will try and guess.
    score_column_name :
        The name of the column with scores. If none, will try
        and guess.

    Returns
    -------
    :
        A dictionary of human gene HGNC IDs to scores.
    """
    return _get_species_scores(
        path=path,
        read_csv_kwargs=read_csv_kwargs,
        gene_symbol_column_name=gene_symbol_column_name,
        score_column_name=score_column_name,
    )


def _get_species_scores(
    path: Union[Path, str, pd.DataFrame],
    read_csv_kwargs: Optional[Dict[str, any]] = None,
    gene_symbol_column_name: Optional[str] = None,
    score_column_name: Optional[str] = None,
    *,
    prefix=None,
    func=None,
) -> Dict[str, float]:
    if isinstance(path, pd.DataFrame):
        df = path
    else:
        df = pd.read_csv(path, **(read_csv_kwargs or {}))
    if gene_symbol_column_name is None:
        gene_symbol_column_name = _guess_symbol_col(df)
    elif gene_symbol_column_name not in df.columns:
        raise ValueError
    if score_column_name is None:
        score_column_name = _guess_score_col(df)
    elif score_column_name not in df.columns:
        raise ValueError

    if prefix is not None and func is not None:
        mapped_gene_symbol_column_name = f"{prefix}_id"
        df[mapped_gene_symbol_column_name] = df[gene_symbol_column_name].map(
            pyobo.get_name_id_mapping(prefix)
        )
        df = df[df[mapped_gene_symbol_column_name].notna()]
    elif prefix is not None or func is not None:
        raise ValueError("If specifying one, must specify both of prefix and func")
    else:
        # If no prefix is given, assume columns are human.
        mapped_gene_symbol_column_name = gene_symbol_column_name
        func = hgnc_client.get_current_hgnc_id

    df["hgnc_id"] = df[mapped_gene_symbol_column_name].map(func)
    df = df.set_index("hgnc_id")
    return df[score_column_name].to_dict()


@autoclient()
def wikipathways_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA with WikiPathways gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_wikipathways(client=client),
        scores=scores,
        directory=directory,
        **kwargs,
    )


@autoclient()
def reactome_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA with Reactome gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_reactome(client=client),
        scores=scores,
        directory=directory,
        **kwargs,
    )


@autoclient()
def phenotype_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA with HPO phenotype gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_phenotype_gene_sets(client=client),
        scores=scores,
        directory=directory,
        **kwargs,
    )


@autoclient()
def go_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA with gene sets for each Gene Ontology term.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_go(client=client),
        scores=scores,
        directory=directory,
        **kwargs,
    )


@autoclient()
def indra_upstream_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    minimum_evidence_count: Optional[int] = None,
    minimum_belief: Optional[float] = None,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA for each entry in the INDRA database and the set
    of human genes that it regulates.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    minimum_evidence_count :
        The minimum number of evidences for a relationship to count it as a regulator.
        Defaults to 1 (i.e., cutoff not applied.
    minimum_belief :
        The minimum belief for a relationship to count it as a regulator.
        Defaults to 0.0 (i.e., cutoff not applied).
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_entity_to_targets(
            client=client,
            minimum_evidence_count=minimum_evidence_count,
            minimum_belief=minimum_belief,
        ),
        scores=scores,
        directory=directory,
        **kwargs,
    )


@autoclient()
def indra_downstream_gsea(
    scores: Dict[str, float],
    directory: Union[None, Path, str] = None,
    *,
    client: Neo4jClient,
    minimum_evidence_count: Optional[int] = None,
    minimum_belief: Optional[float] = None,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA for each entry in the INDRA database and the set
    of human genes that are upstream regulators of it.

    Parameters
    ----------
    client :
        The Neo4j client.
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    minimum_evidence_count :
        The minimum number of evidences for a relationship to count it as a regulator.
        Defaults to 1 (i.e., cutoff not applied.
    minimum_belief :
        The minimum belief for a relationship to count it as a regulator.
        Defaults to 0.0 (i.e., cutoff not applied).
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    return gsea(
        gene_sets=get_entity_to_regulators(
            client=client,
            minimum_evidence_count=minimum_evidence_count,
            minimum_belief=minimum_belief,
        ),
        scores=scores,
        directory=directory,
        **kwargs,
    )


GSEA_RETURN_COLUMNS = [
    "term",
    "name",
    "es",
    "nes",
    "pval",
    "fdr",
    "geneset_size",
    "matched_size",
]


def gsea(
    scores: Dict[str, float],
    gene_sets: Dict[Tuple[str, str], Set[str]],
    directory: Union[None, Path, str] = None,
    alpha: Optional[float] = None,
    keep_insignificant: bool = True,
    **kwargs,
) -> pd.DataFrame:
    """Run GSEA on pre-ranked data.

    Parameters
    ----------
    scores :
        A mapping from HGNC gene identifiers to floating point scores
        (e.g., from a differential gene expression analysis)
    gene_sets :
        A mapping from
    directory :
        Specify the directory if the results should be saved, including
        both a dataframe and plots for each gen set
    alpha :
        The cutoff for significance. Defaults to 0.05
    keep_insignificant :
        If false, removes results with a p value less than alpha.
    kwargs :
        Remaining keyword arguments to pass through to :func:`gseapy.prerank`

    Returns
    -------
    :
        A pandas dataframe with the GSEA results
    """
    if alpha is None:
        alpha = 0.05
    if directory is not None:
        if isinstance(directory, str):
            directory = Path(directory)
        directory.mkdir(exist_ok=True, parents=True)
        directory = directory.as_posix()

    curie_to_name = dict(gene_sets.keys())
    curie_to_gene_sets = {
        curie: hgnc_gene_ids for (curie, _), hgnc_gene_ids in gene_sets.items()
    }

    kwargs.setdefault("permutation_num", 100)
    kwargs.setdefault("format", "svg")
    res = gseapy.prerank(
        rnk=pd.Series(scores),
        gene_sets=curie_to_gene_sets,
        outdir=directory,
        **kwargs,
    )
    res.res2d.index.name = "term"
    rv = res.res2d.reset_index()
    rv["name"] = rv["term"].map(curie_to_name)
    rv = rv[GSEA_RETURN_COLUMNS]
    if not keep_insignificant:
        rv = rv[rv["pval"] < alpha]
    return rv
