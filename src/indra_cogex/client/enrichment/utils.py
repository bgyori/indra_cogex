# -*- coding: utf-8 -*-

"""Utilities for getting gene sets."""

import logging
from collections import defaultdict
from textwrap import dedent
from typing import Dict, Optional, Set, Tuple

from indra_cogex.client.neo4j_client import Neo4jClient, autoclient

__all__ = [
    "collect_gene_sets",
    "get_go",
    "get_wikipathways",
    "get_reactome",
    "get_phenotype_gene_sets",
    "get_entity_to_targets",
    "get_entity_to_regulators",
]

logger = logging.getLogger(__name__)


@autoclient()
def collect_gene_sets(
    query: str,
    *,
    client: Neo4jClient,
) -> Dict[Tuple[str, str], Set[str]]:
    """Collect gene sets based on the given query.

    Parameters
    ----------
    query:
        A cypher query
    client :
        The Neo4j client.

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each queried
        item and whose values are sets of HGNC gene identifiers (as strings)
    """
    curie_to_hgnc_ids = defaultdict(set)
    for result in client.query_tx(query):
        curie = result[0]
        name = result[1]
        hgnc_ids = {
            hgnc_curie.lower().replace("hgnc:", "")
            if hgnc_curie.lower().startswith("hgnc:")
            else hgnc_curie.lower()
            for hgnc_curie in result[2]
        }
        curie_to_hgnc_ids[curie, name].update(hgnc_ids)
    return dict(curie_to_hgnc_ids)


@autoclient(cache=True)
def get_go(*, client: Neo4jClient) -> Dict[Tuple[str, str], Set[str]]:
    """Get GO gene sets.

    Parameters
    ----------
    client :object
        The Neo4j client.

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each GO term
        and whose values are sets of HGNC gene identifiers (as strings)
    """
    query = dedent(
        """\
        MATCH (gene:BioEntity)-[:associated_with]->(term:BioEntity)
        RETURN term.id, term.name, collect(gene.id) as gene_curies;
    """
    )
    logger.info("caching GO with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


@autoclient(cache=True)
def get_wikipathways(*, client: Neo4jClient) -> Dict[Tuple[str, str], Set[str]]:
    """Get WikiPathways gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each WikiPathway
        pathway and whose values are sets of HGNC gene identifiers (as strings)
    """
    query = dedent(
        """\
        MATCH (pathway:BioEntity)-[:haspart]->(gene:BioEntity)
        WHERE pathway.id STARTS WITH "wikipathways" and gene.id STARTS WITH "hgnc"
        RETURN pathway.id, pathway.name, collect(gene.id);
    """
    )
    logger.info("caching WikiPathways with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


@autoclient(cache=True)
def get_reactome(*, client: Neo4jClient) -> Dict[Tuple[str, str], Set[str]]:
    """Get Reactome gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each Reactome
        pathway and whose values are sets of HGNC gene identifiers (as strings)
    """
    query = dedent(
        """\
        MATCH (pathway:BioEntity)-[:haspart]-(gene:BioEntity)
        WHERE pathway.id STARTS WITH "reactome" and gene.id STARTS WITH "hgnc"
        RETURN pathway.id, pathway.name, collect(gene.id);
    """
    )
    logger.info("caching Reactome with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


@autoclient(cache=True)
def get_phenotype_gene_sets(*, client: Neo4jClient) -> Dict[Tuple[str, str], Set[str]]:
    """Get HPO phenotype gene sets.

    Parameters
    ----------
    client :
        The Neo4j client.

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each phenotype
        gene set and whose values are sets of HGNC gene identifiers (as strings)
    """
    query = dedent(
        """\
        MATCH (s:BioEntity)-[:phenotype_has_gene]-(gene:BioEntity)
        WHERE s.id STARTS WITH "HP" and gene.id STARTS WITH "hgnc"
        RETURN s.id, s.name, collect(gene.id);
    """
    )
    logger.info("caching phenotype gene sets with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


@autoclient(cache=True, maxsize=5)
def get_entity_to_targets(
    *,
    client: Neo4jClient,
    minimum_evidence_count: Optional[int] = None,
    minimum_belief: Optional[float] = None,
) -> Dict[Tuple[str, str], Set[str]]:
    """Get a mapping from each entity in the INDRA database to the set of
    human genes that it regulates.

    Parameters
    ----------
    client :
        The Neo4j client.
    minimum_evidence_count :
        The minimum number of evidences for a relationship to count it as a regulator.
        Defaults to 1 (i.e., cutoff not applied.
    minimum_belief :
        The minimum belief for a relationship to count it as a regulator.
        Defaults to 0.0 (i.e., cutoff not applied).

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each entity
        and whose values are sets of HGNC gene identifiers (as strings)
    """
    evidence_line = minimum_evidence_helper(
        minimum_evidence_count=minimum_evidence_count, name="r"
    )
    belief_line = minimum_belief_helper(minimum_belief=minimum_belief, name="r")
    query = dedent(
        f"""\
        MATCH (regulator:BioEntity)-[r:indra_rel]->(gene:BioEntity)
        WHERE
            gene.id STARTS WITH "hgnc"                  // Collecting human genes only
            AND r.stmt_type <> "Complex"                // Ignore complexes since they are non-directional
            AND NOT regulator.id STARTS WITH "uniprot"  // This is a simple way to ignore non-human proteins
            {evidence_line}
            {belief_line}
        RETURN
            regulator.id,
            regulator.name,
            collect(gene.id);
    """
    )
    logger.info("caching entity->targets with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


@autoclient(cache=True, maxsize=5)
def get_entity_to_regulators(
    *,
    client: Neo4jClient,
    minimum_evidence_count: Optional[int] = None,
    minimum_belief: Optional[float] = None,
) -> Dict[Tuple[str, str], Set[str]]:
    """Get a mapping from each entity in the INDRA database to the set of
    human genes that are causally upstream of it.

    Parameters
    ----------
    client :
        The Neo4j client.
    minimum_evidence_count :
        The minimum number of evidences for a relationship to count it as a regulator.
        Defaults to 1 (i.e., cutoff not applied.
    minimum_belief :
        The minimum belief for a relationship to count it as a regulator.
        Defaults to 0.0 (i.e., cutoff not applied).

    Returns
    -------
    :
        A dictionary whose keys that are 2-tuples of CURIE and name of each entity
        and whose values are sets of HGNC gene identifiers (as strings)
    """
    evidence_line = minimum_evidence_helper(
        minimum_evidence_count=minimum_evidence_count, name="r"
    )
    belief_line = minimum_belief_helper(minimum_belief=minimum_belief, name="r")
    query = dedent(
        f"""\
        MATCH (gene:BioEntity)-[r:indra_rel]->(target:BioEntity)
        WHERE
            gene.id STARTS WITH "hgnc"               // Collecting human genes only
            AND r.stmt_type <> "Complex"             // Ignore complexes since they are non-directional
            AND NOT target.id STARTS WITH "uniprot"  // This is a simple way to ignore non-human proteins
            {evidence_line}
            {belief_line}
        RETURN
            target.id,
            target.name,
            collect(gene.id);
    """
    )
    logger.info("caching entity->regulators with Cypher query: %s", query)
    return collect_gene_sets(client=client, query=query)


def minimum_evidence_helper(
    minimum_evidence_count: Optional[float] = None, name: str = "r"
) -> str:
    if minimum_evidence_count is None or minimum_evidence_count == 1:
        return ""
    return f"AND {name}.evidence_count >= {minimum_evidence_count}"


def minimum_belief_helper(
    minimum_belief: Optional[float] = None, name: str = "r"
) -> str:
    if minimum_belief is None or minimum_belief == 0.0:
        return ""
    return f"AND {name}.belief >= {minimum_belief}"
