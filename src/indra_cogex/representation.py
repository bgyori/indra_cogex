# -*- coding: utf-8 -*-

"""Representations for nodes and relations to upload to Neo4j."""


__all__ = ["Node", "Relation", "indra_stmts_from_relations"]

from typing import (
    Any,
    Collection,
    Iterable,
    List,
    Mapping,
    Optional,
    Tuple,
    Dict,
    Union,
)
import json
from indra.databases import identifiers
from indra.ontology.standardize import standardize_name_db_refs
from indra.statements.agent import get_grounding
from indra.statements import stmts_from_json, Statement

NodeJson = Dict[str, Union[Collection[str], Dict[str, Any]]]
RelJson = Dict[str, Union[Mapping[str, Any], Dict]]


class Node:
    """Representation for a node."""

    def __init__(
        self,
        db_ns: str,
        db_id: str,
        labels: Collection[str],
        data: Optional[Mapping[str, Any]] = None,
    ):
        """Initialize the node.

        Parameters
        ----------
        db_ns :
            The namespace associated with the node. Uses the INDRA standard.
        db_id :
            The identifier within the namespace associated with the node.
            Uses the INDRA standard.
        labels :
            A collection of labels for the node.
        data :
            An optional data dictionary associated with the node.
        """
        if not db_ns or not db_id:
            raise ValueError("Missing namespace or ID.")
        self.db_ns = db_ns
        self.db_id = db_id
        self.labels = labels
        self.data = data if data else {}

    @classmethod
    def standardized(
        cls,
        *,
        db_ns: str,
        db_id: str,
        name: Optional[str] = None,
        labels: Collection[str],
    ) -> "Node":
        """Initialize the node, but first standardize the prefix/identifier/name.

        Parameters
        ----------
        db_ns :
            The namespace associated with the node.
        db_id :
            The identifier within the namespace associated with the node.
        name :
            An optional name for the node.
        labels :
            A collection of labels for the node.

        Returns
        -------
        :
            A node with standardized prefix/identifier/name.
        """
        db_ns, db_id, name = standardize(db_ns, db_id, name)
        return cls(
            db_ns,
            db_id,
            labels,
            dict(name=name),
        )

    def grounding(self) -> Tuple[str, str]:
        """Get the grounded namespace and identifier for this node as a tuple

        Returns
        -------
        :
            A tuple of the namespace and identifier for the node.
        """
        return self.db_ns, self.db_id

    def to_json(self) -> NodeJson:
        """Serialize the node to JSON.

        Returns
        -------
        :
            A JSON representation of the node.
        """
        data = {k: v for k, v in self.data.items()}
        data["db_ns"] = self.db_ns
        data["db_id"] = self.db_id
        # Fixme: how to properly serialize labels?
        return {"labels": [lb for lb in self.labels], "data": data}

    def _get_data_str(self) -> str:
        pieces = ["id:'%s:%s'" % (self.db_ns, self.db_id)]
        for k, v in self.data.items():
            if isinstance(v, str):
                value = "'" + v.replace("'", "\\'") + "'"
            elif isinstance(v, (bool, int, float)):
                value = v
            else:
                value = str(v)
            piece = "%s:%s" % (k, value)
            pieces.append(piece)
        data_str = ", ".join(pieces)
        return data_str

    def __str__(self):  # noqa:D105
        data_str = self._get_data_str()
        labels_str = ":".join(self.labels)
        return f"(:{labels_str} {{ {data_str} }})"

    def __repr__(self):  # noqa:D105
        return str(self)


class Relation:
    """Representation for a relation."""

    def __init__(
        self,
        source_ns: str,
        source_id: str,
        target_ns: str,
        target_id: str,
        rel_type: str,
        data: Optional[Mapping[str, Any]] = None,
    ):
        """Initialize the relation.

        Parameters
        ----------
        source_ns :
            The namespace associated with the source node.
        source_id :
            The identifier within the namespace associated with the source node.
        target_ns :
            The namespace associated with the target node.
        target_id :
            The identifier within the namespace associated with the target node.
        rel_type :
            The type of relation.
        data :
            An optional data dictionary associated with the relation.
        """
        self.source_ns = source_ns
        self.source_id = source_id
        self.target_ns = target_ns
        self.target_id = target_id
        self.rel_type = rel_type
        self.data = data if data else {}

    def to_json(self) -> RelJson:
        """Serialize the relation to JSON format.

        Returns
        -------
        :
            A JSON representation of the relation.
        """
        return {
            "source_ns": self.source_ns,
            "source_id": self.source_id,
            "target_ns": self.target_ns,
            "target_id": self.target_id,
            "rel_type": self.rel_type,
            "data": self.data,
        }

    def __str__(self):  # noqa:D105
        data_str = ", ".join(["%s:'%s'" % (k, v) for k, v in self.data.items()])
        return (
            f"({self.source_ns}, {self.source_id})-[:{self.rel_type} {data_str}]->"
            f"({self.target_ns}, {self.target_id})"
        )

    def __repr__(self):  # noqa:D105
        return str(self)


def standardize(
    prefix: str, identifier: str, name: Optional[str] = None
) -> Tuple[str, str, str]:
    """Get a standardized prefix, identifier, and name, if possible.

    Parameters
    ----------
    prefix :
        The prefix to standardize.
    identifier :
        The identifier to standardize.
    name :
        The name to standardize.

    Returns
    -------
    :
        A tuple of the standardized prefix, identifier, and name.
    """
    standard_name, db_refs = standardize_name_db_refs({prefix: identifier})
    name = standard_name if standard_name else name
    db_ns, db_id = get_grounding(db_refs)
    if db_ns is None or db_id is None:
        return prefix, identifier, name
    return db_ns, db_id, name


def norm_id(db_ns, db_id) -> str:
    """Normalize an identifier.

    Parameters
    ----------
    db_ns :
        The namespace of the identifier.
    db_id :
        The identifier.

    Returns
    -------
    :
        The normalized identifier.
    """
    identifiers_ns = identifiers.get_identifiers_ns(db_ns)
    identifiers_id = db_id
    if not identifiers_ns:
        identifiers_ns = db_ns.lower()
    else:
        ns_embedded = identifiers.identifiers_registry.get(identifiers_ns, {}).get(
            "namespace_embedded"
        )
        if ns_embedded:
            identifiers_id = identifiers_id[len(identifiers_ns) + 1 :]
    return f"{identifiers_ns}:{identifiers_id}"


def triple_query(
    source_name: Optional[str] = None,
    source_type: Optional[str] = None,
    source_id: Optional[str] = None,
    relation_name: Optional[str] = None,
    relation_type: Optional[str] = None,
    target_name: Optional[str] = None,
    target_type: Optional[str] = None,
    target_id: Optional[str] = None,
    relation_direction: Optional[str] = "right",
) -> str:
    """Create a Cypher query from the given parameters.

    Parameters
    ----------
    source_name :
        The name of the source node. Optional.
    source_type :
        The type of the source node. Optional.
    source_id :
        The identifier of the source node. Optional.
    relation_name :
        The name of the relation. Optional.
    relation_type :
        The type of the relation. Optional.
    target_name :
        The name of the target node. Optional.
    target_type :
        The type of the target node. Optional.
    target_id :
        The identifier of the target node. Optional.
    relation_direction :
        The direction of the relation, one of 'left', 'right', or 'both'.
        These correspond to <-[]-, -[]->, and -[]-, respectively.

    Returns
    -------
    :
        A Cypher query as a string.
    """
    rel1, rel2 = "-", "-"
    if relation_direction == "left":
        rel1 = "<-"
    elif relation_direction == "right":
        rel2 = "->"
    source = node_query(node_name=source_name, node_type=source_type, node_id=source_id)
    # TODO could later make an alternate function for the relation
    relation = node_query(node_name=relation_name, node_type=relation_type)
    target = node_query(node_name=target_name, node_type=target_type, node_id=target_id)
    return f"({source}){rel1}[{relation}]{rel2}({target})"


def node_query(
    node_name: Optional[str] = None,
    node_type: Optional[str] = None,
    node_id: Optional[str] = None,
) -> str:
    """Create a Cypher node query

    Parameters
    ----------
    node_name :
        The name of the node. Optional.
    node_type :
        The type of the node. Optional.
    node_id :
        The identifier of the node. Optional.

    Returns
    -------
    :
        A Cypher node query as a string.
    """
    if node_name is None:
        node_name = ""
    rv = node_name or ""
    if node_type:
        rv += f":{node_type}"
    if node_id:
        if rv:
            rv += " "
        rv += f"{{id: '{node_id}'}}"
    return rv


def indra_stmts_from_relations(rels: Iterable[Relation]) -> List[Statement]:
    """Convert a list of relations to INDRA Statements.

    Any relations that aren't representing an INDRA Statement are skipped.

    Parameters
    ----------
    rels :
        A list of Relations.

    Returns
    -------
    :
        A list of INDRA Statements.
    """
    stmts_json = [
        json.loads(rel.data["stmt_json"].replace("\\\\", "\\").replace("\\\\", "\\"))
        for rel in rels
    ]
    stmts = stmts_from_json(stmts_json)
    return stmts
