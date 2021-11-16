# -*- coding: utf-8 -*-

"""Base classes for processors."""

import csv
import gzip
import logging
import pickle
from abc import ABC, abstractmethod
from collections import Counter
from pathlib import Path
from typing import ClassVar, Iterable, List, Tuple

import click
import typing

import pystow
from more_click import verbose_option
from tqdm import tqdm

from indra.statements.validate import assert_valid_db_refs
from indra_cogex.representation import Node, Relation, norm_id

__all__ = [
    "Processor",
]

logger = logging.getLogger(__name__)

# deal with importing from wherever with
#  https://stackoverflow.com/questions/36922843/neo4j-3-x-load-csv-absolute-file-path
# Find neo4j conf and comment out this line: dbms.directories.import=import
# /usr/local/Cellar/neo4j/4.1.3/libexec/conf/neo4j.conf for me
# data stored in /usr/local/var/neo4j/data/databases


class Processor(ABC):
    """A processor creates nodes and iterables to upload to Neo4j."""

    name: ClassVar[str]
    description: ClassVar[typing.Mapping[str, str]]
    module: ClassVar[pystow.Module]
    directory: ClassVar[Path]
    nodes_path: ClassVar[Path]
    nodes_indra_path: ClassVar[Path]
    edges_path: ClassVar[Path]
    edges_summary_path: ClassVar[Path]
    importable = True

    def __init_subclass__(cls, **kwargs):
        """Initialize the class attributes."""
        if not getattr(cls, "description"):
            raise TypeError("All processors need description dictionaries")
        cls.module = pystow.module("indra", "cogex", cls.name)
        cls.directory = cls.module.base
        # These are nodes directly in the neo4j encoding
        cls.nodes_path = cls.module.join(name="nodes.tsv.gz")
        # These are nodes in the original INDRA-oriented representation
        # needed for assembly
        cls.nodes_indra_path = cls.module.join(name="nodes.pkl")
        cls.edges_path = cls.module.join(name="edges.tsv.gz")
        cls.edges_summary_path = cls.module.join(name="edges_summary.tsv")

    @abstractmethod
    def get_nodes(self) -> Iterable[Node]:
        """Iterate over the nodes to upload."""

    @abstractmethod
    def get_relations(self) -> Iterable[Relation]:
        """Iterate over the relations to upload."""

    @classmethod
    def get_cli(cls) -> click.Command:
        """Get the CLI for this processor."""

        @click.command()
        @verbose_option
        def _main():
            click.secho(f"Building {cls.name}", fg="green", bold=True)
            processor = cls()
            processor.dump()

        return _main

    @classmethod
    def cli(cls) -> None:
        """Run the CLI for this processor."""
        cls.get_cli()()

    def dump(self) -> Tuple[Path, List[Node], Path, typing.Counter[Tuple[str, str, str]]]:
        """Dump the contents of this processor to CSV files ready for use in ``neo4-admin import``."""
        node_paths, nodes = self._dump_nodes()
        edge_paths, counter = self._dump_edges()
        return node_paths, nodes, edge_paths, counter

    def _dump_nodes(self) -> Tuple[Path, List[Node]]:
        sample_path = self.module.join(name="nodes_sample.tsv")
        nodes = sorted(self.get_nodes(), key=lambda x: (x.db_ns, x.db_id))
        with open(self.nodes_indra_path, "wb") as fh:
            pickle.dump(nodes, fh)
        return self._dump_nodes_to_path(nodes, self.nodes_path, sample_path), nodes

    @staticmethod
    def _dump_nodes_to_path(nodes: Iterable[Node], nodes_path, sample_path=None):
        logger.info(f"Dumping into {nodes_path}...")
        nodes = list(validate_nodes(nodes))
        metadata = sorted(set(key for node in nodes for key in node.data))
        node_rows = (
            (
                norm_id(node.db_ns, node.db_id),
                ";".join(node.labels),
                *[node.data.get(key, "") for key in metadata],
            )
            for node in tqdm(nodes, desc="Nodes", unit_scale=True)
        )

        header = "id:ID", ":LABEL", *metadata
        with gzip.open(nodes_path, mode="wt") as node_file:
            node_writer = csv.writer(node_file, delimiter="\t")  # type: ignore
            node_writer.writerow(header)
            if sample_path:
                with sample_path.open("w") as node_sample_file:
                    node_sample_writer = csv.writer(node_sample_file, delimiter="\t")
                    node_sample_writer.writerow(header)
                    for _, node_row in zip(range(10), node_rows):
                        node_sample_writer.writerow(node_row)
                        node_writer.writerow(node_row)
            # Write remaining nodes
            node_writer.writerows(node_rows)

        return nodes_path

    def _dump_edges(self) -> Tuple[Path, typing.Counter[Tuple[str, str, str]]]:
        sample_path = self.module.join(name="edges_sample.tsv")
        logger.info(f"Dumping into {self.edges_path}...")
        rels = self.get_relations()
        rels = validate_relations(rels)
        rels = sorted(
            rels, key=lambda r: (r.source_ns, r.source_id, r.target_ns, r.target_id)
        )
        metadata = sorted(set(key for rel in rels for key in rel.data))
        edge_rows = (
            (
                norm_id(rel.source_ns, rel.source_id),
                norm_id(rel.target_ns, rel.target_id),
                rel.rel_type,
                *[rel.data.get(key) for key in metadata],
            )
            for rel in tqdm(rels, desc="Edges", unit_scale=True)
        )

        counter = Counter()
        with gzip.open(self.edges_path, "wt") as edge_file:
            edge_writer = csv.writer(edge_file, delimiter="\t")  # type: ignore
            with sample_path.open("w") as edge_sample_file:
                edge_sample_writer = csv.writer(edge_sample_file, delimiter="\t")
                header = ":START_ID", ":END_ID", ":TYPE", *metadata
                edge_sample_writer.writerow(header)
                edge_writer.writerow(header)

                for _, edge_row in zip(range(10), edge_rows):
                    edge_sample_writer.writerow(edge_row)
                    edge_writer.writerow(edge_row)
                    counter[edge_row[2]] += 1

            # Write remaining edges
            for edge_row in edge_rows:
                edge_writer.writerow(edge_row)
                counter[edge_row[2]] += 1

        with self.edges_summary_path.open("w") as file:
            print("source_ns", "rel", "target_ns", "count", sep="\t", file=file)
            for row, count in counter.most_common():
                print(*row, count, sep="\t", file=file)

        return self.edges_path, counter


def validate_nodes(nodes: Iterable[Node]) -> Iterable[Node]:
    for idx, node in enumerate(nodes):
        try:
            assert_valid_db_refs({node.db_ns: node.db_id})
            yield node
        except Exception as e:
            logger.info(f"{idx}: {node} - {e}")
            continue


def validate_relations(relations):
    for idx, rel in enumerate(relations):
        try:
            assert_valid_db_refs({rel.source_ns: rel.source_id})
            assert_valid_db_refs({rel.target_ns: rel.target_id})
            yield rel
        except Exception as e:
            logger.info(f"{idx}: {rel} - {e}")
            continue
