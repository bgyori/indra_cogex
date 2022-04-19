from collections import Counter

from tabulate import tabulate

import pyobo
from biomappings import load_mappings
from indra.databases.identifiers import get_ns_id_from_identifiers
from indra_cogex.representation import Node, standardize

__all__ = [
    "UmlsMapper",
]


class UmlsMapper:
    """A utility class for mapping out of UMLS."""

    prefixes = ["doid", "mesh", "hp", "efo", "mondo"]

    def __init__(self):
        """Prepare the UMLS mappings from PyOBO and Biomappings."""
        #: A dictionary from external prefix to UMLS id to external ID
        self.xrefs = {}

        for prefix in self.prefixes:
            self.xrefs[prefix] = {}
            # Get external to UMLS
            for external_id, umls_id in pyobo.get_filtered_xrefs(
                prefix, "umls"
            ).items():
                self.xrefs[prefix][umls_id] = external_id
            # Get UMLS to external
            for umls_id, external_id in pyobo.get_filtered_xrefs(
                "umls", prefix
            ).items():
                self.xrefs[prefix][umls_id] = external_id

        # Get manually curated UMLS mappings from biomappings
        biomappings_from_umls, biomappings_to_umls = Counter(), Counter()
        for mapping in load_mappings():
            if mapping["source prefix"] == "umls":
                target_prefix = mapping["target prefix"]
                biomappings_from_umls[target_prefix] += 1
                target_id = mapping["target identifier"]
                source_id = mapping["source identifier"]
                if target_prefix in self.xrefs:
                    self.xrefs[target_prefix][target_id] = source_id
                else:
                    self.xrefs[target_prefix] = {
                        target_id: source_id,
                    }
            elif mapping["target prefix"] == "umls":
                source_prefix = mapping["source prefix"]
                biomappings_to_umls[source_prefix] += 1
                source_id = mapping["source identifier"]
                target_id = mapping["target identifier"]
                if source_prefix in self.xrefs:
                    self.xrefs[source_prefix][source_id] = target_id
                else:
                    self.xrefs[source_prefix] = {
                        source_id: target_id,
                    }

        print("Mapping out of UMLS")
        print(tabulate(biomappings_from_umls.most_common()))
        print("Mapping into UMLS")
        print(tabulate(biomappings_to_umls.most_common()))

        print("Total xrefs")
        print(
            tabulate(
                [(prefix, len(self.xrefs[prefix])) for prefix in self.prefixes],
                headers=["Prefix", "Mappings"],
            )
        )

    def lookup(self, umls_id: str):
        for prefix in self.prefixes:
            xrefs = self.xrefs[prefix]
            identifier = xrefs.get(umls_id)
            if identifier is not None:
                return standardize(prefix, identifier)
        return "umls", umls_id, pyobo.get_name("umls", umls_id)

    def standardize(self, umls_id: str):
        prefix, identifier, name = self.lookup(umls_id)
        db_ns, db_id = get_ns_id_from_identifiers(prefix, identifier)
        if db_ns is None:
            db_ns, db_id = prefix, identifier
        return db_ns, db_id, name
