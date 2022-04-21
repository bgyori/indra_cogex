"""
A class to manage curations and the curation cache.

The curation cache should handle:
    - Loading curations to cache in memory
    - Submitting curations to the curation database
    - Update the in memory cache when curations are submitted or at regular
      intervals
"""

from collections import defaultdict
from datetime import datetime, timedelta
from typing import DefaultDict, Dict, List, Optional, Set, Union

import dateutil.parser
import pandas as pd
from indra.sources.indra_db_rest import get_curations, submit_curation

__all__ = [
    "CurationCache",
    "Curation",
    "Curations",
]

Curation = Dict[str, Optional[Union[str, int]]]
Curations = List[Curation]


class CurationCache:
    update_interval: timedelta
    last_update: datetime
    curation_list: Curations
    curations_df: pd.DataFrame

    def __init__(
        self,
        update_interval: timedelta = timedelta(minutes=30),
    ):
        self.update_interval = update_interval
        self.curation_list = []
        self.refresh_curations()

    def refresh_curations(self):
        """Refresh the curation cache"""
        self.curation_list: Curations = [
            self._process_curation(curation) for curation in get_curations()
        ]
        self.curations_df = self._get_curation_df(self.curation_list)
        self.last_update = datetime.utcnow()

    @staticmethod
    def _get_curation_df(curations) -> pd.DataFrame:
        rv = pd.DataFrame(curations).astype(
            dtype={
                "id": pd.Int64Dtype(),
                "pa_hash": pd.Int64Dtype(),
                "source_hash": pd.Int64Dtype(),
            }
        )
        return rv

    @staticmethod
    def _process_curation(curation) -> Curation:
        curation["date"] = dateutil.parser.parse(curation["date"])
        return curation

    def get_curation_cache(
        self, refresh: bool = False, only_most_recent: bool = False
    ) -> Curations:
        """Get the curations in the cache.

        Parameters
        ----------
        refresh :
            Whether to refresh the curation cache
        only_most_recent :
            Set to true to filter out all but the most recent (based on date)
            for all given curator/statement/evidence triples.

        Returns
        -------
        :
            A list of all curations, potentially filtered if
            ``only_most_recent`` is set to true.
        """
        if refresh:
            self.refresh_curations()
        if not only_most_recent:
            return self.curation_list

        # Aggregate all curations by curator/statement/evidence,
        # then only keep the curation with the most recent datetime
        keys = ("curator", "pa_hash", "source_hash")
        aggregator = defaultdict(list)
        for curation in self.curation_list:
            aggregator[self._curation_key(curation)].append(curation)
        return [
            max(values, key=lambda value: value["date"])
            for values in aggregator.values()
        ]

    @staticmethod
    def _curation_key(curation):
        return tuple(curation[key] for key in ("curator", "pa_hash", "source_hash"))

    def get_curations(
        self,
        pa_hash: Optional[int] = None,
        source_hash: Optional[int] = None,
        refresh: bool = False,
    ) -> Curations:
        """Get curations from the cache based on pa_hash and source_hash

        Note: If all curations are needed, it is more efficient to use
        get_curation_cache() method.

        Parameters
        ----------
        pa_hash :
            The statement hash
        source_hash :
            The source hash of the evidence
        refresh :
            Whether to refresh the curation cache

        Returns
        -------
        :
            A list of curations
        """
        if source_hash is not None and pa_hash is None:
            raise ValueError("Must provide a pa_hash if source_hash is provided")

        # Update the curation cache if it is too old or if asked
        if refresh or self.last_update + self.update_interval < datetime.utcnow():
            self.refresh_curations()

        temp_df = self.curations_df
        if pa_hash is not None:
            temp_df = temp_df[temp_df.pa_hash == pa_hash]
        if source_hash is not None:
            temp_df = temp_df[temp_df.source_hash == source_hash]

        return temp_df.copy().to_dict(orient="records")

    def submit_curation(
        self,
        hash_val: int,
        tag: str,
        email: str,
        text: str,
        ev_hash: int,
        source_api: str,
    ) -> int:
        """Submit a curation to the curation database"""

        dbid = submit_curation(
            hash_val=hash_val,
            tag=tag,
            curator_email=email,
            text=text,
            ev_hash=ev_hash,
            source=source_api,
        )

        # Set last update to older than the update interval to force an
        # update on the next call
        self.last_update = (
            datetime.utcnow() - self.update_interval - timedelta(seconds=1)
        )

        return dbid

    def get_correct_evidence_hashes(self, only_most_recent: bool = False) -> Set[int]:
        """Get a set of all evidence hashes marked as correct.

        Parameters
        ----------
        only_most_recent :
            Set to true to filter out all but the most recent (based on date)
            for all given curator/statement/evidence triples.

        Returns
        -------
        :
            A set of evidence hashes (i.e., from the "source_hash" field) that
            have been marked as correct
        """
        d: DefaultDict[int, Curations] = defaultdict(list)
        for curation in self.get_curation_cache(only_most_recent=only_most_recent):
            d[curation["source_hash"]].append(curation)
        return {
            source_hash
            for source_hash, curations in d.items()
            if any(curation["tag"] == "correct" for curation in curations)
        }

    def get_incorrect_evidence_hashes(self, only_most_recent: bool = False) -> Set[int]:
        """Get a set of all evidence hashes marked as incorrect (undisputed).

        Parameters
        ----------
        only_most_recent :
            Set to true to filter out all but the most recent (based on date)
            for all given curator/statement/evidence triples.

        Returns
        -------
        :
            A set of evidence hashes (i.e., from the "source_hash" field) that
            have been not been marked incorrect
        """
        d: DefaultDict[int, Curations] = defaultdict(list)
        for curation in self.get_curation_cache(only_most_recent=only_most_recent):
            d[curation["source_hash"]].append(curation)
        return {
            source_hash
            for source_hash, curations in d.items()
            if all(curation["tag"] != "correct" for curation in curations)
        }

    def get_curated_evidence_hashes(self, only_most_recent: bool = False) -> Set[int]:
        """Get a set of all evidence hashes.

        Parameters
        ----------
        only_most_recent :
            Set to true to filter out all but the most recent (based on date)
            for all given curator/statement/evidence triples.

        Returns
        -------
        :
            A set of all evidence hashes (i.e., from the "source_hash" field) that
            have been curated
        """
        return {
            curation["source_hash"]
            for curation in self.get_curation_cache(only_most_recent=only_most_recent)
        }

    def get_correct_statement_hashes(self, only_most_recent: bool = False) -> Set[int]:
        """Get a set of all statement hashes marked as correct.

        Parameters
        ----------
        only_most_recent :
            Set to true to filter out all but the most recent (based on date)
            for all given curator/statement/evidence triples.

        Returns
        -------
        :
            A set of statement hashes (i.e., from the "pa_hash" field) that
            have been marked as correct by any curator, for any evidence
        """
        return {
            curation["pa_hash"]
            for curation in self.get_curation_cache(only_most_recent=only_most_recent)
            if curation["tag"] == "correct"
        }
