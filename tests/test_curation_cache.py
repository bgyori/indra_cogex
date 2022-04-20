import unittest

from indra_cogex.apps.curator.curation_cache import CurationCache, Curations


class MockCurationCache(CurationCache):
    def __init__(self, curations: Curations):
        super().__init__()
        self.curation_list.extend(curations)

    def refresh_curations(self):
        pass


class TestCurationCache(unittest.TestCase):
    def test_get_incorrect_source_hashes(self):
        curations = [
            dict(source_hash=1, tag="correct"),
            dict(source_hash=1, tag="incorrect"),
            dict(source_hash=2, tag="correct"),
            dict(source_hash=3, tag="incorrect"),
        ]
        curation_cache = MockCurationCache(curations)
        incorrect_hashes = curation_cache.get_incorrect_source_hashes()
        self.assertEqual({3}, incorrect_hashes)
