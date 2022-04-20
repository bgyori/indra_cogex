import unittest

from indra_cogex.apps.curator.curation_cache import CurationCache, Curations


class MockCurationCache(CurationCache):
    def __init__(self, curations: Curations):
        super().__init__()
        self.curation_list.extend(curations)

    def refresh_curations(self):
        pass


class TestCurationCache(unittest.TestCase):
    def test_get_incorrect_evidence_hashes(self):
        curations = [
            dict(source_hash=1, tag="correct"),
            dict(source_hash=1, tag="incorrect"),
            dict(source_hash=2, tag="correct"),
            dict(source_hash=3, tag="incorrect"),
        ]
        curation_cache = MockCurationCache(curations)
        evidence_hashes = curation_cache.get_incorrect_evidence_hashes()
        self.assertEqual({3}, evidence_hashes)

    def test_get_correct_statement_hashes(self):
        curations = [
            dict(pa_hash=1, tag="correct"),
            dict(pa_hash=1, tag="incorrect"),
            dict(pa_hash=2, tag="correct"),
            dict(pa_hash=3, tag="incorrect"),
        ]
        curation_cache = MockCurationCache(curations)
        statement_hashes = curation_cache.get_correct_statement_hashes()
        self.assertEqual({1, 2}, statement_hashes)
