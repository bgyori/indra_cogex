import unittest

from indra_cogex.apps.curator.curation_cache import CurationCache, Curations


class MockCurationCache(CurationCache):
    def __init__(self, curations: Curations):
        super().__init__()
        self.curation_list.extend(
            sorted(
                (self._process_curation(curation) for curation in curations),
                key=self._curation_key,
            )
        )

    def refresh_curations(self):
        pass


TEST_DATE = "Thu, 29 Nov 2018 18:00:08 GMT"
LATER_TEST_DATE = "Thu, 29 Nov 2018 19:00:08 GMT"


def _curation(date=TEST_DATE, curator="charlie", **kwargs):
    return dict(curator=curator, date=date, **kwargs)


class TestCurationCache(unittest.TestCase):
    def test_get_incorrect_evidence_hashes(self):
        curations = [
            _curation(pa_hash=1, source_hash=1, tag="correct"),
            _curation(pa_hash=1, source_hash=1, tag="incorrect"),
            _curation(pa_hash=1, source_hash=2, tag="correct"),
            _curation(pa_hash=1, source_hash=3, tag="incorrect"),
        ]
        curation_cache = MockCurationCache(curations)
        self.assertEqual({3}, curation_cache.get_incorrect_evidence_hashes())
        self.assertEqual({1, 2}, curation_cache.get_correct_evidence_hashes())
        self.assertEqual({1, 2, 3}, curation_cache.get_curated_evidence_hashes())

    def test_get_correct_statement_hashes(self):
        curations = [
            _curation(pa_hash=1, source_hash=1, tag="correct"),
            _curation(pa_hash=1, source_hash=1, tag="incorrect"),
            _curation(pa_hash=2, source_hash=2, tag="correct"),
            _curation(pa_hash=3, source_hash=3, tag="incorrect"),
        ]
        curation_cache = MockCurationCache(curations)
        statement_hashes = curation_cache.get_correct_statement_hashes()
        self.assertEqual({1, 2}, statement_hashes)

    def test_get_recent_curations(self):
        # this simulates the scenario when a curation is later amended
        input_curations = [
            dict(source_hash=1, pa_hash=1, curator="ben", tag="nope", date=TEST_DATE),
            dict(
                source_hash=1,
                pa_hash=1,
                curator="charlie",
                tag="correct",
                date=TEST_DATE,
            ),
            dict(
                source_hash=1,
                pa_hash=1,
                curator="charlie",
                tag="nope",
                date=LATER_TEST_DATE,
            ),
            dict(
                source_hash=1,
                pa_hash=2,
                curator="charlie",
                tag="correct",
                date=TEST_DATE,
            ),
        ]
        curation_cache = MockCurationCache(input_curations)
        expected = [
            dict(source_hash=1, pa_hash=1, curator="ben", tag="nope", date=TEST_DATE),
            dict(
                source_hash=1,
                pa_hash=1,
                curator="charlie",
                tag="nope",
                date=LATER_TEST_DATE,
            ),
            dict(
                source_hash=1,
                pa_hash=2,
                curator="charlie",
                tag="correct",
                date=TEST_DATE,
            ),
        ]
        expected = [curation_cache._process_curation(curation) for curation in expected]
        self.assertEqual(
            expected, curation_cache.get_curation_cache(only_most_recent=True)
        )
