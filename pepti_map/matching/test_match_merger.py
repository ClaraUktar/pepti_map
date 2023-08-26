from pepti_map.matching.match_merger import MatchMerger


class TestMatchMerger:
    def test_empty_matches_return_empty_result(self):
        merge_result = MatchMerger([]).merge_matches()
        assert merge_result == ([], [])

    def test_single_match(self):
        merge_result = MatchMerger([{11, 12, 13, 14}]).merge_matches()
        assert merge_result == ([{11, 12, 13, 14}], [[0]])

    # TODO: Add more tests
