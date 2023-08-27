from pepti_map.matching.match_merger import MatchMerger


class TestMatchMerger:
    def test_empty_matches_return_empty_result(self):
        merge_result = MatchMerger([]).merge_matches()
        assert merge_result == ([], [])

    def test_single_match(self):
        merge_result = MatchMerger([{11, 12, 13, 14}]).merge_matches()
        assert merge_result == ([{11, 12, 13, 14}], [[0]])

    def test_multiple_matches(self):
        merge_results = MatchMerger(
            [
                {1, 2, 3, 4},
                {5, 6, 7},
                {1, 2, 3, 4, 5},
                {8, 9, 10},
                {9, 2},
                {7, 8, 9, 10},
            ]
        ).merge_matches()
        assert merge_results == (
            [{1, 2, 3, 4, 5}, {5, 6, 7}, {7, 8, 9, 10}, {9, 2}],
            [[0, 2], [1], [3, 5], [4]],
        )

    # TODO: Add more tests
    # e.g. for matches list with None in it, different jaccard index
