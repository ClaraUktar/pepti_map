import unittest
from pepti_map.matching.match_merger_alternative import MatchMerger

# TODO: Use new match merger


class TestMatchMerger(unittest.TestCase):
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
        expected_result = (
            [{1, 2, 3, 4, 5}, {5, 6, 7}, {7, 8, 9, 10}, {9, 2}],
            [[0, 2], [1], [3, 5], [4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_matches_with_empty_sets(self):
        merge_results = MatchMerger(
            [
                {1, 2, 3},
                None,
                {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
                {6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21},
                None,
            ]
        ).merge_matches()
        expected_result = (
            [{1, 2, 3}, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21}],
            [[0], [2, 3]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_low_jaccard_index(self):
        merge_results = MatchMerger(
            [{1, 2, 3}, {3, 4, 5, 6}, {5, 6, 7, 8}, {8, 9, 10, 11}, {7, 8, 10, 12}],
            jaccard_index_threshold=0.3,
        ).merge_matches()
        expected_result = (
            [{1, 2, 3}, {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}],
            [[0], [1, 2, 3, 4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_late_merge_detection(self):
        merge_results = MatchMerger(
            [
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 20},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
            jaccard_index_threshold=0.1,
        ).merge_matches()
        expected_result = (
            [{17, 18, 19, 20}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}],
            [[3], [0, 1, 2, 4, 5]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_multiple_late_merge_detections(self):
        merge_results = MatchMerger(
            [
                {19, 20, 21, 22},
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 10},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
            jaccard_index_threshold=0.1,
        ).merge_matches()
        expected_result = (
            [
                {
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10,
                    11,
                    12,
                    13,
                    14,
                    15,
                    16,
                    17,
                    18,
                    19,
                    20,
                    21,
                    22,
                }
            ],
            [[0, 1, 2, 3, 4, 5, 6]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    # TODO: Add more tests
    # e.g. for matches list with None in it, different jaccard index
