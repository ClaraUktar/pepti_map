from unittest.mock import MagicMock
import numpy as np

from pepti_map.matching.match_merger import MatchMerger


class TestMatchMerger:
    def test_filter_out_none_entries(self):
        test_matches = [{1, 2, 3, 4}, None, None, {5, 6, 7, 8}, {9, 10, 11, 12}, None]
        filtered_matches = [{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}]
        match_merger = MatchMerger(test_matches)
        assert match_merger.matches == filtered_matches

    def test_filter_out_none_entries_precomputed_intersections(self):
        test_matches = [
            {1, 2, 3, 4},
            {20, 30},
            None,
            {3, 4, 7, 8},
            {8, 10, 11, 12},
            None,
        ]
        filtered_matches = [
            {1, 2, 3, 4},
            {20, 30},
            {3, 4, 7, 8},
            {8, 10, 11, 12},
        ]
        test_precomputed_intersections = np.array(
            [
                [4, 0, 0, 2, 0, 0],
                [0, 2, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [2, 0, 0, 4, 1, 0],
                [0, 0, 0, 1, 4, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        )
        filtered_intersections = np.array(
            [
                [4, 0, 2, 0],
                [0, 2, 0, 0],
                [2, 0, 4, 1],
                [0, 0, 1, 4],
            ]
        )
        MatchMerger._create_exact_jaccard_calculator = MagicMock()
        match_merger = MatchMerger(test_matches, 0.7, test_precomputed_intersections)
        assert match_merger.matches == filtered_matches
        np.testing.assert_array_equal(
            MatchMerger._create_exact_jaccard_calculator.call_args.args[0],
            filtered_intersections,
        )
