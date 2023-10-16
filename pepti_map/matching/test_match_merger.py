from typing import Set
from unittest.mock import patch
from datasketch import LeanMinHash, MinHash
import numpy as np
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)

from pepti_map.matching.match_merger import NUM_BYTES_FOR_MIN_HASH_VALUES, MatchMerger


def generate_min_hash_for_set(test_set: Set) -> LeanMinHash:
    test_min_hash = MinHash()
    test_min_hash.update_batch(
        [elem.to_bytes(NUM_BYTES_FOR_MIN_HASH_VALUES, "big") for elem in test_set]
    )
    return LeanMinHash(test_min_hash)


class TestMatchMerger:
    def test_filter_out_none_entries(self):
        test_matches = [{1, 2, 3, 4}, None, None, {5, 6, 7, 8}, {9, 10, 11, 12}, None]
        filtered_matches = [{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}]
        match_merger = MatchMerger(test_matches)
        assert match_merger.matches == filtered_matches

    @patch(
        "pepti_map.matching.match_merger.MatchMerger._create_exact_jaccard_calculator"
    )
    def test_filter_out_none_entries_precomputed_intersections(
        self, mocked_create_exact_jaccard_calculator
    ):
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
        match_merger = MatchMerger(test_matches, 0.7, test_precomputed_intersections)
        assert match_merger.matches == filtered_matches
        np.testing.assert_array_equal(
            mocked_create_exact_jaccard_calculator.call_args.args[0],
            filtered_intersections,
        )

    @patch(
        (
            "pepti_map.matching.jaccard_index_calculation.min_hash_calculator"
            ".MinHashCalculator.__init__"
        ),
        return_value=None,
    )
    def test_min_hash_creation(self, mock_min_hash_calculator_init):
        test_matches = [
            {1, 2, 3, 4},
            {20, 30},
            None,
            {3, 4, 7, 8},
            {8, 10, 11, 12},
            None,
        ]
        expected_min_hashes = [
            generate_min_hash_for_set(match)
            for match in test_matches
            if match is not None
        ]
        MatchMerger(test_matches, 0.7)
        mock_min_hash_calculator_init.assert_called_once()
        assert mock_min_hash_calculator_init.call_args.args[0] == expected_min_hashes

    @patch(
        (
            "pepti_map.matching.jaccard_index_calculation.exact_jaccard_calculator"
            ".ExactJaccardCalculator.__init__"
        ),
        return_value=None,
    )
    def test_exact_jaccard_calculation(self, mock_exact_jaccard_calculator_init):
        test_matches = [
            {1, 2, 3, 4},
            {20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
            None,
            {3, 4, 7, 8},
            {8, 10, 11, 12},
            None,
        ]
        test_precomputed_intersections = np.array(
            [
                [4, 0, 0, 2, 0, 0],
                [0, 10, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [2, 0, 0, 4, 1, 0],
                [0, 0, 0, 1, 4, 0],
                [0, 0, 0, 0, 0, 0],
            ],
            dtype=np.uint32,
        )
        expected_jaccard_indexes = np.round(
            np.array(
                [
                    [1.0, 0.0, 1 / 3, 0.0],
                    [0.0, 1.0, 0.0, 0.0],
                    [1 / 3, 0.0, 1.0, 1 / 7],
                    [0.0, 0.0, 1 / 7, 1.0],
                ],
            )
            * IJaccardIndexCalculator.JACCARD_INT_MULTIPLICATION_FACTOR
        ).astype(np.uint16)
        MatchMerger(test_matches, 0.7, test_precomputed_intersections)
        mock_exact_jaccard_calculator_init.assert_called_once()
        np.testing.assert_allclose(
            mock_exact_jaccard_calculator_init.call_args.args[0],
            expected_jaccard_indexes,
            atol=1,
        )
