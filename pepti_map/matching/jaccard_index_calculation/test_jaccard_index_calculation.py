from typing import Set
from math import isclose
from datasketch import LeanMinHash, MinHash
import numpy as np

from pepti_map.matching.jaccard_index_calculation.exact_jaccard_calculator import (
    ExactJaccardCalculator,
)
from pepti_map.matching.jaccard_index_calculation.min_hash_calculator import (
    MinHashCalculator,
)
from pepti_map.matching.match_merger import NUM_BYTES_FOR_MIN_HASH_VALUES


def generate_min_hash_for_set(test_set: Set) -> LeanMinHash:
    test_min_hash = MinHash()
    test_min_hash.update_batch(
        [elem.to_bytes(NUM_BYTES_FOR_MIN_HASH_VALUES, "big") for elem in test_set]
    )
    return LeanMinHash(test_min_hash)


class TestExactJaccardCalculator:
    def test_get_jaccard_index_matrix(self):
        jaccard_indexes = np.array(
            [
                [1.0, 0.0, 1 / 3, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [1 / 3, 0.0, 1.0, 1 / 7],
                [0.0, 0.0, 1 / 7, 1.0],
            ]
        )
        exact_jaccard_calculator = ExactJaccardCalculator(jaccard_indexes)
        np.testing.assert_array_equal(
            exact_jaccard_calculator.get_jaccard_index_matrix(), jaccard_indexes
        )

    def test_get_single_value(self):
        jaccard_indexes = np.array(
            [
                [1.0, 0.0, 1 / 3, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [1 / 3, 0.0, 1.0, 1 / 7],
                [0.0, 0.0, 1 / 7, 1.0],
            ]
        )
        exact_jaccard_calculator = ExactJaccardCalculator(jaccard_indexes)
        assert exact_jaccard_calculator.get_jaccard_index(0, 2) == 1 / 3
        assert exact_jaccard_calculator.get_jaccard_index(1, 1) == 1.0
        assert exact_jaccard_calculator.get_jaccard_index(3, 0) == 0.0
        assert exact_jaccard_calculator.get_jaccard_index(2, 3) == 1 / 7


class TestMinHashCalculator:
    def test_get_jaccard_index_matrix(self):
        test_matches = [
            {1, 2, 3, 4},
            {20, 30},
            {3, 4, 7, 8},
            {8, 10, 11, 12},
        ]
        test_min_hashes = [generate_min_hash_for_set(match) for match in test_matches]
        expected_jaccard_indexes = np.array(
            [
                [1.0, 0.0, 1 / 3, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [1 / 3, 0.0, 1.0, 1 / 7],
                [0.0, 0.0, 1 / 7, 1.0],
            ],
            dtype=np.float16,
        )
        min_hash_calculator = MinHashCalculator(test_min_hashes)
        np.testing.assert_allclose(
            min_hash_calculator.get_jaccard_index_matrix(),
            expected_jaccard_indexes,
            atol=0.015,
            rtol=0.0,
        )

    def test_get_single_value(self):
        test_matches = [
            {1, 2, 3, 4},
            {20, 30},
            {3, 4, 7, 8},
            {8, 10, 11, 12},
        ]
        test_min_hashes = [generate_min_hash_for_set(match) for match in test_matches]
        min_hash_calculator = MinHashCalculator(test_min_hashes)
        assert isclose(
            min_hash_calculator.get_jaccard_index(0, 2),
            1 / 3,
            rel_tol=0.0,
            abs_tol=0.015,
        )
        assert isclose(
            min_hash_calculator.get_jaccard_index(1, 1),
            1.0,
            rel_tol=0.0,
            abs_tol=0.015,
        )
        assert isclose(
            min_hash_calculator.get_jaccard_index(3, 0),
            0.0,
            rel_tol=0.0,
            abs_tol=0.015,
        )
        assert isclose(
            min_hash_calculator.get_jaccard_index(2, 3),
            1 / 7,
            rel_tol=0.0,
            abs_tol=0.015,
        )
