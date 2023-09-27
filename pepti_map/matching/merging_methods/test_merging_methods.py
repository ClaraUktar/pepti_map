import unittest
import numpy as np

from pepti_map.matching.jaccard_index_calculation.exact_jaccard_calculator import (
    ExactJaccardCalculator,
)
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)
from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.full_matrix_merging import (
    FullMatrixMergingMethod,
)

# TODO: Refactor to reduce code duplication between test classes


def get_jaccard_array_int_representation(jaccard_array):
    return np.round(
        jaccard_array * IJaccardIndexCalculator.JACCARD_INT_MULTIPLICATION_FACTOR
    ).astype(np.uint16)


class TestAgglomerativeClusteringMergingMethod(unittest.TestCase):
    def test_empty_matches_return_empty_result(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([]))
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator
        ).generate_merged_result([], [])
        assert merge_results == ([], [])

    def test_single_match(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([[1.0]]))
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator
        ).generate_merged_result([0], [{11, 12, 13, 14}])
        assert merge_results == ([{11, 12, 13, 14}], [[0]])

    def test_floating_point(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([[1.0, 0.7], [0.7, 1.0]]))
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator
        ).generate_merged_result(
            [0, 1], [{1, 2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 8, 9, 10}]
        )
        expected_result = (
            [{1, 2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 8, 9, 10}],
            [[0], [1]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_multiple_matches(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 0.8, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 6],
                        [0.8, 1 / 7, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0, 0.25, 0.75],
                        [0.0, 0.0, 0.0, 0.25, 1.0, 0.2],
                        [0.0, 1 / 6, 0.0, 0.75, 0.2, 1.0],
                    ]
                )
            )
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5],
            [
                {1, 2, 3, 4},
                {5, 6, 7},
                {1, 2, 3, 4, 5},
                {8, 9, 10},
                {9, 2},
                {7, 8, 9, 10},
            ],
        )
        expected_result = (
            [{1, 2, 3, 4, 5}, {5, 6, 7}, {7, 8, 9, 10}, {9, 2}],
            [[0, 2], [1], [3, 5], [4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_low_jaccard_index(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 1 / 6, 0.0, 0.0, 0.0],
                        [1 / 6, 1.0, 1 / 3, 0.0, 0.0],
                        [0.0, 1 / 3, 1.0, 1 / 7, 1 / 3],
                        [0.0, 0.0, 1 / 7, 1.0, 1 / 3],
                        [0.0, 0.0, 1 / 3, 1 / 3, 1.0],
                    ]
                )
            )
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.3
        ).generate_merged_result(
            [0, 1, 2, 3, 4],
            [{1, 2, 3}, {3, 4, 5, 6}, {5, 6, 7, 8}, {8, 9, 10, 11}, {7, 8, 10, 12}],
        )
        expected_result = (
            [{1, 2, 3}, {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}],
            [[0], [1, 2, 3, 4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_late_merge_detection(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 1 / 7, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 7],
                        [1 / 7, 1 / 7, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 1.0, 1 / 7],
                        [0.0, 1 / 7, 0.0, 0.0, 1 / 7, 1.0],
                    ]
                )
            )
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.1
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5],
            [
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 20},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
        )
        expected_result = (
            [{17, 18, 19, 20}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}],
            [[3], [0, 1, 2, 4, 5]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_multiple_late_merge_detections(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 0.0, 0.0, 1 / 7, 0.0, 0.0],
                        [0.0, 1.0, 0.0, 1 / 7, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 7],
                        [0.0, 1 / 7, 1 / 7, 1.0, 1 / 7, 0.0, 0.0],
                        [1 / 7, 0.0, 0.0, 1 / 7, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1 / 7],
                        [0.0, 0.0, 1 / 7, 0.0, 0.0, 1 / 7, 1.0],
                    ]
                )
            )
        )
        merge_results = AgglomerativeClusteringMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.1
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5, 6],
            [
                {19, 20, 21, 22},
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 10},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
        )
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


class TestFullMatrixMergingMethod(unittest.TestCase):
    def test_empty_matches_return_empty_result(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([]))
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator
        ).generate_merged_result([], [])
        assert merge_results == ([], [])

    def test_single_match(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([[1.0]]))
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator
        ).generate_merged_result([0], [{11, 12, 13, 14}])
        assert merge_results == ([{11, 12, 13, 14}], [[0]])

    def test_floating_point(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(np.array([[1.0, 0.7], [0.7, 1.0]]))
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator
        ).generate_merged_result(
            [0, 1], [{1, 2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 8, 9, 10}]
        )
        expected_result = (
            [{1, 2, 3, 4, 5, 6, 7, 8}, {2, 3, 4, 5, 6, 7, 8, 9, 10}],
            [[0], [1]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_multiple_matches(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 0.8, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 6],
                        [0.8, 1 / 7, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0, 0.25, 0.75],
                        [0.0, 0.0, 0.0, 0.25, 1.0, 0.2],
                        [0.0, 1 / 6, 0.0, 0.75, 0.2, 1.0],
                    ]
                )
            )
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5],
            [
                {1, 2, 3, 4},
                {5, 6, 7},
                {1, 2, 3, 4, 5},
                {8, 9, 10},
                {9, 2},
                {7, 8, 9, 10},
            ],
        )
        expected_result = (
            [{1, 2, 3, 4, 5}, {5, 6, 7}, {7, 8, 9, 10}, {9, 2}],
            [[0, 2], [1], [3, 5], [4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_low_jaccard_index(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 1 / 6, 0.0, 0.0, 0.0],
                        [1 / 6, 1.0, 1 / 3, 0.0, 0.0],
                        [0.0, 1 / 3, 1.0, 1 / 7, 1 / 3],
                        [0.0, 0.0, 1 / 7, 1.0, 1 / 3],
                        [0.0, 0.0, 1 / 3, 1 / 3, 1.0],
                    ]
                )
            )
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.3
        ).generate_merged_result(
            [0, 1, 2, 3, 4],
            [{1, 2, 3}, {3, 4, 5, 6}, {5, 6, 7, 8}, {8, 9, 10, 11}, {7, 8, 10, 12}],
        )
        expected_result = (
            [{1, 2, 3}, {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}],
            [[0], [1, 2, 3, 4]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_late_merge_detection(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 1 / 7, 0.0, 0.0, 0.0],
                        [0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 7],
                        [1 / 7, 1 / 7, 1.0, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 1.0, 1 / 7],
                        [0.0, 1 / 7, 0.0, 0.0, 1 / 7, 1.0],
                    ]
                )
            )
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.1
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5],
            [
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 20},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
        )
        expected_result = (
            [{17, 18, 19, 20}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}],
            [[3], [0, 1, 2, 4, 5]],
        )
        self.assertCountEqual(
            zip(merge_results[0], merge_results[1]),
            zip(expected_result[0], expected_result[1]),
        )

    def test_multiple_late_merge_detections(self):
        jaccard_calculator = ExactJaccardCalculator(
            get_jaccard_array_int_representation(
                np.array(
                    [
                        [1.0, 0.0, 0.0, 0.0, 1 / 7, 0.0, 0.0],
                        [0.0, 1.0, 0.0, 1 / 7, 0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0, 1 / 7, 0.0, 0.0, 1 / 7],
                        [0.0, 1 / 7, 1 / 7, 1.0, 1 / 7, 0.0, 0.0],
                        [1 / 7, 0.0, 0.0, 1 / 7, 1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1 / 7],
                        [0.0, 0.0, 1 / 7, 0.0, 0.0, 1 / 7, 1.0],
                    ]
                )
            )
        )
        merge_results = FullMatrixMergingMethod(
            jaccard_calculator, jaccard_index_threshold=0.1
        ).generate_merged_result(
            [0, 1, 2, 3, 4, 5, 6],
            [
                {19, 20, 21, 22},
                {1, 2, 3, 4},
                {5, 6, 7, 8},
                {1, 5, 9, 10},
                {17, 18, 19, 10},
                {11, 12, 13, 14},
                {8, 11, 15, 16},
            ],
        )
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
