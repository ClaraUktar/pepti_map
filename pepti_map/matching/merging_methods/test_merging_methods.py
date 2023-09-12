from typing import Set
import unittest

from datasketch import LeanMinHash, MinHash
from pepti_map.matching.match_merger import NUM_BYTES_FOR_MIN_HASH_VALUES

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)


def generate_min_hash_for_set(test_set: Set) -> LeanMinHash:
    test_min_hash = MinHash()
    test_min_hash.update_batch(
        [elem.to_bytes(NUM_BYTES_FOR_MIN_HASH_VALUES, "big") for elem in test_set]
    )
    return LeanMinHash(test_min_hash)


class TestAgglomerativeClustering(unittest.TestCase):
    # TODO
    pass
    # def test_empty_matches_return_empty_result(self):
    #     merge_results = AgglomerativeClusteringMergingMethod([]).generate_merged_result(
    #         [], []
    #     )
    #     assert merge_results == ([], [])

    # def test_single_match(self):
    #     test_set = {11, 12, 13, 14}
    #     test_min_hash = generate_min_hash_for_set(test_set)
    #     merge_result = AgglomerativeClusteringMergingMethod(
    #         [test_min_hash]
    #     ).generate_merged_result([0], [test_set])
    #     assert merge_result == ([{11, 12, 13, 14}], [[0]])

    # def test_multiple_matches(self):
    #     test_sets = [
    #         {1, 2, 3, 4},
    #         {5, 6, 7},
    #         {1, 2, 3, 4, 5},
    #         {8, 9, 10},
    #         {9, 2},
    #         {7, 8, 9, 10},
    #     ]
    #     merge_results = AgglomerativeClusteringMergingMethod(
    #         [generate_min_hash_for_set(test_set) for test_set in test_sets], 1.0
    #     ).generate_merged_result([i for i in range(0, len(test_sets))], test_sets)
    #     expected_result = (
    #         [{1, 2, 3, 4, 5}, {5, 6, 7}, {7, 8, 9, 10}, {9, 2}],
    #         [[0, 2], [1], [3, 5], [4]],
    #     )
    #     self.assertCountEqual(
    #         zip(merge_results[0], merge_results[1]),
    #         zip(expected_result[0], expected_result[1]),
    #     )

    # def test_matches_with_empty_sets(self):
    #     merge_results = MatchMerger(
    #         [
    #             {1, 2, 3},
    #             None,
    #             {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17},
    #             {6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21},
    #             None,
    #         ]
    #     ).merge_matches()
    #     expected_result = (
    #         [{1, 2, 3}, {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20, 21}],
    #         [[0], [2, 3]],
    #     )
    #     self.assertCountEqual(
    #         zip(merge_results[0], merge_results[1]),
    #         zip(expected_result[0], expected_result[1]),
    #     )

    # def test_low_jaccard_index(self):
    #     merge_results = MatchMerger(
    #         [{1, 2, 3}, {3, 4, 5, 6}, {5, 6, 7, 8}, {8, 9, 10, 11}, {7, 8, 10, 12}],
    #         jaccard_index_threshold=0.3,
    #     ).merge_matches()
    #     expected_result = (
    #         [{1, 2, 3}, {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}],
    #         [[0], [1, 2, 3, 4]],
    #     )
    #     self.assertCountEqual(
    #         zip(merge_results[0], merge_results[1]),
    #         zip(expected_result[0], expected_result[1]),
    #     )

    # def test_late_merge_detection(self):
    #     merge_results = MatchMerger(
    #         [
    #             {1, 2, 3, 4},
    #             {5, 6, 7, 8},
    #             {1, 5, 9, 10},
    #             {17, 18, 19, 20},
    #             {11, 12, 13, 14},
    #             {8, 11, 15, 16},
    #         ],
    #         jaccard_index_threshold=0.1,
    #     ).merge_matches()
    #     expected_result = (
    #         [{17, 18, 19, 20}, {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16}],
    #         [[3], [0, 1, 2, 4, 5]],
    #     )
    #     self.assertCountEqual(
    #         zip(merge_results[0], merge_results[1]),
    #         zip(expected_result[0], expected_result[1]),
    #     )

    # def test_multiple_late_merge_detections(self):
    #     merge_results = MatchMerger(
    #         [
    #             {19, 20, 21, 22},
    #             {1, 2, 3, 4},
    #             {5, 6, 7, 8},
    #             {1, 5, 9, 10},
    #             {17, 18, 19, 10},
    #             {11, 12, 13, 14},
    #             {8, 11, 15, 16},
    #         ],
    #         jaccard_index_threshold=0.1,
    #     ).merge_matches()
    #     expected_result = (
    #         [
    #             {
    #                 1,
    #                 2,
    #                 3,
    #                 4,
    #                 5,
    #                 6,
    #                 7,
    #                 8,
    #                 9,
    #                 10,
    #                 11,
    #                 12,
    #                 13,
    #                 14,
    #                 15,
    #                 16,
    #                 17,
    #                 18,
    #                 19,
    #                 20,
    #                 21,
    #                 22,
    #             }
    #         ],
    #         [[0, 1, 2, 3, 4, 5, 6]],
    #     )
    #     self.assertCountEqual(
    #         zip(merge_results[0], merge_results[1]),
    #         zip(expected_result[0], expected_result[1]),
    #     )
