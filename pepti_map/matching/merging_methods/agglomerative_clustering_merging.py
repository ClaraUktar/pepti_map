import logging
from typing import List, Set, Tuple

import numpy as np
import numpy.typing as npt
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash
from sklearn.cluster import AgglomerativeClustering

# TODO: Is using numpy.float16 sufficient?


class AgglomerativeClusteringMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        super(AgglomerativeClusteringMergingMethod, self).__init__(
            min_hashes, jaccard_index_threshold
        )
        self.distance_matrix: npt.NDArray[np.float16] = np.empty(
            shape=(len(min_hashes), len(min_hashes)), dtype=np.float16
        )

    def _get_jaccard_index_approximation(
        self,
        current_index: int,
        current_min_hash: LeanMinHash,
        index_to_compare: int,
        min_hash_to_compare: LeanMinHash,
    ) -> float:
        if current_index == index_to_compare:
            return 0.0
        if index_to_compare < current_index:
            return self.distance_matrix[index_to_compare][current_index]
        return current_min_hash.jaccard(min_hash_to_compare)

    def _generate_row_for_index(self, current_index: int) -> None:
        current_min_hash = self.min_hashes[current_index]
        self.distance_matrix[current_index] = np.array(
            [
                self._get_jaccard_index_approximation(
                    current_index, current_min_hash, min_hash_index, min_hash
                )
                for min_hash_index, min_hash in enumerate(self.min_hashes)
            ],
            dtype=np.float16,
        )

    def _generate_distance_matrix(self) -> None:
        for current_index in range(0, len(self.min_hashes)):
            self._generate_row_for_index(current_index)

    def _generate_result_from_labels(
        self,
        n_clusters: int,
        cluster_labels: npt.NDArray,
        peptide_indexes: List[int],
        matches: List[Set[int]],
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_matches = [set() for _ in range(0, n_clusters)]
        peptide_mappings = [[] for _ in range(0, n_clusters)]
        for label_index, cluster_label in enumerate(cluster_labels):
            merged_matches[cluster_label].update(matches[label_index])
            peptide_mappings[cluster_label].append(peptide_indexes[label_index])
        return merged_matches, peptide_mappings

    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        if len(peptide_indexes) != len(matches):
            error_message = (
                "The length of the given peptide_indexes and matches were expected "
                "to be the same, but differed (len(peptide_indexes) == "
                f"{len(peptide_indexes)}, len(matches) == {len(matches)})."
            )
            logging.error(error_message)
            raise ValueError(error_message)
        if len(peptide_indexes) < 2:
            peptide_mappings = []
            for peptide in peptide_indexes:
                peptide_mappings.append([peptide])
            return (matches, peptide_mappings)

        self._generate_distance_matrix()
        clustering = AgglomerativeClustering(
            metric="precomputed",
            compute_full_tree=True,
            n_clusters=None,
            distance_threshold=(1.0 - self.jaccard_index_threshold),
            linkage="single",
        )
        clustering.fit(self.distance_matrix)
        return self._generate_result_from_labels(
            clustering.n_clusters_, clustering.labels_, peptide_indexes, matches
        )
