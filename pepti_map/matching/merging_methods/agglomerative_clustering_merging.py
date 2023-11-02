import logging
from typing import List, Set, Tuple

import numpy as np
import numpy.typing as npt
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from sklearn.cluster import AgglomerativeClustering

# TODO: Is using numpy.float16 sufficient?
# Already has some rounding problems at second decimal
# Better to use float32 or int (as first digit will always be 0 or 1)


class AgglomerativeClusteringMergingMethod(IMergingMethod):
    def __init__(
        self,
        jaccard_calculator: IJaccardIndexCalculator,
        jaccard_index_threshold: float = 0.5,
    ):
        self._distance_matrix: npt.NDArray[np.uint16]
        super(AgglomerativeClusteringMergingMethod, self).__init__(
            jaccard_calculator, jaccard_index_threshold
        )

    def _process_jaccard_calculator(
        self, jaccard_calculator: IJaccardIndexCalculator
    ) -> None:
        # TODO: Does this result in two matrices?
        self._distance_matrix = np.subtract(
            round(1.0 * IJaccardIndexCalculator.JACCARD_INT_MULTIPLICATION_FACTOR),
            jaccard_calculator.get_jaccard_index_matrix(),
            dtype=np.uint16,
        )

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

        clustering = AgglomerativeClustering(
            metric="precomputed",
            compute_full_tree=True,
            n_clusters=None,
            distance_threshold=(
                round(1.0 * IJaccardIndexCalculator.JACCARD_INT_MULTIPLICATION_FACTOR)
                - self.jaccard_index_threshold
            ),
            linkage="single",
        )
        clustering.fit(self._distance_matrix)
        return self._generate_result_from_labels(
            clustering.n_clusters_, clustering.labels_, peptide_indexes, matches
        )
