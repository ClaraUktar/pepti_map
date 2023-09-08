from typing import List, Set, Tuple

import numpy as np
import numpy.typing as npt
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash


class FullMatrixMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        super(FullMatrixMergingMethod, self).__init__(
            min_hashes, jaccard_index_threshold
        )
        self._merge_indication_matrix: npt.NDArray[np.bool_] = np.empty(
            shape=(len(min_hashes), len(min_hashes)), dtype=bool
        )

    def _should_merge_sets(
        self,
        current_index: int,
        current_min_hash: LeanMinHash,
        index_to_compare: int,
        min_hash_to_compare: LeanMinHash,
    ) -> bool:
        if current_index == index_to_compare:
            return True
        if index_to_compare < current_index:
            return self._merge_indication_matrix[index_to_compare][current_index]
        return (
            current_min_hash.jaccard(min_hash_to_compare) > self.jaccard_index_threshold
        )

    def _fill_matrix_row(self, current_index: int) -> None:
        current_min_hash = self.min_hashes[current_index]
        self._merge_indication_matrix[current_index] = np.array(
            [
                self._should_merge_sets(
                    current_index, current_min_hash, min_hash_index, min_hash
                )
                for min_hash_index, min_hash in enumerate(self.min_hashes)
            ],
            dtype=bool,
        )

    def _fill_merge_indication_matrix(self) -> None:
        for current_index in range(0, len(self.min_hashes)):
            self._fill_matrix_row(current_index)

    def _generate_result_from_matrix(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        # TODO
        raise NotImplementedError

    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        self._fill_merge_indication_matrix()
        return self._generate_result_from_matrix(peptide_indexes, matches)
