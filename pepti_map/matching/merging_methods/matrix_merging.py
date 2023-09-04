from typing import Any, List, Set, Tuple

import numpy as np
import numpy.typing as npt
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash


class SymmetricMatrix:
    def __init__(self, size: int, dtype: npt.DTypeLike):
        self.matrix = np.empty(
            shape=(SymmetricMatrix._gaussian_sum(size),), dtype=dtype
        )
        self.size = size

    def _get_true_index_for_cell(self, row: int, column: int) -> int:
        if row > column:
            switch = column
            column = row
            row = switch
        return self.size * row - SymmetricMatrix._gaussian_sum(row - 1) + column - row

    def get_entry(self, row: int, column: int) -> Any:
        if row >= self.size or column >= self.size:
            raise IndexError
        return self.matrix[self._get_true_index_for_cell(row, column)]

    def set_entry(self, row: int, column: int, entry: Any) -> None:
        if row >= self.size or column >= self.size:
            raise IndexError
        self.matrix[self._get_true_index_for_cell(row, column)] = entry

    @staticmethod
    def _gaussian_sum(n: int) -> int:
        return int((n * (n + 1)) / 2)


class MatrixMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        super(MatrixMergingMethod, self).__init__(min_hashes, jaccard_index_threshold)
        self.merge_indication_matrix = SymmetricMatrix(size=len(min_hashes), dtype=bool)

    def _fill_matrix_row(self, current_index: int) -> None:
        self.merge_indication_matrix.set_entry(
            current_index, current_index, True
        )  # TODO: Is this line needed?
        current_min_hash = self.min_hashes[current_index]
        for index_to_compare in range(current_index + 1, len(self.min_hashes)):
            if (
                current_min_hash.jaccard(self.min_hashes[index_to_compare])
                > self.jaccard_index_threshold
            ):
                self.merge_indication_matrix.set_entry(
                    current_index, index_to_compare, True
                )
                continue
            self.merge_indication_matrix.set_entry(
                current_index, index_to_compare, False
            )

    def _fill_merge_indication_matrix(self) -> None:
        for current_index in range(0, len(self.min_hashes)):
            self._fill_matrix_row(current_index)

    def _generate_result_from_matrix(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        raise NotImplementedError

    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        self._fill_merge_indication_matrix()
        return self._generate_result_from_matrix(peptide_indexes, matches)
