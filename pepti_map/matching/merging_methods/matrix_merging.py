import logging
from typing import Any, List, Set, Tuple, Union

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
        use_standard_numpy_array=True,
    ):
        super(MatrixMergingMethod, self).__init__(min_hashes, jaccard_index_threshold)
        self._merge_indication_matrix: Union[npt.NDArray[np.bool_], SymmetricMatrix]
        if use_standard_numpy_array:
            self._merge_indication_matrix = np.empty(
                shape=(len(min_hashes), len(min_hashes)), dtype=bool
            )
        else:
            self._merge_indication_matrix = SymmetricMatrix(
                size=len(min_hashes), dtype=bool
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
            return self._merge_indication_matrix[index_to_compare][  # type: ignore
                current_index
            ]
        return (
            current_min_hash.jaccard(min_hash_to_compare) > self.jaccard_index_threshold
        )

    def _fill_matrix_row_np(self, current_index: int) -> None:
        current_min_hash = self.min_hashes[current_index]
        self._merge_indication_matrix[current_index] = np.array(  # type: ignore
            [
                self._should_merge_sets(
                    current_index, current_min_hash, min_hash_index, min_hash
                )
                for min_hash_index, min_hash in enumerate(self.min_hashes)
            ],
            dtype=np.float16,
        )

    def _fill_matrix_row_symmetric(self, current_index: int) -> None:
        self._merge_indication_matrix.set_entry(  # type: ignore
            current_index, current_index, True
        )  # TODO: Is this line needed?
        current_min_hash = self.min_hashes[current_index]
        for index_to_compare in range(current_index + 1, len(self.min_hashes)):
            if (
                current_min_hash.jaccard(self.min_hashes[index_to_compare])
                > self.jaccard_index_threshold
            ):
                self._merge_indication_matrix.set_entry(  # type: ignore
                    current_index, index_to_compare, True
                )
                continue
            self._merge_indication_matrix.set_entry(  # type: ignore
                current_index, index_to_compare, False
            )

    def _fill_merge_indication_matrix_np(self) -> None:
        for current_index in range(0, len(self.min_hashes)):
            self._fill_matrix_row_np(current_index)

    def _fill_merge_indication_matrix_symmetric(self) -> None:
        for current_index in range(0, len(self.min_hashes)):
            self._fill_matrix_row_symmetric(current_index)

    def _generate_result_from_matrix(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        # TODO
        raise NotImplementedError

    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        # TODO: Not super beautiful
        # Is there a better implementation option?
        # Are both variants even needed?
        if isinstance(self._merge_indication_matrix, np.ndarray):
            self._fill_merge_indication_matrix_np()
        elif isinstance(self._merge_indication_matrix, SymmetricMatrix):
            self._fill_merge_indication_matrix_symmetric()
        else:
            assertion_error = (
                "Merge indication matrix is neither a numpy.ndarray "
                "nor a SymmetricMatrix instance."
            )
            logging.error(assertion_error)
            raise AssertionError(assertion_error)

        return self._generate_result_from_matrix(peptide_indexes, matches)
