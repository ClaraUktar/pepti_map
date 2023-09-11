from typing import List
from datasketch import LeanMinHash
import numpy as np
import numpy.typing as npt

from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)


class MinHashCalculator(IJaccardIndexCalculator):
    def __init__(self, min_hashes: List[LeanMinHash]):
        self.min_hashes: List[LeanMinHash] = min_hashes

    def get_jaccard_index(self, first_index: int, second_index: int) -> float:
        return self.min_hashes[first_index].jaccard(self.min_hashes[second_index])

    def _get_matrix_value(
        self,
        current_index: int,
        current_min_hash: LeanMinHash,
        index_to_compare: int,
        min_hash_to_compare: LeanMinHash,
        jaccard_index_matrix: npt.NDArray,
    ) -> float:
        if current_index == index_to_compare:
            return 1.0
        if index_to_compare < current_index:
            return jaccard_index_matrix[index_to_compare][current_index]
        return current_min_hash.jaccard(min_hash_to_compare)

    def get_jaccard_index_matrix(self) -> npt.NDArray[np.float16]:
        jaccard_index_matrix = np.empty(
            shape=(len(self.min_hashes), len(self.min_hashes)), dtype=np.float16
        )
        for row_index in range(0, len(self.min_hashes)):
            current_min_hash = self.min_hashes[row_index]
            jaccard_index_matrix[row_index] = np.array(
                [
                    self._get_matrix_value(
                        row_index,
                        current_min_hash,
                        min_hash_index,
                        min_hash,
                        jaccard_index_matrix,
                    )
                    for min_hash_index, min_hash in enumerate(self.min_hashes)
                ],
                dtype=np.float16,
            )
        return jaccard_index_matrix
