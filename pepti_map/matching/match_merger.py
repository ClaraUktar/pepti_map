from typing import List, Literal, Set, Tuple, Union
from datasketch import LeanMinHash, MinHash
import numpy as np
import numpy.typing as npt
from pepti_map.matching.jaccard_index_calculation.exact_jaccard_calculator import (
    ExactJaccardCalculator,
)
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)
from pepti_map.matching.jaccard_index_calculation.min_hash_calculator import (
    MinHashCalculator,
)

from pepti_map.matching.merging_methods.merging_method_helper import get_merging_method

NUM_BYTES_FOR_MIN_HASH_VALUES = 4


class MatchMerger:
    def __init__(
        self,
        matches: List[Union[Set[int], None]],
        jaccard_index_threshold: float = 0.7,
        # TODO: Add array dtype?
        precomputed_intersections: Union[npt.NDArray, None] = None,
    ):
        # TODO: Want to delete matches after merging
        self.jaccard_index_threshold: float = jaccard_index_threshold
        self.peptide_indexes: List[int] = []
        self.matches: List[Set[int]] = []
        deleted_indexes = self._postprocess_matches(matches)
        self.jaccard_calculator: IJaccardIndexCalculator
        if precomputed_intersections is not None:
            precomputed_intersections = self._postprocess_precomputed_intersections(
                precomputed_intersections, deleted_indexes
            )
            self._create_exact_jaccard_calculator(precomputed_intersections)
        else:
            self._create_min_hash_calculator()

    def _postprocess_matches(
        self, uncleaned_matches: List[Union[Set[int], None]]
    ) -> List[int]:
        deleted_indexes = []
        for peptide_index, uncleaned_match in enumerate(uncleaned_matches):
            if uncleaned_match is None:
                deleted_indexes.append(peptide_index)
                continue
            self.matches.append(uncleaned_match)
            self.peptide_indexes.append(peptide_index)
        return deleted_indexes

    def _postprocess_precomputed_intersections(
        self, precomputed_intersections: npt.NDArray, deleted_indexes: List[int]
    ) -> npt.NDArray:
        return np.delete(
            np.delete(precomputed_intersections, deleted_indexes, axis=0),
            deleted_indexes,
            axis=1,
        )

    def _create_exact_jaccard_calculator(
        self, precomputed_intersections: npt.NDArray
    ) -> None:
        set_sizes = np.array([len(match) for match in self.matches], dtype=np.uint32)
        row_index, column_index = np.ix_(
            np.arange(precomputed_intersections.shape[0]),
            np.arange(precomputed_intersections.shape[1]),
        )
        # Broadcasted cell-wise Jaccard-Index computation
        precomputed_intersections = precomputed_intersections / (
            set_sizes[row_index] + set_sizes[column_index] - precomputed_intersections
        )
        self.jaccard_calculator = ExactJaccardCalculator(precomputed_intersections)

    def _create_min_hash_calculator(self) -> None:
        min_hashes: List[LeanMinHash] = []
        for match in self.matches:
            min_hash = MinHash()
            min_hash.update_batch(
                [
                    element_to_add.to_bytes(NUM_BYTES_FOR_MIN_HASH_VALUES, "big")
                    for element_to_add in match
                ]
            )
            min_hashes.append(LeanMinHash(min_hash))

        self.jaccard_calculator = MinHashCalculator(min_hashes)

    def merge_matches(
        self,
        method: Literal["agglomerative-clustering", "full-matrix"],
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        # TODO: Add options to parameterize methods?
        return get_merging_method(
            method, self.jaccard_calculator, self.jaccard_index_threshold
        ).generate_merged_result(self.peptide_indexes, self.matches)
