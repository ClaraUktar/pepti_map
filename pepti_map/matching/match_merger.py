from typing import List, Literal, Set, Tuple, Union
from datasketch import LeanMinHash, MinHash
import numpy.typing as npt
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
        precomputed_intersections: Union[npt.NDArray, None] = None,
    ):
        # TODO: Want to delete matches after merging
        self.jaccard_index_threshold: float = jaccard_index_threshold
        self.peptide_indexes: List[int] = []
        self.matches: List[Set[int]] = []
        self._postprocess_matches(matches)
        # TODO: Also postprocess precomputed intersections to remove entries
        # corresponding with the None entries
        self.jaccard_calculator: IJaccardIndexCalculator
        if precomputed_intersections is not None:
            self._create_exact_jaccard_calculator()
        else:
            self._create_min_hash_calculator()

    def _postprocess_matches(
        self, uncleaned_matches: List[Union[Set[int], None]]
    ) -> None:
        for peptide_index, uncleaned_match in enumerate(uncleaned_matches):
            if uncleaned_match is None:
                continue
            self.matches.append(uncleaned_match)
            self.peptide_indexes.append(peptide_index)

    def _create_exact_jaccard_calculator(self) -> None:
        # TODO
        pass

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
