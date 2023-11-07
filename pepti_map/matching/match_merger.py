from typing import List, Literal, Set, Tuple, Union
from datasketch import LeanMinHash, MinHash
from pathlib import Path
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
from pepti_map.constants import PATH_TO_MERGED_INDEXES, PATH_TO_MERGED_MATCHES

NUM_BYTES_FOR_MIN_HASH_VALUES = 4


class MatchMerger:
    def __init__(
        self,
        matches: List[Union[Set[int], None]],
        jaccard_index_threshold: float = 0.5,
        precomputed_intersections: Union[npt.NDArray[np.uint32], None] = None,
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
            del precomputed_intersections
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
        self,
        precomputed_intersections: npt.NDArray[np.uint32],
        deleted_indexes: List[int],
    ) -> npt.NDArray[np.uint32]:
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
        precomputed_intersections_large = np.floor_divide(
            precomputed_intersections
            * ExactJaccardCalculator.JACCARD_INT_MULTIPLICATION_FACTOR,
            (
                set_sizes[row_index]
                + set_sizes[column_index]
                - precomputed_intersections
            ),
            dtype=np.uint32,
        )
        precomputed_intersections = precomputed_intersections_large.astype(np.uint16)
        del precomputed_intersections_large

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

    # TODO: Add test
    @staticmethod
    def save_merged_result(
        merged_matches: List[Set[int]],
        peptide_indexes: List[List[int]],
        path_to_sets_file: Path = PATH_TO_MERGED_MATCHES,
        path_to_indexes_file: Path = PATH_TO_MERGED_INDEXES,
    ) -> None:
        with open(path_to_sets_file, "wt", encoding="utf-8") as merged_matches_file:
            merged_matches_file.writelines(
                [
                    (",".join([str(match_elem) for match_elem in merged_match]) + "\n")
                    for merged_match in merged_matches
                ]
            )

        with open(path_to_indexes_file, "wt", encoding="utf-8") as peptide_indexes_file:
            peptide_indexes_file.writelines(
                [
                    (",".join([str(index_elem) for index_elem in peptide_index]) + "\n")
                    for peptide_index in peptide_indexes
                ]
            )

    # TODO: Add test
    @staticmethod
    def load_merged_result(
        path_to_sets_file: Path = PATH_TO_MERGED_MATCHES,
        path_to_indexes_file: Path = PATH_TO_MERGED_INDEXES,
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_matches: List[Set[int]] = []
        with open(path_to_sets_file, "rt", encoding="utf-8") as merged_matches_file:
            for line in merged_matches_file:
                line = line.strip()
                merged_matches.append(
                    set([int(match_elem) for match_elem in line.split(",")])
                )

        peptide_indexes: List[List[int]] = []
        with open(path_to_indexes_file, "rt", encoding="utf-8") as peptide_indexes_file:
            for line in peptide_indexes_file:
                line = line.strip()
                peptide_indexes.append(
                    [int(match_elem) for match_elem in line.split(",")]
                )

        return merged_matches, peptide_indexes
