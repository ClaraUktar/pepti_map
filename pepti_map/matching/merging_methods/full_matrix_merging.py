import logging
from typing import List, Set, Tuple

import numpy as np
import numpy.typing as npt
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod


class FullMatrixMergingMethod(IMergingMethod):
    def __init__(
        self,
        jaccard_calculator: IJaccardIndexCalculator,
        jaccard_index_threshold: float = 0.7,
    ):
        self._merge_indication_matrix: npt.NDArray[np.bool_]
        super(FullMatrixMergingMethod, self).__init__(
            jaccard_calculator, jaccard_index_threshold
        )

    def _process_jaccard_calculator(
        self, jaccard_calculator: IJaccardIndexCalculator
    ) -> None:
        self._merge_indication_matrix = (
            jaccard_calculator.get_jaccard_index_matrix() > self.jaccard_index_threshold
        )

    def _update_set_with_merges(
        self,
        set_to_update: set,
        set_index: int,
        previous_indexes: List[int],
        peptide_indexes: List[int],
        matches: List[Set[int]],
    ) -> Tuple[Set, List[int]]:
        merge_indications = np.where(
            self._merge_indication_matrix[set_index] == True  # noqa: E712
        )[0]
        previous_indexes.append(peptide_indexes[set_index])
        matches[set_index].clear()
        if len(merge_indications) == 0:
            return (set_to_update, previous_indexes)
        for index_to_merge in merge_indications:
            set_to_merge = matches[index_to_merge]
            if len(set_to_merge) == 0:
                continue
            updated_set, previous_indexes = self._update_set_with_merges(
                set_to_merge, index_to_merge, previous_indexes, peptide_indexes, matches
            )
            set_to_update.update(updated_set)
        return (set_to_update, previous_indexes)

    def _generate_result_from_matrix(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_matches: List[Set[int]] = []
        peptide_mappings: List[List[int]] = []
        for current_index, current_set in enumerate(matches):
            merged_set, merged_indexes = self._update_set_with_merges(
                current_set, current_index, [], peptide_indexes, matches
            )
            merged_indexes.sort()
            merged_matches.append(merged_set)
            peptide_mappings.append(merged_indexes)

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
        return self._generate_result_from_matrix(peptide_indexes, matches)
