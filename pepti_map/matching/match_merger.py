from typing import List, Literal, Set, Tuple, Union
from datasketch import LeanMinHash, MinHash

from pepti_map.matching.merging_methods.merging_method_helper import get_merging_method

NUM_BYTES_FOR_MIN_HASH_VALUES = 4


class MatchMerger:
    def __init__(
        self, matches: List[Union[Set[int], None]], jaccard_index_threshold: float = 0.7
    ):
        # TODO: Want to delete matches after merging
        self.jaccard_index_threshold: float = jaccard_index_threshold
        self.peptide_indexes: List[int] = []
        self.matches: List[Set[int]] = []
        self._postprocess_matches(matches)

        self.min_hashes: List[LeanMinHash] = []
        self._create_min_hashes_from_matches()

    def _postprocess_matches(
        self, uncleaned_matches: List[Union[Set[int], None]]
    ) -> None:
        for peptide_index, uncleaned_match in enumerate(uncleaned_matches):
            if uncleaned_match is None:
                continue
            self.matches.append(uncleaned_match)
            self.peptide_indexes.append(peptide_index)

    def _create_min_hashes_from_matches(self) -> None:
        for match in self.matches:
            min_hash = MinHash()
            min_hash.update_batch(
                [
                    element_to_add.to_bytes(NUM_BYTES_FOR_MIN_HASH_VALUES, "big")
                    for element_to_add in match
                ]
            )
            self.min_hashes.append(LeanMinHash(min_hash))

    def _generate_merged_result_from_indexes(
        self, merged_indexes: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_matches = []
        peptide_mappings = []
        for merged_index_set in merged_indexes:
            current_set = set()
            current_peptides = []

            for peptide_id in merged_index_set:
                true_peptide_id = self.peptide_indexes[peptide_id]
                set_to_merge = self.matches[true_peptide_id]
                current_set.update(set_to_merge)
                current_peptides.append(true_peptide_id)

            merged_matches.append(current_set)
            peptide_mappings.append(current_peptides)

        return (merged_matches, peptide_mappings)

    def merge_matches(
        self, method: Literal["agglomerative-clustering", "distance-matrix", "simple"]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_indexes = get_merging_method(
            method, self.min_hashes, self.jaccard_index_threshold
        ).generate_merged_indexes()
        return self._generate_merged_result_from_indexes(merged_indexes)
