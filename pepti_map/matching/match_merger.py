from typing import List, Literal, Set, Tuple, Union
from datasketch import LeanMinHash, MinHash

from pepti_map.matching.merging_methods.merging_method import IMergingMethod

NUM_BYTES_FOR_MIN_HASH_VALUES = 4


class MatchMerger:
    def __init__(
        self, matches: List[Union[Set[int], None]], jaccard_index_threshold: float = 0.7
    ):
        # TODO: Want to delete matches after merging
        self.matches: List[Union[Set[int], None]] = matches
        self.jaccard_index_threshold: float = jaccard_index_threshold
        self.min_hashes: List[Union[LeanMinHash, None]] = []
        self._create_min_hashes_from_matches()

    def _create_min_hashes_from_matches(self) -> None:
        for match in self.matches:
            if match is None:
                self.min_hashes.append(None)
                continue
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
        # TODO
        pass

    def merge_matches(
        self, method: Literal["agglomerative-clustering", "distance-matrix", "simple"]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        merged_indexes = IMergingMethod.initialize_method(
            method, self.min_hashes, self.jaccard_index_threshold
        ).generate_merged_indexes()
        return self._generate_merged_result_from_indexes(merged_indexes)
