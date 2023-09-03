from typing import List, Set, Union
from datasketch import LeanMinHash, MinHash

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
