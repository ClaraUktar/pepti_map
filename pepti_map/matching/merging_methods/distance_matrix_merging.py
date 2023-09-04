from typing import List, Set
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash


class DistanceMatrixMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        super(DistanceMatrixMergingMethod, self).__init__(
            min_hashes, jaccard_index_threshold
        )

    def generate_merged_indexes(self) -> List[Set[int]]:
        raise NotImplementedError
