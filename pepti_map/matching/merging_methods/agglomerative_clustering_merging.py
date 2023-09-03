from typing import List, Set, Union
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash


class AgglomerativeClusteringMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[Union[LeanMinHash, None]],
        jaccard_index_threshold: float = 0.7,
    ):
        super(AgglomerativeClusteringMergingMethod, self).__init__(
            min_hashes, jaccard_index_threshold
        )

    def generate_merged_indexes(self) -> List[Set[int]]:
        raise NotImplementedError