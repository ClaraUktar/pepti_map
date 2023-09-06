from typing import List, Set, Tuple
from pepti_map.matching.merging_methods.merging_method import IMergingMethod

from datasketch import LeanMinHash


class SimpleMergingMethod(IMergingMethod):
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        super(SimpleMergingMethod, self).__init__(min_hashes, jaccard_index_threshold)

    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        raise NotImplementedError
