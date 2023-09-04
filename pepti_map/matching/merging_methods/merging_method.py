from abc import ABC, abstractmethod
from typing import List, Set, Tuple

from datasketch import LeanMinHash


class IMergingMethod(ABC):
    @abstractmethod
    def __init__(
        self,
        min_hashes: List[LeanMinHash],
        jaccard_index_threshold: float = 0.7,
    ):
        self.min_hashes: List[LeanMinHash] = min_hashes
        self.jaccard_index_threshold = jaccard_index_threshold

    @abstractmethod
    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        """
        Returns a list of merged sets plus a list of corresponding
        peptide indexes per merged set.
        """
        raise NotImplementedError
