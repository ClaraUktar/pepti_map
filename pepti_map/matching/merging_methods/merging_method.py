from abc import ABC, abstractmethod
from typing import List, Set, Union

from datasketch import LeanMinHash


class IMergingMethod(ABC):
    @abstractmethod
    def __init__(
        self,
        min_hashes: List[Union[LeanMinHash, None]],
        jaccard_index_threshold: float = 0.7,
    ):
        self.min_hashes: List[Union[LeanMinHash, None]] = min_hashes
        self.jaccard_index_threshold = jaccard_index_threshold

    @abstractmethod
    def generate_merged_indexes(self) -> List[Set[int]]:
        """
        Returns a list of merged sets containing the indexes corresponding to
        the original sets/MinHashes
        (= peptide indexes)
        """
        raise NotImplementedError
