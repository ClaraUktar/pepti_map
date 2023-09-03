from abc import ABC, abstractmethod
from typing import List, Literal, Set, Union

from datasketch import LeanMinHash

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.distance_matrix_merging import (
    DistanceMatrixMergingMethod,
)
from pepti_map.matching.merging_methods.simple_merging import SimpleMergingMethod


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

    @staticmethod
    def initialize_method(
        method: Literal["agglomerative-clustering", "distance-matrix", "simple"],
        min_hashes: List[Union[LeanMinHash, None]],
        jaccard_index_threshold: float = 0.7,
    ) -> "IMergingMethod":
        if method == "agglomerative-clustering":
            return AgglomerativeClusteringMergingMethod(
                min_hashes, jaccard_index_threshold
            )
        elif method == "distance-matrix":
            return DistanceMatrixMergingMethod(min_hashes, jaccard_index_threshold)
        elif method == "simple":
            return SimpleMergingMethod(min_hashes, jaccard_index_threshold)
        else:
            raise ValueError(
                (
                    "Method must be one of 'agglomerative-clustering', "
                    "'distance-matrix', 'simple'."
                )
            )
