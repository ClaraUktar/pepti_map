from typing import List, Literal
from datasketch import LeanMinHash

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.matrix_merging import MatrixMergingMethod
from pepti_map.matching.merging_methods.merging_method import IMergingMethod


def get_merging_method(
    method: Literal["agglomerative-clustering", "matrix"],
    min_hashes: List[LeanMinHash],
    jaccard_index_threshold: float = 0.7,
) -> IMergingMethod:
    if method == "agglomerative-clustering":
        return AgglomerativeClusteringMergingMethod(min_hashes, jaccard_index_threshold)
    elif method == "matrix":
        return MatrixMergingMethod(min_hashes, jaccard_index_threshold)
    else:
        raise ValueError(
            ("Method must be one of 'agglomerative-clustering', " "'distance-matrix'.")
        )
