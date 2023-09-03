from typing import List, Literal, Union
from datasketch import LeanMinHash

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.distance_matrix_merging import (
    DistanceMatrixMergingMethod,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod
from pepti_map.matching.merging_methods.simple_merging import SimpleMergingMethod


def get_merging_method(
    method: Literal["agglomerative-clustering", "distance-matrix", "simple"],
    min_hashes: List[Union[LeanMinHash, None]],
    jaccard_index_threshold: float = 0.7,
) -> IMergingMethod:
    if method == "agglomerative-clustering":
        return AgglomerativeClusteringMergingMethod(min_hashes, jaccard_index_threshold)
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
