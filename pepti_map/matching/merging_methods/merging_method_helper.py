import logging
from typing import List, Literal
from datasketch import LeanMinHash

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.full_matrix_merging import (
    FullMatrixMergingMethod,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod


def get_merging_method(
    method: Literal["agglomerative-clustering", "full-matrix"],
    min_hashes: List[LeanMinHash],
    jaccard_index_threshold: float = 0.7,
) -> IMergingMethod:
    if method == "agglomerative-clustering":
        return AgglomerativeClusteringMergingMethod(min_hashes, jaccard_index_threshold)
    elif method == "full-matrix":
        return FullMatrixMergingMethod(min_hashes, jaccard_index_threshold)
    else:
        error_message = (
            "Method must be one of 'agglomerative-clustering', 'full-matrix'."
        )
        logging.error(error_message)
        raise ValueError(error_message)
