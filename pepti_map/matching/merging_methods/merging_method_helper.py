from typing import List, Literal
from datasketch import LeanMinHash

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.full_matrix_merging import (
    FullMatrixMergingMethod,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod
from pepti_map.matching.merging_methods.symmetric_matrix_merging import (
    SymmetricMatrixMergingMethod,
)


def get_merging_method(
    method: Literal["agglomerative-clustering", "full-matrix", "symmetric-matrix"],
    min_hashes: List[LeanMinHash],
    jaccard_index_threshold: float = 0.7,
) -> IMergingMethod:
    if method == "agglomerative-clustering":
        return AgglomerativeClusteringMergingMethod(min_hashes, jaccard_index_threshold)
    # TODO: Are both matrix merging methods needed?
    elif method == "full-matrix":
        return FullMatrixMergingMethod(min_hashes, jaccard_index_threshold)
    elif method == "symmetric-matrix":
        return SymmetricMatrixMergingMethod(min_hashes, jaccard_index_threshold)
    else:
        raise ValueError(
            (
                "Method must be one of 'agglomerative-clustering', "
                "'full-matrix', 'symmetric-matrix."
            )
        )
