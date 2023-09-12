import logging
from typing import Literal
from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)

from pepti_map.matching.merging_methods.agglomerative_clustering_merging import (
    AgglomerativeClusteringMergingMethod,
)
from pepti_map.matching.merging_methods.full_matrix_merging import (
    FullMatrixMergingMethod,
)
from pepti_map.matching.merging_methods.merging_method import IMergingMethod


def get_merging_method(
    method: Literal["agglomerative-clustering", "full-matrix"],
    jaccard_calculator: IJaccardIndexCalculator,
    jaccard_index_threshold: float = 0.7,
) -> IMergingMethod:
    if method == "agglomerative-clustering":
        return AgglomerativeClusteringMergingMethod(
            jaccard_calculator, jaccard_index_threshold
        )
    elif method == "full-matrix":
        return FullMatrixMergingMethod(jaccard_calculator, jaccard_index_threshold)
    else:
        error_message = (
            "Method must be one of 'agglomerative-clustering', 'full-matrix'."
        )
        logging.error(error_message)
        raise ValueError(error_message)
