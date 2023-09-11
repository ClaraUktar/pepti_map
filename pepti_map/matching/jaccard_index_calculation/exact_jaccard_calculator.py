import numpy as np
import numpy.typing as npt

from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)


class ExactJaccardCalculator(IJaccardIndexCalculator):
    def __init__(self, exact_jaccard_indexes: npt.NDArray[np.float16]):
        self.exact_jaccard_indexes: npt.NDArray[np.float16] = exact_jaccard_indexes

    def get_jaccard_index(self, first_index: int, second_index: int) -> float:
        return self.exact_jaccard_indexes[first_index][second_index]

    def get_jaccard_index_matrix(self) -> npt.NDArray[np.float16]:
        return self.exact_jaccard_indexes
