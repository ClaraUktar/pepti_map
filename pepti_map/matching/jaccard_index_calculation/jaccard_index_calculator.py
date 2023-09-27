from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt


class IJaccardIndexCalculator(ABC):
    JACCARD_INT_MULTIPLICATION_FACTOR = 10000

    @abstractmethod
    def get_jaccard_index(self, first_index: int, second_index: int) -> float:
        """
        Returns the Jaccard Index value for the two sets corresponding
        to the gievn indexes.
        """
        raise NotImplementedError

    @abstractmethod
    def get_jaccard_index_matrix(self) -> npt.NDArray[np.uint16]:
        """
        Returns a matrix of pairwise Jaccard index values for all sets,
        with the indexes corresponding to the set indexes.
        The values are represented as integers, so that for actual value j,
        stored integer value i: j = i / JACCARD_INT_MULTIPLICATION_FACTOR.
        """
        raise NotImplementedError
