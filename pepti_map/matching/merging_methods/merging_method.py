from abc import ABC, abstractmethod
from typing import List, Set, Tuple

from pepti_map.matching.jaccard_index_calculation.jaccard_index_calculator import (
    IJaccardIndexCalculator,
)


class IMergingMethod(ABC):
    def __init__(
        self,
        jaccard_calculator: IJaccardIndexCalculator,
        jaccard_index_threshold: float = 0.7,
    ):
        self.jaccard_index_threshold = jaccard_index_threshold
        self._process_jaccard_calculator(jaccard_calculator)

    @abstractmethod
    def _process_jaccard_calculator(
        self, jaccard_calculator: IJaccardIndexCalculator
    ) -> None:
        """
        Process the given Jaccard index calculator as needed for the merging method.
        """
        raise NotImplementedError

    @abstractmethod
    def generate_merged_result(
        self, peptide_indexes: List[int], matches: List[Set[int]]
    ) -> Tuple[List[Set[int]], List[List[int]]]:
        """
        Returns a list of merged sets plus a list of corresponding
        peptide indexes per merged set.
        """
        raise NotImplementedError
