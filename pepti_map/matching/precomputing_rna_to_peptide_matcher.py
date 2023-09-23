from pathlib import Path
from typing import List
from itertools import combinations
import numpy as np
import numpy.typing as npt
from pepti_map.constants import PATH_TO_PRECOMPUTED_INTERSECTIONS

from pepti_map.matching.rna_to_peptide_matcher import RNAToPeptideMatcher
from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex


class PrecomputingRNAToPeptideMatcher(RNAToPeptideMatcher):
    def __init__(
        self,
        kmer_index: PeptideKmerIndex,
        number_of_clusters: int,
        peptide_to_cluster_mapping: List[int],
    ):
        super(PrecomputingRNAToPeptideMatcher, self).__init__(
            kmer_index, number_of_clusters, peptide_to_cluster_mapping
        )
        self._precomputed_intersections = np.zeros(
            shape=(number_of_clusters, number_of_clusters), dtype=np.uint32
        )

    def add_peptide_matches_for_rna_read(
        self, rna_read_id: int, rna_read_sequence: str
    ) -> None:
        matched_peptides = set()
        matched_clusters = set()
        for kmer in self._process_rna_read_to_kmers(rna_read_sequence):
            peptide_matches = self._kmer_index.getEntryForKmer(kmer)
            for match in peptide_matches:
                cluster_match = self._peptide_to_cluster_mapping[match]
                if self._matches[cluster_match] is None:
                    self._matches[cluster_match] = set()
                self._matches[
                    cluster_match
                ].add(  # pyright: ignore[reportOptionalMemberAccess]
                    rna_read_id
                )
                matched_peptides.add(match)
                matched_clusters.add(cluster_match)
        for match in matched_peptides:
            self._matches_per_peptide[match] += 1
        for intersecting_set_pair in list(combinations(matched_clusters, 2)):
            # TODO: Can this be improved to not set both entries?
            self._precomputed_intersections[intersecting_set_pair[0]][
                intersecting_set_pair[1]
            ] += 1
            self._precomputed_intersections[intersecting_set_pair[1]][
                intersecting_set_pair[0]
            ] += 1

    def get_precomputed_intersections(self) -> npt.NDArray[np.uint32]:
        for cluster_index in range(0, self._precomputed_intersections.shape[0]):
            if self._matches[cluster_index] is None:
                continue
            self._precomputed_intersections[cluster_index][cluster_index] = len(
                self._matches[cluster_index]  # pyright: ignore[reportGeneralTypeIssues]
            )
        return self._precomputed_intersections

    def save_precomputed_intersections(
        self, path_to_intersections: Path = PATH_TO_PRECOMPUTED_INTERSECTIONS
    ) -> None:
        np.savez_compressed(
            path_to_intersections,
            intersections=self.get_precomputed_intersections(),
        )

    @staticmethod
    def load_precomputed_intersections(
        path_to_intersections: Path = PATH_TO_PRECOMPUTED_INTERSECTIONS,
    ) -> npt.NDArray[np.uint32]:
        return np.load(path_to_intersections)["intersections"]
