from typing import List
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
        # TOOD: Init precomputed intersections array

    def add_peptide_matches_for_rna_read(
        self, rna_read_id: int, rna_read_sequence: str
    ) -> None:
        matched_peptides = set()
        matched_clusters = set()
        for kmer in self._process_rna_read_to_kmers(rna_read_sequence):
            peptide_matches = self._kmer_index.getEntryForKmer(kmer)
            for match in peptide_matches:
                cluster_match = self._peptide_to_cluster_mapping[match]
                if self.matches[cluster_match] is None:
                    self.matches[cluster_match] = set()
                self.matches[
                    cluster_match
                ].add(  # pyright: ignore[reportOptionalMemberAccess]
                    rna_read_id
                )
                matched_peptides.add(match)
                matched_clusters.add(cluster_match)
        for match in matched_peptides:
            self._matches_per_peptide[match] += 1
        # TODO: Process matched_clusters for precomputed intersections
        # TODO: Leave out symmetrical entries (e.g. (1, 1)
        # and process later by adding set size?)
