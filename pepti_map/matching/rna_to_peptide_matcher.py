from typing import Generator, List, Set

from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.util.k_mer import split_into_kmer
from pepti_map.util.three_frame_translation import get_three_frame_translations


class RNAToPeptideMatcher:
    def __init__(self, kmer_index: PeptideKmerIndex, number_of_peptides: int):
        self.kmer_index: PeptideKmerIndex = kmer_index
        self.matches: List[Set[str]] = [set() for _ in range(0, number_of_peptides)]

    def _process_rna_read_to_kmers(self, sequence: str) -> Generator[str, None, None]:
        # TODO: Exchange all T for U? (inplace?)
        # TODO: Construct reverse complement here or during file reading?
        for translation in get_three_frame_translations(sequence):
            for kmer in split_into_kmer(translation[0], self.kmer_index.kmer_length):
                yield kmer[0]

    def add_peptide_matches_for_rna_read(
        self, rna_read_id: str, rna_read_sequence: str
    ) -> None:
        for kmer in self._process_rna_read_to_kmers(rna_read_sequence):
            peptide_matches = self.kmer_index.getEntryForKmer(kmer)
            for match in peptide_matches:
                self.matches[match].add(rna_read_id)
