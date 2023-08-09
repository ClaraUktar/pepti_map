import logging
from typing import List, Union
import pandas as pd

from pepti_map.rna_data.rna_kmer_index import RNAKmerIndex
from pepti_map.util.k_mer import split_into_kmer


class PeptideToRNAMatcher:
    peptides: pd.DataFrame
    kmer_index: Union[RNAKmerIndex, None]
    kmer_length: int

    def __init__(self, peptides: pd.DataFrame):
        self.peptides = peptides
        self.kmer_index = None
        self.kmer_length = 0

    def _find_rna_read_matches_for_peptide(self, peptide: str) -> List[str]:
        try:
            assert self.kmer_index is not None and self.kmer_length != 0
        except AssertionError as assertion_error:
            logging.error(
                "A k-mer index needs to be provided to enable the matching. "
                "Expected a RNAKmerIndex object, but received None"
            )
            raise assertion_error

        matches = []
        for kmer in split_into_kmer(peptide, self.kmer_length):
            # TODO: Add the sequences from the index, we probably don't need
            # the tuples there anymore, only list of ids as value
            pass
        return matches

    def find_rna_read_matches_for_peptides(
        self, kmer_index: RNAKmerIndex, kmer_length: int
    ) -> pd.DataFrame:
        self.kmer_index = kmer_index
        # TODO: Extract k-mer length from index instead?
        # (Important for already existing index)
        # Might make sense to use a standard dict after all
        self.kmer_length = kmer_length

        # TODO: Apply for each row in df

        self.kmer_index = None
        return self.peptides
