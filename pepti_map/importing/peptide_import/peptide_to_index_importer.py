from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.util.k_mer import split_into_kmer


class PeptideToIndexImporter:
    def __init__(self, kmer_length: int = 7):
        self.kmer_length: int = kmer_length
        self.kmer_index: PeptideKmerIndex = PeptideKmerIndex(kmer_length)

    def reset(self) -> None:
        self.kmer_index.clear()

    def import_file_to_index(self, file_path: str) -> PeptideKmerIndex:
        # TODO: Should reset not happen automatically
        # to enable import of multiple files
        self.reset()
        number_of_peptides = 0

        with open(file_path, "rt", encoding="utf-8") as peptide_file:
            for index, line in enumerate(peptide_file):
                # TODO: Exchange all I for L?
                sequence = line.strip()
                for kmer in split_into_kmer(sequence, self.kmer_length):
                    self.kmer_index.appendToEntryForKmer(kmer[0], index)
                if sequence != "":
                    number_of_peptides += 1

        # TODO: Do this in a prettier way
        self.kmer_index.number_of_peptides = number_of_peptides

        return self.kmer_index
