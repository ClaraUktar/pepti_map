from collections import defaultdict
import gzip
from typing import List


class PeptideKmerIndex:
    def __init__(self, kmer_length: int = 7):
        self.kmer_index: "defaultdict[str, List[int]]" = defaultdict(list)
        self.kmer_length: int = kmer_length
        # TODO: Do this in a prettier way
        self.number_of_peptides: int = -1

    def clear(self) -> None:
        self.kmer_index.clear()

    def getEntryForKmer(self, kmer: str) -> List[int]:
        return self.kmer_index[kmer]

    def appendToEntryForKmer(
        self, kmer: str, entry: int, check_for_duplicates=False
    ) -> None:
        # TODO: Is there a better option for dedplication? Use sets instead?
        if check_for_duplicates:
            if entry in self.kmer_index[kmer]:
                return
        self.kmer_index[kmer].append(entry)

    # TODO: This does not check for duplicates at the moment! Should it?
    def extendEntryForKmer(self, kmer: str, entry: List[int]) -> None:
        self.kmer_index[kmer].extend(entry)

    def dump_index_to_file(self, filepath: str) -> None:
        with gzip.open(filepath, "wt", encoding="utf-8") as index_file:
            index_file.write(f"# n_peptides={self.number_of_peptides}\n")
            for item in self.kmer_index.items():
                file_entry = f"{item[0]}\t"
                for index, value in enumerate(item[1]):
                    file_entry += str(value)
                    if index != len(item[1]) - 1:
                        file_entry += ";"
                index_file.write(file_entry)
                index_file.write("\n")

    @classmethod
    def load_index_from_file(cls, filepath: str) -> "PeptideKmerIndex":
        kmer_index: cls
        with gzip.open(filepath, "rt", encoding="utf-8") as index_file:
            number_of_peptides = int(index_file.readline().replace("# n_peptides=", ""))
            kmer_length = len(index_file.readline().split(sep="\t")[0])
            kmer_index = cls(kmer_length)
            kmer_index.number_of_peptides = number_of_peptides
            index_file.seek(0)
            next(index_file)

            for line in index_file:
                kmer, values = line.split(sep="\t")
                kmer_index.extendEntryForKmer(
                    kmer, [int(value) for value in values.split(";")]
                )
        return kmer_index
