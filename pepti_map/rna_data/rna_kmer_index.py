from collections import defaultdict
import gzip
from typing import List, Tuple


class RNAKmerIndex:
    kmer_index: "defaultdict[str, List[Tuple[str, int, int]]]"

    def __init__(self):
        self.kmer_index = defaultdict(list)

    def clear(self) -> None:
        self.kmer_index.clear()

    def getEntryForKmer(self, kmer: str) -> List[Tuple[str, int, int]]:
        return self.kmer_index[kmer]

    def appendToEntryForKmer(self, kmer: str, entry: Tuple[str, int, int]) -> None:
        self.kmer_index[kmer].append(entry)

    def extendEntryForKmer(self, kmer: str, entry: List[Tuple[str, int, int]]) -> None:
        self.kmer_index[kmer].extend(entry)

    def dump_index_to_file(self, file_path: str) -> None:
        with gzip.open(file_path, "wt", encoding="utf-8") as index_file:
            for item in self.kmer_index.items():
                file_entry = f"{item[0]}\t"
                for index, value in enumerate(item[1]):
                    file_entry += ",".join([value[0], str(value[1]), str(value[2])])
                    if index != len(item[1]) - 1:
                        file_entry += ";"
                index_file.write(file_entry)
                index_file.write("\n")

    def load_index_from_file(self, file_path: str) -> "RNAKmerIndex":
        self.clear()
        with gzip.open(file_path, "rt", encoding="utf-8") as index_file:
            for line in index_file:
                kmer, values = line.split(sep="\t")
                for value in values.split(sep=";"):
                    sequence_id, frame, position = value.split(sep=",")
                    self.appendToEntryForKmer(
                        kmer, (sequence_id, int(frame), int(position))
                    )
        return self
