from collections import defaultdict
from typing import List, Tuple


class RNAKmerIndex:
    kmer_index: "defaultdict[str, List[Tuple[str, int, int]]]" = defaultdict(list)

    def clear(self) -> None:
        self.kmer_index.clear()

    def getEntryForKmer(self, kmer: str) -> List[Tuple[str, int, int]]:
        return self.kmer_index[kmer]

    def appendToEntryForKmer(self, kmer: str, entry: Tuple[str, int, int]) -> None:
        self.kmer_index[kmer].append(entry)

    def extendEntryForKmer(self, kmer: str, entry: List[Tuple[str, int, int]]) -> None:
        self.kmer_index[kmer].extend(entry)
