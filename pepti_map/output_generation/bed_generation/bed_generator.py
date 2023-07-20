from csv import writer as csv_writer
from typing import List


class BEDEntry:
    chrom: str
    chromStart: int
    chromEnd: int
    name: str

    def __init__(self, chrom: str, chromStart: int, chromEnd: int, name: str):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.name = name

    def columnsToWritableList(self) -> List[str]:
        return [self.chrom, str(self.chromStart), str(self.chromEnd), self.name]


class BEDGenerator:
    entries: List[BEDEntry]

    # TODO: Is the list copied here? -> Optimize
    def __init__(self, entries: List[BEDEntry]):
        self.entries = entries

    def writeEntriesToFile(self, output_filepath: str) -> None:
        with open(output_filepath, "w", newline="") as output_file:
            writer = csv_writer(output_file, delimiter="\t")
            writer.writerow(["# chrom chromStart chromEnd name"])
            for entry in self.entries:
                writer.writerow(entry.columnsToWritableList())


# TODO: DELETE
def main():
    e1 = BEDEntry("chr1", 0, 20, "map1")
    e2 = BEDEntry("chr2", 5, 17, "map2")
    e3 = BEDEntry("chr3", 22, 99, "map3 bla bla")
    gen = BEDGenerator([e1, e2, e3])
    gen.writeEntriesToFile("./bed_test.bed")


if __name__ == "__main__":
    main()
