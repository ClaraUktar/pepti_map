import gzip
import logging
import pyfastx
from pathlib import Path
from typing import List, TextIO, Tuple

from Bio.Seq import MutableSeq


class RNAReadsRetriever:
    def __init__(self, filepaths: List[Path], cutoff: int = -1):
        self._filepaths: List[Path] = filepaths
        self._cutoff: int = cutoff

        if len(self._filepaths) > 2 or len(self._filepaths) < 1:
            error_message = (
                "Only one file (single-end sequencing) "
                "or two files (pairend-end sequencing) expected "
                "for the RNA-seq data. "
                f"Received {len(self._filepaths)} files."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        self._first_file_index = pyfastx.Fastq(self._filepaths[0].as_posix())
        self._second_file_index = pyfastx.Fastq(self._filepaths[1].as_posix())

    # @staticmethod
    # def _is_gzip(filepath: Path) -> bool:
    #     return filepath.name.endswith(".gz")

    @staticmethod
    def _sort_paired_end_read_ids(
        read_ids: List[int],
    ) -> Tuple[Tuple[List[int], List[int]], List[int]]:
        first_file_reads = []
        second_file_reads = []
        sorted_original_first_file = []
        sorted_original_second_file = []
        for read_id in read_ids:
            actual_id = int(str(read_id)[0:-1])
            if str(read_id)[-1] == "1":
                first_file_reads.append(actual_id)
                sorted_original_first_file.append(read_id)
            else:
                second_file_reads.append(actual_id)
                sorted_original_second_file.append(read_id)
        first_file_reads.sort()
        second_file_reads.sort()
        sorted_original_first_file.extend(sorted_original_second_file)
        return (first_file_reads, second_file_reads), sorted_original_first_file

    def _process_line(self, line: str, is_reverse_complement: bool) -> str:
        line = line.strip()
        if self._cutoff > 0:
            line = line[0 : self._cutoff]  # noqa: E203
        if is_reverse_complement:
            line = str(MutableSeq(line).reverse_complement(inplace=True))
        return line

    def _append_sequences_for_file(
        self,
        read_ids_for_file: List[int],
        rna_file: TextIO,
        is_reverse_complement: bool,
        sequences: List[str],
    ) -> None:
        current_read_id = read_ids_for_file.pop(0)
        for line_index, line in enumerate(rna_file):
            if (line_index // 4) + 1 != current_read_id:
                continue

            if line_index != (current_read_id - 1) * 4 + 1:
                continue

            sequences.append(self._process_line(line, is_reverse_complement))
            try:
                current_read_id = read_ids_for_file.pop(0)
            except IndexError:
                break

    def get_read_sequences_for_ids(
        self, read_ids: List[int]
    ) -> Tuple[List[int], List[str]]:
        if len(read_ids) == 0:
            return read_ids, []

        sequences = []

        for read_id in read_ids:
            actual_id = int(str(read_id)[0:-1])
            if str(read_id)[-1] == "1":
                # Read id is from first file
                sequence = self._first_file_index[actual_id - 1]
                sequences.append(self._process_line(sequence, False))
            else:
                # Read id is from second file
                sequence = self._second_file_index[actual_id - 1]
                sequences.append(self._process_line(sequence, True))

        return read_ids, sequences

        # (
        #     sorted_read_ids,
        #     original_sorted_ids,
        # ) = RNAReadsRetriever._sort_paired_end_read_ids(read_ids)

        # sequences = []
        # for file_index, filepath in enumerate(self._filepaths):
        #     read_ids_for_file = sorted_read_ids[file_index]
        #     if len(read_ids_for_file) == 0:
        #         continue

        #     if RNAReadsRetriever._is_gzip(filepath):
        #         with gzip.open(filepath, "rt", encoding="utf-8") as rna_file:
        #             self._append_sequences_for_file(
        #                 read_ids_for_file, rna_file, file_index == 1, sequences
        #             )
        #     else:
        #         with open(filepath, "rt", encoding="utf-8") as rna_file:
        #             self._append_sequences_for_file(
        #                 read_ids_for_file, rna_file, file_index == 1, sequences
        #             )

        # return original_sorted_ids, sequences
