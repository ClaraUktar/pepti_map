import logging
import pyfastx
from pathlib import Path
from typing import List, Tuple

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
        if len(self._filepaths) == 2:
            self._second_file_index = pyfastx.Fastq(self._filepaths[1].as_posix())

    def _process_line(self, line: str, is_reverse_complement: bool) -> str:
        line = line.strip()
        if self._cutoff > 0:
            line = line[0 : self._cutoff]  # noqa: E203
        if is_reverse_complement:
            line = str(MutableSeq(line).reverse_complement(inplace=True))
        return line

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
                sequence = self._first_file_index[actual_id - 1].seq
                sequences.append(self._process_line(sequence, False))
            else:
                # Read id is from second file
                sequence = self._second_file_index[actual_id - 1].seq
                sequences.append(self._process_line(sequence, True))

        return read_ids, sequences
