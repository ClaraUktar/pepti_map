import gzip
import logging
from pathlib import Path
from typing import List, Tuple, Union

from Bio.Seq import MutableSeq


class LazyRNAReader(object):
    def __init__(self, filepaths: List[Path], cutoff: int = -1):
        if len(filepaths) > 2 or len(filepaths) < 1:
            error_message = (
                "Only one file (single-read sequencing) "
                "or two files (pairend-end sequencing) expected "
                "for the RNA-seq data. "
                f"Received {len(filepaths)} files."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        self._filepaths = filepaths

        self._open_filehandle = None

        self._cutoff: int = cutoff

        self._line_count_for_current_sequence: int = 0
        self._sequence_id: str = ""
        self._sequence: str = ""

    @staticmethod
    def _is_gzip(filepath: Path) -> bool:
        with gzip.open(filepath, "rt", encoding="utf-8") as rna_data_gzipped:
            try:
                rna_data_gzipped.read(1)
                is_gzip = True
            except gzip.BadGzipFile:
                is_gzip = False

        return is_gzip

    def _process_line(
        self, line: str, is_reverse_complement: bool
    ) -> Union[Tuple[str, str], None]:
        if self._line_count_for_current_sequence == 0:
            self._sequence_id = line.strip()
        elif self._line_count_for_current_sequence == 1:
            self._sequence = line.strip()

            # TODO: Exchange all T for U? (inplace?)

            if self._cutoff > 0:
                self._sequence = self._sequence[0 : self._cutoff]  # noqa: E203

            if is_reverse_complement:
                self._sequence = str(
                    MutableSeq(self._sequence).reverse_complement(inplace=True)
                )
        # Information from field 2 (line 3) is not needed

        # Always read 4 lines per sequence, as per FASTQ format
        # For now skip quality info (line 4), getting cutoff value supplied by user
        elif self._line_count_for_current_sequence == 3:
            # yield (sequence_id, sequence)
            self._line_count_for_current_sequence = 0
            # sequence_id = ""
            # sequence = ""
            return (self._sequence_id, self._sequence)

        self._line_count_for_current_sequence = (
            self._line_count_for_current_sequence + 1
        )

    def __iter__(self):
        for index, filepath in enumerate(self._filepaths):
            # if self._is_gzip(filepath):
            #     self._open_filehandle = gzip.open(filepath, "rt", encoding="utf-8")
            # else:
            self._open_filehandle = open(filepath, "rt", encoding="utf-8")

            for line in self._open_filehandle:
                processed_line = self._process_line(line, index == 1)
                if processed_line is not None:
                    yield processed_line

            self._open_filehandle.close()

    def __del__(self):
        if self._open_filehandle is not None:
            self._open_filehandle.close()
