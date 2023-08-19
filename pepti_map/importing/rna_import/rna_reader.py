import gzip
import logging
from Bio.Seq import MutableSeq
from typing import Generator, List, TextIO, Tuple


# TODO: Write docstrings and test
class RNAReader:
    def __init__(self):
        self.cutoff: int = -1

    def _get_ids_and_sequences(
        self, rna_data: TextIO, is_reverse_complement: bool = False
    ) -> Generator[Tuple[str, str], None, None]:
        line_count_for_current_sequence: int = 0
        sequence_id = ""
        sequence = ""

        for line in rna_data:
            if line_count_for_current_sequence == 0:
                sequence_id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()

                # TODO: Exchange all T for U? (inplace?)

                if self.cutoff > 0:
                    sequence = sequence[0 : self.cutoff]  # noqa: E203

                if is_reverse_complement:
                    sequence = str(
                        MutableSeq(sequence).reverse_complement(inplace=True)
                    )
            # Information from field 2 (line 3) is not needed

            # Always read 4 lines per sequence, as per FASTQ format
            # For now skip quality info (line 4), getting cutoff value supplied by user
            elif line_count_for_current_sequence == 3:
                yield (sequence_id, sequence)

                line_count_for_current_sequence = 0
                sequence_id = ""
                sequence = ""
                continue

            line_count_for_current_sequence = line_count_for_current_sequence + 1

    def _is_gzip(self, file_path: str) -> bool:
        with gzip.open(file_path, "rt", encoding="utf-8") as rna_data_gzipped:
            try:
                rna_data_gzipped.read(1)
                is_gzip = True
            except gzip.BadGzipFile:
                is_gzip = False

        return is_gzip

    def _read_file_in_correct_format(
        self, file_path: str, is_reverse_complement: bool = False
    ) -> Generator[Tuple[str, str], None, None]:
        if self._is_gzip(file_path):
            logging.info(
                f"Detected gzip file: {file_path}. Reading in compressed format..."
            )
            with gzip.open(file_path, "rt", encoding="utf-8") as rna_data_gzipped:
                for sequence_entry in self._get_ids_and_sequences(
                    rna_data_gzipped, is_reverse_complement
                ):
                    yield sequence_entry
        else:
            logging.info(
                (
                    f"File {file_path} is not a gzip file. "
                    "Trying to read as uncompressed file..."
                )
            )
            with open(file_path, "rt", encoding="utf-8") as rna_data:
                for sequence_entry in self._get_ids_and_sequences(
                    rna_data, is_reverse_complement
                ):
                    yield sequence_entry

    def read_lines(self, file_paths: List[str], cutoff: int = -1):
        self.cutoff = cutoff

        if len(file_paths) > 2 or len(file_paths) < 1:
            error_message = (
                "Only one file (single-read sequencing) "
                "or two files (pairend-end sequencing) expected "
                "for the RNA-seq data. "
                f"Received {len(file_paths)} files."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        for index, file_path in enumerate(file_paths):
            for sequence_entry in self._read_file_in_correct_format(
                file_path, index == 1
            ):
                yield sequence_entry
