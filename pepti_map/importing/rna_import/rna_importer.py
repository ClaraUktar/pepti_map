import logging
from typing import Dict, List, TextIO, Tuple
from Bio.Seq import MutableSeq
import pandas as pd
import gzip

from pepti_map.util.three_frame_translation import get_three_frame_translations


class RNAImporter:
    kmer_length: int
    _cutoff: int
    # TODO: Use objects for values instead? (Better readability)
    _rna_dict: Dict[str, Tuple[List[str], str, int]] = {}
    _rna_dict_translated: Dict[str, Tuple[List[str], str, int, List[int]]] = {}

    def __init__(self, kmer_length: int = 6):
        """
        :param int kmer_length: The k-mer size used during the mapping of peptides to
        RNA. As the RNA is 3-frame translated for the mapping, the k-mer size refers to
        amino acids.
        """
        self.kmer_length = kmer_length

    def reset(self) -> None:
        self._rna_dict = {}
        self._rna_dict_translated = {}
        self._cutoff = -1

    def set_kmer_length(self, kmer_length) -> None:
        self.kmer_length = kmer_length

    def _add_entry_to_dict(self, id: str, sequence: str) -> None:
        duplicate = self._rna_dict.get(sequence)

        if duplicate is not None:
            self._rna_dict[sequence] = (
                duplicate[0] + [id],
                duplicate[1],
                duplicate[2] + 1,
            )
        else:
            self._rna_dict[sequence] = ([id], sequence, 1)

    def _add_entry_to_translated_dict(self, id: str, sequence: str) -> None:
        for translation in get_three_frame_translations(sequence):
            duplicate = self._rna_dict_translated.get(translation[0])

            if duplicate is not None:
                self._rna_dict_translated[translation[0]] = (
                    duplicate[0] + [id],
                    duplicate[1],
                    duplicate[2] + 1,
                    duplicate[3] + [translation[1]],
                )
            else:
                self._rna_dict_translated[translation[0]] = (
                    [id],
                    translation[0],
                    1,
                    [translation[1]],
                )

    def _add_rna_data_to_dict(
        self,
        rna_data: TextIO,
        is_reverse_complement: bool = False,
        should_translate: bool = True,
    ) -> None:
        line_count_for_current_sequence: int = 0
        id = ""
        sequence = ""

        for line in rna_data:
            if line_count_for_current_sequence == 0:
                id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()

                # TODO: Exchange all T for U? (inplace?)

                if self._cutoff > 0:
                    sequence = sequence[0 : self._cutoff]  # noqa: E203

                if is_reverse_complement:
                    sequence = str(
                        MutableSeq(sequence).reverse_complement(inplace=True)
                    )
            # Information from field 2 (line 3) is not needed

            # Always read 4 lines per sequence, as per FASTQ format
            # For now skip quality info (line 4), getting cutoff value supplied by user
            elif line_count_for_current_sequence == 3:
                if should_translate:
                    self._add_entry_to_translated_dict(id, sequence)
                else:
                    self._add_entry_to_dict(id, sequence)

                line_count_for_current_sequence = 0
                id = ""
                sequence = ""
                continue

            line_count_for_current_sequence = line_count_for_current_sequence + 1

    def _fill_dict_from_file(
        self,
        file_path: str,
        is_reverse_complement: bool = False,
        should_translate: bool = True,
    ) -> None:
        with gzip.open(file_path, "rt") as rna_data_gzipped:
            try:
                rna_data_gzipped.read(1)
                rna_data_gzipped.seek(0)
                logging.info(
                    f"Detected gzip file: {file_path}. Reading in compressed format..."
                )
                return self._add_rna_data_to_dict(
                    rna_data_gzipped, is_reverse_complement, should_translate
                )
            except gzip.BadGzipFile:
                logging.info(
                    (
                        f"File {file_path} is not a gzip file. "
                        "Trying to read as uncompressed file..."
                    )
                )

        with open(file_path, "rt") as rna_data:
            return self._add_rna_data_to_dict(
                rna_data, is_reverse_complement, should_translate
            )

    # TODO: Use numpy instead?
    def import_files(
        self, file_paths: List[str], cutoff: int = -1, should_translate: bool = True
    ) -> pd.DataFrame:
        """
        Reads the file(s) given and transforms them into a pandas DataFrame,
        with columns `ids`, `sequence`, and `count`.

        :param List[str] file_paths: A list containing the paths to the files that
        should be imported. In case of single-end sequencing, only one file path
        is expected. In case of paired-end sequencing, two file paths are expected.
        For the reads from the second file, the reverse complement is constructed
        and then the reads are merged with those from the first file into one output.
        :param int cutoff: The position of the last base in the reads after which a
        cutoff should be performed, starting with 1. If given, the value should be > 0.
        :returns A pandas DataFrame. The column `ids` contains all ids with duplicate
        read sequences after cutoff. The column `sequence` contains the corresponding
        sequence after cutoff. The column `count` contains the number of duplicates
        for the sequence.
        :rtype pandas.DataFrame
        :raises ValueError: Raised if the list of file paths does not contain exactly
        1 or 2 entries.
        :raises FileNotFoundError: Raised if no file could be found for a given path.
        """
        self.reset()
        self._cutoff = cutoff

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
            self._fill_dict_from_file(file_path, index == 1, should_translate)

        rna_df: pd.DataFrame
        if should_translate:
            rna_df = pd.DataFrame(
                list(self._rna_dict_translated.values()),
                columns=["ids", "sequence", "count", "frames"],
            ).astype(
                {
                    "ids": "object",
                    "sequence": "string",
                    "count": "uint32",
                    "frames": "object",
                }
            )
        else:
            rna_df = pd.DataFrame(
                list(self._rna_dict.values()), columns=["ids", "sequence", "count"]
            ).astype({"ids": "object", "sequence": "string", "count": "uint32"})
        print(rna_df)
        print(rna_df.info(verbose=True))
        return rna_df
