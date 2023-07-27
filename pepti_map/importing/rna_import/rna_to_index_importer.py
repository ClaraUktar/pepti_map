import logging
from typing import List, TextIO
from Bio.Seq import MutableSeq
import gzip

from pepti_map.util.three_frame_translation import get_three_frame_translations


class RNARead:
    id: str
    sequence: str
    is_reverse_complement: bool

    def __init__(self, id: str, sequence: str, is_reverse_complement: bool):
        self.id = id
        self.sequence = sequence
        self.is_reverse_complement = is_reverse_complement


class RNAToIndexImporter:
    kmer_length: int
    _cutoff: int
    _rna_reads: List[RNARead] = []

    def __init__(self, kmer_length: int = 6):
        """
        :param int kmer_length: The k-mer size used during the mapping of peptides to
        RNA. As the RNA is 3-frame translated for the mapping, the k-mer size refers to
        amino acids.
        """
        self.kmer_length = kmer_length

    def reset(self) -> None:
        self._cutoff = -1
        self._rna_reads = []

    def set_kmer_length(self, kmer_length) -> None:
        self.kmer_length = kmer_length

    def _construct_index(self) -> None:
        if len(self._rna_reads) == 0:
            error_message = (
                "No RNA reads available to index. "
                "Please make sure to import a valid, non-empty FASTQ file "
                "before building an index."
            )
            logging.error(error_message)
            raise ValueError(error_message)
        # TODO: Build index

    def _add_rna_data_to_list(
        self,
        rna_data: TextIO,
        is_reverse_complement: bool = False,
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

                # TODO: Construct reverse complement here or in parallelizable part?
                # if is_reverse_complement:
                #     sequence = str(
                #         MutableSeq(sequence).reverse_complement(inplace=True)
                #     )
            # Information from field 2 (line 3) is not needed

            # Always read 4 lines per sequence, as per FASTQ format
            # For now skip quality info (line 4), getting cutoff value supplied by user
            elif line_count_for_current_sequence == 3:
                self._rna_reads.append(RNARead(id, sequence, is_reverse_complement))

                line_count_for_current_sequence = 0
                id = ""
                sequence = ""
                continue

            line_count_for_current_sequence = line_count_for_current_sequence + 1

    def _fill_reads_list_from_file(
        self,
        file_path: str,
        is_reverse_complement: bool = False,
    ) -> None:
        with gzip.open(file_path, "rt") as rna_data_gzipped:
            try:
                rna_data_gzipped.read(1)
                rna_data_gzipped.seek(0)
                logging.info(
                    f"Detected gzip file: {file_path}. Reading in compressed format..."
                )
                return self._add_rna_data_to_list(
                    rna_data_gzipped, is_reverse_complement
                )
            except gzip.BadGzipFile:
                logging.info(
                    (
                        f"File {file_path} is not a gzip file. "
                        "Trying to read as uncompressed file..."
                    )
                )

        with open(file_path, "rt") as rna_data:
            return self._add_rna_data_to_list(rna_data, is_reverse_complement)

    def import_files_to_index(self, file_paths: List[str], cutoff: int = -1) -> None:
        """
        Reads the FASTQ file(s) given and constructs a k-mer index (n-gram index)
        from them.

        :param List[str] file_paths: A list containing the paths to the files that
        should be imported. In case of single-end sequencing, only one file path
        is expected. In case of paired-end sequencing, two file paths are expected.
        For the reads from the second file, the reverse complement is constructed
        and then the reads are merged with those from the first file into one index.
        :param int cutoff: The position of the last base in the reads after which a
        cutoff should be performed, starting with 1. If given, the value should be > 0.
        :raises ValueError: Raised if the list of file paths does not contain exactly
        1 or 2 entries.
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
            self._fill_reads_list_from_file(file_path, index == 1)

        print(self._rna_reads)

        self._construct_index()
