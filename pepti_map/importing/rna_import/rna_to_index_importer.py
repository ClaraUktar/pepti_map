import itertools
import logging
from typing import List, TextIO, Tuple
from Bio.Seq import MutableSeq
import gzip
from pepti_map.rna_data.rna_kmer_index import RNAKmerIndex
from pepti_map.util.k_mer import split_into_kmer

from pepti_map.util.three_frame_translation import get_three_frame_translations


class RNAToIndexImporter:
    kmer_length: int
    _cutoff: int

    def __init__(self, kmer_length: int = 7):
        """
        :param int kmer_length: The k-mer size used during the mapping of peptides to
        RNA. As the RNA is 3-frame translated for the mapping, the k-mer size refers to
        amino acids.
        """
        self.kmer_length = kmer_length
        self.kmer_index = RNAKmerIndex(kmer_length)

    def reset(self) -> None:
        self._cutoff = -1
        self.kmer_index.clear()

    def set_kmer_length(self, kmer_length: int) -> None:
        self.kmer_length = kmer_length

    def _generate_kmers(
        self, sequence: str, frame: int, sequence_id: str
    ) -> List[Tuple[str, str, int, int]]:
        kmers = split_into_kmer(sequence, self.kmer_length)
        return [(kmer[0], sequence_id, frame, kmer[1]) for kmer in kmers]

    def _process_rna_read_to_kmer(
        self, sequence_id: str, sequence: str, is_reverse_complement: bool
    ) -> List[Tuple[str, str, int, int]]:
        """
        TODO
        """
        # TODO: Exchange all T for U? (inplace?)
        # TODO: Construct reverse complement here or during file reading?
        if is_reverse_complement:
            sequence = str(MutableSeq(sequence).reverse_complement(inplace=True))

        translations = get_three_frame_translations(sequence)
        return list(
            itertools.chain.from_iterable(
                [
                    self._generate_kmers(translation[0], translation[1], sequence_id)
                    for translation in translations
                ]
            )
        )

    def _add_rna_read_to_index(
        self, sequence_id: str, sequence: str, is_reverse_complement: bool
    ) -> None:
        kmers = self._process_rna_read_to_kmer(
            sequence_id, sequence, is_reverse_complement
        )
        for kmer in kmers:
            self.kmer_index.appendToEntryForKmer(kmer[0], (kmer[1], kmer[2], kmer[3]))

    def _add_rna_data_to_index(
        self,
        rna_data: TextIO,
        is_reverse_complement: bool = False,
    ) -> None:
        line_count_for_current_sequence: int = 0
        sequence_id = ""
        sequence = ""

        for line in rna_data:
            if line_count_for_current_sequence == 0:
                sequence_id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()

                if self._cutoff > 0:
                    sequence = sequence[0 : self._cutoff]  # noqa: E203

            # Information from field 2 (line 3) is not needed

            # Always read 4 lines per sequence, as per FASTQ format
            # For now skip quality info (line 4), getting cutoff value supplied by user
            elif line_count_for_current_sequence == 3:
                self._add_rna_read_to_index(
                    sequence_id, sequence, is_reverse_complement
                )

                line_count_for_current_sequence = 0
                sequence_id = ""
                sequence = ""
                continue

            line_count_for_current_sequence = line_count_for_current_sequence + 1

    def _fill_index_from_file(
        self,
        file_path: str,
        is_reverse_complement: bool = False,
    ) -> None:
        with gzip.open(file_path, "rt", encoding="utf-8") as rna_data_gzipped:
            try:
                rna_data_gzipped.read(1)
                rna_data_gzipped.seek(0)
                logging.info(
                    f"Detected gzip file: {file_path}. Reading in compressed format..."
                )
                return self._add_rna_data_to_index(
                    rna_data_gzipped, is_reverse_complement
                )
            except gzip.BadGzipFile:
                logging.info(
                    (
                        f"File {file_path} is not a gzip file. "
                        "Trying to read as uncompressed file..."
                    )
                )

        with open(file_path, "rt", encoding="utf-8") as rna_data:
            return self._add_rna_data_to_index(rna_data, is_reverse_complement)

    def import_files_to_index(
        self, file_paths: List[str], cutoff: int = -1
    ) -> RNAKmerIndex:
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
        :returns An `RNAKmerIndex` object containing the constructed index.
        :rtype RNAKmerIndex
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
            self._fill_index_from_file(file_path, index == 1)

        return self.kmer_index
