import gzip
from pathlib import Path
from unittest.mock import mock_open, patch
from io import StringIO
import pandas as pd

import pytest

from pepti_map.importing.rna_import.rna_importer import RNAImporter
from pepti_map.importing.rna_import.lazy_rna_reader import LazyRNAReader
from pepti_map.importing.rna_import.rna_reads_retriever import RNAReadsRetriever
from pepti_map.importing.rna_import.rna_to_index_importer import RNAToIndexImporter
from pepti_map.importing.rna_import.testdata_rna_importer import (
    EXPECTED_RESULT_DF_PAIRED_END,
    EXPECTED_RESULT_DF_SINGLE_END,
    EXPECTED_RESULT_DF_SINGLE_END_CUTOFF,
    EXPECTED_RESULT_INDEX_PAIRED_END,
    EXPECTED_RESULT_INDEX_SINGLE_END,
    EXPECTED_RESULT_INDEX_SINGLE_END_CUTOFF,
    EXPECTED_RESULT_INDEX_SINGLE_END_SHORT_KMERS,
    EXPECTED_RESULT_LINES_PAIRED_END,
    EXPECTED_RESULT_LINES_SINGLE_END,
    EXPECTED_RESULT_LINES_SINGLE_END_CUTOFF,
    MOCK_FILE_1_CONTENT,
    MOCK_FILE_2_CONTENT,
)


class TestRNAImporter:
    @pytest.fixture(autouse=True)
    def _init_rna_importer(self):
        self.rna_importer = RNAImporter()

    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files([])

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files(["file1", "file2", "file3"])

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_import_single_end_file(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file"], should_translate=False
        )
        pd.testing.assert_frame_equal(result_df, EXPECTED_RESULT_DF_SINGLE_END)

    @patch(
        "gzip.open",
        side_effect=[
            StringIO(MOCK_FILE_1_CONTENT),
            StringIO(MOCK_FILE_2_CONTENT),
        ],
    )
    def test_import_paired_end_file(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file1", "path/to/file2"], should_translate=False
        )
        pd.testing.assert_frame_equal(result_df, EXPECTED_RESULT_DF_PAIRED_END)

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_cutoff(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file"], cutoff=10, should_translate=False
        )
        pd.testing.assert_frame_equal(result_df, EXPECTED_RESULT_DF_SINGLE_END_CUTOFF)


class TestRNAToIndexImporter:
    @pytest.fixture(autouse=True)
    def _init_rna_importer(self):
        self.rna_importer = RNAToIndexImporter(kmer_length=7)

    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files_to_index([])

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files_to_index(["file1", "file2", "file3"])

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_import_single_end_file(self, _):
        resulting_index = self.rna_importer.import_files_to_index(["path/to/file"])
        assert resulting_index.kmer_index == EXPECTED_RESULT_INDEX_SINGLE_END

    @patch(
        "gzip.open",
        side_effect=[
            StringIO(MOCK_FILE_1_CONTENT),
            StringIO(MOCK_FILE_2_CONTENT),
        ],
    )
    def test_import_paired_end_file(self, _):
        resulting_index = self.rna_importer.import_files_to_index(
            ["path/to/file1", "path/to/file2"]
        )
        assert resulting_index.kmer_index == EXPECTED_RESULT_INDEX_PAIRED_END

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_cutoff(self, _):
        resulting_index = self.rna_importer.import_files_to_index(
            ["path/to/file"], cutoff=30
        )
        assert resulting_index.kmer_index == EXPECTED_RESULT_INDEX_SINGLE_END_CUTOFF

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_kmer_len(self, _):
        resulting_index = RNAToIndexImporter(kmer_length=3).import_files_to_index(
            ["path/to/file"]
        )
        assert (
            resulting_index.kmer_index == EXPECTED_RESULT_INDEX_SINGLE_END_SHORT_KMERS
        )


class TestLazyRNAReader:
    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            for _ in LazyRNAReader([]):
                continue

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            for _ in LazyRNAReader([Path("file1"), Path("file2"), Path("file3")]):
                continue

    @patch("gzip.open", return_value=StringIO(MOCK_FILE_1_CONTENT))
    def test_read_single_end_file_gzipped(self, _):
        read_lines = []
        for line in LazyRNAReader([Path("path/to/file.gz")]):
            read_lines.append(line)
        assert read_lines == EXPECTED_RESULT_LINES_SINGLE_END

    @patch(
        "builtins.open",
        side_effect=[
            StringIO(MOCK_FILE_1_CONTENT),
            StringIO(MOCK_FILE_2_CONTENT),
        ],
    )
    def test_read_paired_end_file(self, _):
        read_lines = []
        for line in LazyRNAReader([Path("path/to/file1"), Path("path/to/file2")]):
            read_lines.append(line)
        assert read_lines == EXPECTED_RESULT_LINES_PAIRED_END

    @patch("builtins.open", mock_open(read_data=MOCK_FILE_1_CONTENT))
    def test_read_single_end_file_cutoff(self):
        read_lines = []
        for line in LazyRNAReader([Path("path/to/file")], cutoff=10):
            read_lines.append(line)
        assert read_lines == EXPECTED_RESULT_LINES_SINGLE_END_CUTOFF


class TestRNAReadsRetriever:
    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            RNAReadsRetriever([])

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            RNAReadsRetriever([Path("file1"), Path("file2"), Path("file3")])

    def test_single_end_file_gzipped(self, tmp_path):
        with gzip.open(tmp_path / "test_file.gz", "wt", encoding="utf-8") as test_file:
            test_file.write(
                MOCK_FILE_1_CONTENT  # pyright: ignore[reportGeneralTypeIssues]
            )

        reads_retriever = RNAReadsRetriever([tmp_path / "test_file.gz"])
        expected_ids = [11, 41, 51]
        expected_reads = [
            (
                "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTT"
                "GGTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
            ),
            (
                "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGT"
                "ACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
            ),
            (
                "ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGT"
                "GTACAAACAGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
            ),
        ]
        assert reads_retriever.get_read_sequences_for_ids([11, 41, 51]) == (
            expected_ids,
            expected_reads,
        )

    def test_paired_end_file(self, tmp_path):
        with gzip.open(
            tmp_path / "test_file_1.gz", "wt", encoding="utf-8"
        ) as test_file_1:
            test_file_1.write(
                MOCK_FILE_1_CONTENT  # pyright: ignore[reportGeneralTypeIssues]
            )
        with gzip.open(
            tmp_path / "test_file_2.gz", "wt", encoding="utf-8"
        ) as test_file_2:
            test_file_2.write(
                MOCK_FILE_2_CONTENT  # pyright: ignore[reportGeneralTypeIssues]
            )

        reads_retriever = RNAReadsRetriever(
            [tmp_path / "test_file_1.gz", tmp_path / "test_file_2.gz"]
        )
        expected_ids = [12, 21, 31, 32, 41, 52]
        expected_reads = [
            (
                "GAGCTCGATTTGCAATGAGCCGCCCCTCTTCCATAAAATCTTGCAAATTC"
                "CGCGTCTCGGGCTCTGTCCAACACGTATGCCTCCCTTCACGATGACTTAAG"
            ),
            (
                "ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCT"
                "ACGAGGCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC"
            ),
            (
                "GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGT"
                "TTTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC"
            ),
            (
                "AGCACCTGGTGCATTACTTCCTACCAACATGGATACAAGATGGGCCTTGCG"
                "CTCTTTAAGTCCGCTGAGATTGCTCTATTCACTAAATCGTACGTACTGGA"
            ),
            (
                "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGT"
                "ACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
            ),
            (
                "GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTT"
                "TTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC"
            ),
        ]
        assert reads_retriever.get_read_sequences_for_ids([12, 21, 31, 32, 41, 52]) == (
            expected_ids,
            expected_reads,
        )

    def test_cutoff(self, tmp_path):
        with open(tmp_path / "test_file_1.fq", "wt", encoding="utf-8") as test_file:
            test_file.write(MOCK_FILE_1_CONTENT)
        reads_retriever = RNAReadsRetriever([tmp_path / "test_file_1.fq"], cutoff=10)
        expected_ids = [11, 41, 51]
        expected_reads = [
            "GCGTGTAATG",
            "GCGTGTAATG",
            "ACCGCCACGC",
        ]
        assert reads_retriever.get_read_sequences_for_ids([11, 41, 51]) == (
            expected_ids,
            expected_reads,
        )
