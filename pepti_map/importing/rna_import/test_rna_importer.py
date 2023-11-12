import gzip
from pathlib import Path
from unittest.mock import mock_open, patch
from io import StringIO
import pandas as pd

import pytest

from pepti_map.importing.rna_import.lazy_rna_reader import LazyRNAReader
from pepti_map.importing.rna_import.rna_reads_retriever import RNAReadsRetriever
from pepti_map.importing.rna_import.testdata_rna_importer import (
    EXPECTED_RESULT_LINES_PAIRED_END,
    EXPECTED_RESULT_LINES_SINGLE_END,
    EXPECTED_RESULT_LINES_SINGLE_END_CUTOFF,
    MOCK_FILE_1_CONTENT,
    MOCK_FILE_2_CONTENT,
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
