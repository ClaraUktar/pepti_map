from pathlib import Path
from unittest.mock import mock_open, patch
from io import StringIO
import pandas as pd

import pytest

from pepti_map.importing.rna_import.rna_importer import RNAImporter
from pepti_map.importing.rna_import.lazy_rna_reader import LazyRNAReader
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
