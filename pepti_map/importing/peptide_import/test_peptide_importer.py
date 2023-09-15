from pathlib import Path
from unittest.mock import patch, mock_open
import pandas as pd
import pytest

from pepti_map.importing.peptide_import.peptide_importer import PeptideImporter
from pepti_map.importing.peptide_import.peptide_to_index_importer import (
    PeptideToIndexImporter,
)
from pepti_map.importing.peptide_import.testdata_peptide_importer import (
    EXPECTED_PEPTIDE_MAPPING,
    EXPECTED_PEPTIDE_MAPPING_PROTEIN_GROUPS,
    EXPECTED_RESULT_DF_DEDUPLICATED,
    EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED,
    EXPECTED_RESULT_LIST,
    MOCK_FILE_CONTENT,
    PROTEIN_GROUPS_MOCK_FILE_CONTENT,
)


class TestPeptideSimpleFormatImport:
    @pytest.fixture(autouse=True)
    def _init_peptide_importer(self):
        self.peptide_importer = PeptideImporter()

    def test_import_as_df_with_duplicates(self):
        with patch(
            "builtins.open", mock_open(read_data=MOCK_FILE_CONTENT)
        ) as peptide_file_mock:
            result_df = self.peptide_importer.import_file_with_deduplication(
                "path/to/file"
            )
            pd.testing.assert_frame_equal(result_df, EXPECTED_RESULT_DF_DEDUPLICATED)
        peptide_file_mock.assert_called_with("path/to/file", "rt", encoding="utf-8")

    @patch("builtins.open", mock_open(read_data=MOCK_FILE_CONTENT))
    def test_basic_import(self):
        result_list = self.peptide_importer.import_file(Path("path/to/file"))
        assert result_list == EXPECTED_RESULT_LIST


class TestPeptideToIndexImporter:
    @patch("builtins.open", mock_open(read_data=MOCK_FILE_CONTENT))
    @patch(
        (
            "pepti_map.importing.peptide_import.peptide_to_index_importer"
            ".PeptideToIndexImporter._write_peptide_to_cluster_mapping_file"
        )
    )
    def test_import_file_to_index(self, mock_write_peptide_to_cluster_mapping_file):
        resulting_index, _ = PeptideToIndexImporter().import_file_to_index(
            Path("path/to/file")
        )
        assert resulting_index.kmer_index == EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED
        assert resulting_index.number_of_peptides == 7
        mock_write_peptide_to_cluster_mapping_file.assert_called_once_with(
            EXPECTED_PEPTIDE_MAPPING
        )

    @patch("builtins.open", mock_open(read_data=PROTEIN_GROUPS_MOCK_FILE_CONTENT))
    @patch(
        (
            "pepti_map.importing.peptide_import.peptide_to_index_importer"
            ".PeptideToIndexImporter._write_peptide_to_cluster_mapping_file"
        )
    )
    def test_import_file_to_index_with_protein_groups(
        self, mock_write_peptide_to_cluster_mapping_file
    ):
        resulting_index, _ = PeptideToIndexImporter().import_file_to_index(
            Path("path/to/file")
        )
        assert resulting_index.kmer_index == EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED
        assert resulting_index.number_of_peptides == 5
        mock_write_peptide_to_cluster_mapping_file.assert_called_once_with(
            EXPECTED_PEPTIDE_MAPPING_PROTEIN_GROUPS
        )
