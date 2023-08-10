from unittest.mock import patch, mock_open
import pandas as pd
import pytest

from pepti_map.importing.peptide_import.peptide_importer import PeptideImporter


class TestPeptideSimpleFormatImport:
    # TODO: Extract into separate data file
    # Test data was randomly generated
    mock_file_content = """GQLDR
    NCYQKAQHLYTPEGRKGMVHLTWDRTVLPPPCMDIVDRHRSSRY
    QIPCI
    LSVRCEVQCHWDYDLPVMRHKTFLAPCHTWM
    YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE
    FNQGVRDQR
    ISNVGQYSQMLSEFEFCKKVYCSCEFDVDRGSQICSDDHVWSPKKSPKMC
    SHDTY
    YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE
    FMRPCQFFHCWMFSMDDH"""

    expected_result_df_deduplicated = pd.DataFrame(
        {
            "ids": [[0], [1], [2], [3], [4, 8], [5], [6], [7], [9]],
            "sequence": [
                "GQLDR",
                "NCYQKAQHLYTPEGRKGMVHLTWDRTVLPPPCMDIVDRHRSSRY",
                "QIPCI",
                "LSVRCEVQCHWDYDLPVMRHKTFLAPCHTWM",
                "YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE",
                "FNQGVRDQR",
                "ISNVGQYSQMLSEFEFCKKVYCSCEFDVDRGSQICSDDHVWSPKKSPKMC",
                "SHDTY",
                "FMRPCQFFHCWMFSMDDH",
            ],
            "count": [1, 1, 1, 1, 2, 1, 1, 1, 1],
        }
    ).astype(dtype={"ids": "object", "sequence": "string", "count": "uint32"})

    expected_result_basic = pd.DataFrame(
        {
            "sequence": [
                "GQLDR",
                "NCYQKAQHLYTPEGRKGMVHLTWDRTVLPPPCMDIVDRHRSSRY",
                "QIPCI",
                "LSVRCEVQCHWDYDLPVMRHKTFLAPCHTWM",
                "YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE",
                "FNQGVRDQR",
                "ISNVGQYSQMLSEFEFCKKVYCSCEFDVDRGSQICSDDHVWSPKKSPKMC",
                "SHDTY",
                "YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE",
                "FMRPCQFFHCWMFSMDDH",
            ]
        }
    ).astype(dtype={"sequence": "string"})

    @pytest.fixture(autouse=True)
    def _init_peptide_importer(self):
        self.peptide_importer = PeptideImporter()

    def test_import_as_df_with_duplicates(self):
        with patch(
            "builtins.open", mock_open(read_data=self.mock_file_content)
        ) as peptide_file_mock:
            result_df = self.peptide_importer.import_file_with_deduplication(
                "path/to/file"
            )
            pd.testing.assert_frame_equal(
                result_df, self.expected_result_df_deduplicated
            )
        peptide_file_mock.assert_called_with("path/to/file", "rt", encoding="utf-8")

    @patch("builtins.open", mock_open(read_data=mock_file_content))
    def test_basic_import(self):
        result_df = self.peptide_importer.import_file("path/to/file")
        pd.testing.assert_frame_equal(result_df, self.expected_result_basic)
