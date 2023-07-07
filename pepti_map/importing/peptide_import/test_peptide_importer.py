from unittest.mock import patch, mock_open
import pandas as pd

from pepti_map.importing.peptide_import import peptide_importer


class TestPeptideSimpleFormatImport:
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

    expected_result_df = pd.DataFrame(
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
                "FMRPCQFFHCWMFSMDDH",
            ],
            "count": [1, 1, 1, 1, 2, 1, 1, 1, 1],
        }
    )

    def test_import_as_df_with_duplicates(self):
        with patch(
            "builtins.open", mock_open(read_data=self.mock_file_content)
        ) as peptide_file_mock:
            result_df = peptide_importer.import_file("path/to/file")
            pd.testing.assert_frame_equal(result_df, self.expected_result_df)
        peptide_file_mock.assert_called_with("path/to/file")
