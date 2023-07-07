"""
- test single file import
- test import of two files
- DONE test error for incorrect number of files given
- test cutoff
- test gzip vs uncompressed?
"""


from unittest.mock import patch
import pandas as pd
import io

import pytest

from pepti_map.importing.rna_import import rna_importer


class TestRNAImporter:
    mock_file_1_content = """@ABC1234567.1 1/1
    GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA
    +
    A*}?_]?+i%(Bh48IfmITF;Wf!@i:q.$$G~H[]lo)d5/Xl;m"I*UDY@uy]_{]2}QFM%2uAd4-K>2{tF,.kVWjcE-f*f?57?)h_flau
    @ABC1234567.2 2/1
    ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCTACGAGGCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC
    +
    EQn;SjL`t$a]&GE%S).<9-Hl!#xtwmup/r!LJPP%2b_@SXtaHc'6m>6AmeR#pP~P$;FL1m0XYvI)zKr%2-XlWO\\Y8c.OW!v2./N^\\
    @ABC1234567.3 3/1
    GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTTTTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC
    +
    h=[U[nDhr,8VmiVsx[F=gZ>UIk5>^4Z/0=N}%vA]U$|{byQAb}:YUH~FJn_92*-Uv$mk*kg,p}SV3>zT/pIK3o-)v[B/%3XcSH&'<
    @ABC1234567.4 4/1
    GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA
    +
    '~~3QhIkE}yGC<jzj.Lx1|,oV/{dQMZ/ksH>/iH5uK`3c;K]=a-O/N{3.G+ym=\\_8!-;A`JBA:~-tYt7l]5J3%ZMlf\\},a"NcE}I$
    @ABC1234567.5 5/1
    CATACAAGGAATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAACAGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC
    +
    <zpLIB(|@OT;0R`<5n5hbT5<|kSBZo!5kI1)1aL~7qX3\\2MIbj*/'+}pk)kVQdlcK1\\Fe|M68b2`2Y/Zz}"NUhOLzay097Y.$'`]/
    """

    mock_file_2_content = """
    """

    expected_result_df_single_end = pd.DataFrame(
        {
            "ids": [
                "@ABC1234567.1 1/1,@ABC1234567.4 4/1",
                "@ABC1234567.2 2/1",
                "@ABC1234567.3 3/1",
                "@ABC1234567.5 5/1",
            ],
            "sequence": [
                (
                    "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTG"
                    "GTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
                ),
                (
                    "ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCTACGAG"
                    "GCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC"
                ),
                (
                    "GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTT"
                    "TTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC"
                ),
                (
                    "CATACAAGGAATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAAC"
                    "AGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
                ),
            ],
            "count": [2, 1, 1, 1],
        }
    )

    expected_result_df_paired_end = pd.DataFrame({})

    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            rna_importer.import_file([])

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            rna_importer.import_file(["file1", "file2", "file3"])

    @patch("gzip.open", return_value=io.StringIO(mock_file_1_content))
    def test_import_single_end_file(self, _):
        result_df = rna_importer.import_file(["path/to/file"])
        pd.testing.assert_frame_equal(result_df, self.expected_result_df_single_end)

    # @Mock(
    #     "gzip",
    #     side_effect=[str.encode(mock_file_1_content), str.encode(mock_file_2_content)],
    # )
    # def test_import_paired_end_file(self):
    #     result_df = rna_importer.import_file(["path/to/file1", "path/to/file2"])
    #     pd.testing.assert_frame_equal(result_df, self.expected_result_df_paired_end)
