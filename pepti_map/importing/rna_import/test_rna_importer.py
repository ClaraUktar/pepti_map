from unittest.mock import patch
from io import StringIO
import pandas as pd

import pytest

from pepti_map.importing.rna_import.rna_importer import RNAImporter


class TestRNAImporter:
    rna_importer = RNAImporter()

    # Test data was randomly generated
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
    ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAACAGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC
    +
    <zpLIB(|@OT;0R`<5n5hbT5<|kSBZo!5kI1)1aL~7qX3\\2MIbj*/'+}pk)kVQdlcK1\\Fe|M68b2`2Y/Zz}"NUhOLzay097Y.$'`]/
    """

    mock_file_2_content = """@ABC1234567.1 1/2
    CTTAAGTCATCGTGAAGGGAGGCATACGTGTTGGACAGAGCCCGAGACGCGGAATTTGCAAGATTTTATGGAAGAGGGGCGGCTCATTGCAAATCGAGCTC
    +
    mp^]Ro,hU^,R1[\\)2BRU)9,9ksW"UCqAHjE3_W#0*9:`!X.GkWpO{*1^+MLkLRz<gh6c_#db*Olo9j.!2BhNPbjPhVC)l9]hSs#E8
    @ABC1234567.2 2/2
    AGGCTGCACGGTTTCCGCCCAGGGGATCCGCTTTCGGGTATTTCTCTCACGGTATTTCTCCAAAAAAGAAAATTCTCGTCAGACGCTGGTCAGGATCCAGA
    +
    goP*&{}8<FDW\\Hk]K2{eAAF{%%:nP}6Woh@"gPqW8`kojH2CciLy0*Vxix7&Zj3V1D"O[A&Uv-v(q,6n-$1dZ6E)#wXiF}!dp")W'
    @ABC1234567.3 3/2
    TCCAGTACGTACGATTTAGTGAATAGAGCAATCTCAGCGGACTTAAAGAGCGCAAGGCCCATCTTGTATCCATGTTGGTAGGAAGTAATGCACCAGGTGCT
    +
    52tSeIApSGDrD.gjaVT!'xC([GIVHybb@'OKIc9#&4u%Bp*w;H~JCcm:s,v[vddko@`T'F1,D5X,\\ue.xL6JY=qX>4)Fu?M4+DU7@
    @ABC1234567.4 4/2
    CACCGCGAATTACAATGCTTGTCTGGGGCAGACGCATTACCAGCTGCNAGAGTCCAGAGTTAAGTTGTGTACCGTTGCCCGTTGGTAGAACTCGCACAAGC
    +
    {~I&KDf`}L$.=XLjiI)]o/-t5nn$]?mOEB&+`+3-QqMwrG/Y5Jt}X!Gz+b``KzJ!0eJ7#4Y;vnT(>bMfLrwRrlK%gDzh"",0XkU*_
    @ABC1234567.5 5/2
    GATAGGAGAGGATGAGTGCGCTTACTCCGGGGTGTTGCTATAAGATGAAAACAGGGGTATAAAACAGTTGGGAGCGTACTCCCGCGCGTATCTCGAGGGCC
    +
    .6qJoOJ-^"e~@U=ZU#gZIYj^Mlp{r)X*$+]"{mUbGl16N(Z,:=4P!09=Z3s5w<tj7V|(9l7,7-E-\\nZ$IIGt*]:)j5R@>/7quD[2O
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
                    "ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAAC"
                    "AGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
                ),
            ],
            "count": [2, 1, 1, 1],
        }
    ).astype({"ids": "string", "sequence": "string", "count": "uint32"})

    expected_result_df_paired_end = pd.DataFrame(
        {
            "ids": [
                "@ABC1234567.1 1/1,@ABC1234567.4 4/1",
                "@ABC1234567.2 2/1",
                "@ABC1234567.3 3/1,@ABC1234567.5 5/2",
                "@ABC1234567.5 5/1",
                "@ABC1234567.1 1/2",
                "@ABC1234567.2 2/2",
                "@ABC1234567.3 3/2",
                "@ABC1234567.4 4/2",
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
                    "ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAAC"
                    "AGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
                ),
                (
                    "GAGCTCGATTTGCAATGAGCCGCCCCTCTTCCATAAAATCTTGCAAATTC"
                    "CGCGTCTCGGGCTCTGTCCAACACGTATGCCTCCCTTCACGATGACTTAAG"
                ),
                (
                    "TCTGGATCCTGACCAGCGTCTGACGAGAATTTTCTTTTTTGGAGAAATAC"
                    "CGTGAGAGAAATACCCGAAAGCGGATCCCCTGGGCGGAAACCGTGCAGCCT"
                ),
                (
                    "AGCACCTGGTGCATTACTTCCTACCAACATGGATACAAGATGGGCCTTGCG"
                    "CTCTTTAAGTCCGCTGAGATTGCTCTATTCACTAAATCGTACGTACTGGA"
                ),
                (
                    "GCTTGTGCGAGTTCTACCAACGGGCAACGGTACACAACTTAACTCTGGACTCT"
                    "NGCAGCTGGTAATGCGTCTGCCCCAGACAAGCATTGTAATTCGCGGTG"
                ),
            ],
            "count": [2, 1, 2, 1, 1, 1, 1, 1],
        }
    ).astype({"ids": "string", "sequence": "string", "count": "uint32"})

    expected_result_df_single_end_cutoff = pd.DataFrame(
        {
            "ids": [
                "@ABC1234567.1 1/1,@ABC1234567.4 4/1",
                "@ABC1234567.2 2/1,@ABC1234567.5 5/1",
                "@ABC1234567.3 3/1",
            ],
            "sequence": ["GCGTGTAATG", "ACCGCCACGC", "GGCCCTCGAG"],
            "count": [2, 2, 1],
        }
    ).astype({"ids": "string", "sequence": "string", "count": "uint32"})

    def test_raises_error_when_no_file_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files([])

    def test_raises_error_when_more_than_two_files_given(self):
        with pytest.raises(ValueError):
            self.rna_importer.import_files(["file1", "file2", "file3"])

    @patch("gzip.open", return_value=StringIO(mock_file_1_content))
    def test_import_single_end_file(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file"], should_translate=False
        )
        pd.testing.assert_frame_equal(result_df, self.expected_result_df_single_end)

    @patch(
        "gzip.open",
        side_effect=[
            StringIO(mock_file_1_content),
            StringIO(mock_file_2_content),
        ],
    )
    def test_import_paired_end_file(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file1", "path/to/file2"], should_translate=False
        )
        pd.testing.assert_frame_equal(result_df, self.expected_result_df_paired_end)

    @patch("gzip.open", return_value=StringIO(mock_file_1_content))
    def test_cutoff(self, _):
        result_df = self.rna_importer.import_files(
            ["path/to/file"], cutoff=10, should_translate=False
        )
        pd.testing.assert_frame_equal(
            result_df, self.expected_result_df_single_end_cutoff
        )

    # TODO: Test gzip vs uncompressed?
