from collections import defaultdict
import csv
from pathlib import Path
from typing import List, Set, Union
import numpy as np
import pytest
from pepti_map.importing.peptide_import.testdata_peptide_importer import (
    EXPECTED_PEPTIDE_MAPPING,
    EXPECTED_PEPTIDE_MAPPING_PROTEIN_GROUPS,
    EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED,
    EXPECTED_RESULT_LIST,
)
from pepti_map.matching.precomputing_rna_to_peptide_matcher import (
    PrecomputingRNAToPeptideMatcher,
)

from pepti_map.matching.rna_to_peptide_matcher import (
    PEPTIDE_READ_QUANT_FILENAME,
    RNAToPeptideMatcher,
)
from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex


class TestRNAToPeptideMatcher:
    @pytest.fixture(autouse=True)
    def _init_matcher(self):
        self.kmer_index = PeptideKmerIndex()
        self.kmer_index.kmer_index = defaultdict(
            list, EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED.copy()
        )
        self.matcher = RNAToPeptideMatcher(self.kmer_index, 7, EXPECTED_PEPTIDE_MAPPING)

    def test_add_single_match(self):
        EXPECTED_MATCHING_RESULT: List[Union[Set[int], None]] = [
            None for _ in range(0, 7)
        ]
        EXPECTED_MATCHING_RESULT[2] = set([1])
        EXPECTED_MATCHING_RESULT[3] = set([1])
        EXPECTED_MATCHING_RESULT[5] = set([1])

        self.matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        assert self.matcher.get_matches() == EXPECTED_MATCHING_RESULT

    def test_add_multiple_matches(self):
        EXPECTED_MATCHING_RESULT: List[Union[Set[int], None]] = [
            None for _ in range(0, 7)
        ]

        EXPECTED_MATCHING_RESULT[0] = set([2])
        EXPECTED_MATCHING_RESULT[2] = set([1, 2])
        EXPECTED_MATCHING_RESULT[3] = set([1])
        EXPECTED_MATCHING_RESULT[5] = set([1, 2])

        self.matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        self.matcher.add_peptide_matches_for_rna_read(
            2,
            "GATGTAAGTTGATATCGTAGACCGGCATCGCAAGGCACAACCCTGGCGTGAACCGA",
        )
        assert self.matcher.get_matches() == EXPECTED_MATCHING_RESULT

    def test_peptide_quant_file(self, tmp_path):
        # TODO: Move to data file?
        EXPECTED_FILE_CONTENTS = [
            ["peptide_sequence", "group_id", "n_reads_peptide", "n_reads_group"],
            ["GQLDR", "-1", "0", "0"],
            ["NCYQKAQHLYTPEGRKGMVHLTWDRTVLPPPCMDIVDRHRSSRY", "0", "1", "1"],
            ["QIPCI", "-1", "0", "0"],
            ["LSVRCEVQCHWDYDLPVMRHKTFLAPCHTWM", "1", "0", "0"],
            ["YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE", "2", "2", "2"],
            ["FNQGVRDQR", "3", "1", "1"],
            ["ISNVGQYSQMLSEFEFCKKVYCSCEFDVDRGSQICSDDHVWSPKKSPKMC", "4", "0", "0"],
            ["SHDTY", "-1", "0", "0"],
            ["YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE", "5", "2", "2"],
            ["FMRPCQFFHCWMFSMDDH", "6", "0", "0"],
        ]

        self.matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        self.matcher.add_peptide_matches_for_rna_read(
            2,
            "GATGTAAGTTGATATCGTAGACCGGCATCGCAAGGCACAACCCTGGCGTGAACCGA",
        )
        temp_directory: Path = tmp_path / "peptide_quant"
        temp_directory.mkdir()
        self.matcher.write_peptide_read_quant_file(temp_directory, EXPECTED_RESULT_LIST)
        file_entries = []
        with open(
            temp_directory / PEPTIDE_READ_QUANT_FILENAME, "rt", encoding="utf-8"
        ) as peptide_quant_file:
            reader = csv.reader(peptide_quant_file, delimiter="\t", lineterminator="\n")
            for row in reader:
                file_entries.append(row)
        assert file_entries == EXPECTED_FILE_CONTENTS

    def test_peptide_quant_file_groups(self, tmp_path):
        group_matcher = RNAToPeptideMatcher(
            self.kmer_index, 5, EXPECTED_PEPTIDE_MAPPING_PROTEIN_GROUPS
        )
        # TODO: Move to data file?
        EXPECTED_FILE_CONTENTS = [
            ["peptide_sequence", "group_id", "n_reads_peptide", "n_reads_group"],
            ["GQLDR", "-1", "0", "0"],
            ["NCYQKAQHLYTPEGRKGMVHLTWDRTVLPPPCMDIVDRHRSSRY", "0", "1", "1"],
            ["QIPCI", "-1", "0", "0"],
            ["LSVRCEVQCHWDYDLPVMRHKTFLAPCHTWM", "1", "0", "0"],
            ["YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE", "2", "2", "2"],
            ["FNQGVRDQR", "2", "1", "2"],
            ["ISNVGQYSQMLSEFEFCKKVYCSCEFDVDRGSQICSDDHVWSPKKSPKMC", "3", "0", "0"],
            ["SHDTY", "-1", "0", "0"],
            ["YCYRSEDLLKISQQCARHRKAQPWETYVCIRTFTPHTMMHE", "4", "2", "2"],
            ["FMRPCQFFHCWMFSMDDH", "0", "0", "1"],
        ]

        group_matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        group_matcher.add_peptide_matches_for_rna_read(
            2,
            "GATGTAAGTTGATATCGTAGACCGGCATCGCAAGGCACAACCCTGGCGTGAACCGA",
        )
        temp_directory: Path = tmp_path / "peptide_quant"
        temp_directory.mkdir()
        group_matcher.write_peptide_read_quant_file(
            temp_directory, EXPECTED_RESULT_LIST
        )
        file_entries = []
        with open(
            temp_directory / PEPTIDE_READ_QUANT_FILENAME, "rt", encoding="utf-8"
        ) as peptide_quant_file:
            reader = csv.reader(peptide_quant_file, delimiter="\t", lineterminator="\n")
            for row in reader:
                file_entries.append(row)
        assert file_entries == EXPECTED_FILE_CONTENTS


class TestPrecomputingRNAToPeptideMatcher:
    @pytest.fixture(autouse=True)
    def _init_matcher(self):
        self.kmer_index = PeptideKmerIndex()
        self.kmer_index.kmer_index = defaultdict(
            list, EXPECTED_RESULT_INDEX_ISOLEUCINE_REPLACED.copy()
        )
        self.matcher = PrecomputingRNAToPeptideMatcher(
            self.kmer_index, 7, EXPECTED_PEPTIDE_MAPPING
        )

    def test_add_single_match_with_precomputation(self):
        EXPECTED_MATCHING_RESULT: List[Union[Set[int], None]] = [
            None for _ in range(0, 7)
        ]
        EXPECTED_MATCHING_RESULT[2] = set([1])
        EXPECTED_MATCHING_RESULT[3] = set([1])
        EXPECTED_MATCHING_RESULT[5] = set([1])
        expected_precomputed_intersections = np.array(
            [
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 0, 1, 0],
                [0, 0, 1, 1, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [0, 0, 1, 1, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=np.uint32,
        )

        self.matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        assert self.matcher.get_matches() == EXPECTED_MATCHING_RESULT
        np.testing.assert_array_equal(
            self.matcher.get_precomputed_intersections(),
            expected_precomputed_intersections,
        )

    def test_add_multiple_matches_wth_precomputation(self):
        EXPECTED_MATCHING_RESULT: List[Union[Set[int], None]] = [
            None for _ in range(0, 7)
        ]

        EXPECTED_MATCHING_RESULT[0] = set([2])
        EXPECTED_MATCHING_RESULT[2] = set([1, 2])
        EXPECTED_MATCHING_RESULT[3] = set([1])
        EXPECTED_MATCHING_RESULT[5] = set([1, 2])
        expected_precomputed_intersections = np.array(
            [
                [1, 0, 1, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [1, 0, 2, 1, 0, 2, 0],
                [0, 0, 1, 1, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0],
                [1, 0, 2, 1, 0, 2, 0],
                [0, 0, 0, 0, 0, 0, 0],
            ],
            dtype=np.uint32,
        )

        self.matcher.add_peptide_matches_for_rna_read(
            1,
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        self.matcher.add_peptide_matches_for_rna_read(
            2,
            "GATGTAAGTTGATATCGTAGACCGGCATCGCAAGGCACAACCCTGGCGTGAACCGA",
        )
        assert self.matcher.get_matches() == EXPECTED_MATCHING_RESULT
        np.testing.assert_array_equal(
            self.matcher.get_precomputed_intersections(),
            expected_precomputed_intersections,
        )
