from collections import defaultdict
import pytest
from pepti_map.importing.peptide_import.testdata_peptide_importer import (
    EXPECTED_RESULT_INDEX,
)

from pepti_map.matching.rna_to_peptide_matcher import RNAToPeptideMatcher
from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex


class TestRNAToPeptideMatcher:
    @pytest.fixture(autouse=True)
    def _init_matcher(self):
        self.kmer_index = PeptideKmerIndex()
        self.kmer_index.kmer_index = defaultdict(list, EXPECTED_RESULT_INDEX.copy())
        self.matcher = RNAToPeptideMatcher(self.kmer_index, 10)

    def test_add_single_match(self):
        EXPECTED_MATCHING_RESULT = [set() for _ in range(0, 10)]
        EXPECTED_MATCHING_RESULT[4].add("read1")
        EXPECTED_MATCHING_RESULT[5].add("read1")
        EXPECTED_MATCHING_RESULT[8].add("read1")

        self.matcher.add_peptide_matches_for_rna_read(
            "read1",
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        assert self.matcher.matches == EXPECTED_MATCHING_RESULT

    def test_add_multiple_matches(self):
        EXPECTED_MATCHING_RESULT = [set() for _ in range(0, 10)]
        # TODO Finish test
        self.matcher.add_peptide_matches_for_rna_read(
            "read1",
            "AGCTTTCACGCCGCATACGATGATGCACGAATTTAATCAGGGGGTCCGAGATCCAG",
        )
        self.matcher.add_peptide_matches_for_rna_read(
            "read2",
            "",
        )
        assert self.matcher.matches == EXPECTED_MATCHING_RESULT