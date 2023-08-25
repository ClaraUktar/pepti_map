from pepti_map.util.jaccard_index import jaccard_index
from pepti_map.util.k_mer import split_into_kmer
from pepti_map.util.three_frame_translation import get_three_frame_translations


class TestThreeFrameTranslation:
    def test_three_frame_translation_dna(self):
        test_sequence = "ATGTCCGACGGGACTTGAC"
        expected_result = [("MSDGT*", 0), ("CPTGLD", 1), ("VRRDL", 2)]
        assert list(get_three_frame_translations(test_sequence)) == expected_result


class TestKmerSplitting:
    def test_split_into_kmers(self):
        test_sequence = "WHQVRNWCKHVEIEQCLECV"
        expected_result = [
            ("WHQVRNWCKH", 0),
            ("HQVRNWCKHV", 1),
            ("QVRNWCKHVE", 2),
            ("VRNWCKHVEI", 3),
            ("RNWCKHVEIE", 4),
            ("NWCKHVEIEQ", 5),
            ("WCKHVEIEQC", 6),
            ("CKHVEIEQCL", 7),
            ("KHVEIEQCLE", 8),
            ("HVEIEQCLEC", 9),
            ("VEIEQCLECV", 10),
        ]
        assert list(split_into_kmer(test_sequence, 10)) == expected_result

    def test_split_by_stop_codons(self):
        test_sequence = "ACNVMILCLF*SARFFG*VLP"
        expected_result = [
            ("ACNVMI", 0),
            ("CNVMIL", 1),
            ("NVMILC", 2),
            ("VMILCL", 3),
            ("MILCLF", 4),
            ("SARFFG", 11),
        ]
        assert list(split_into_kmer(test_sequence, 6)) == expected_result


class TestJaccardIndex:
    def test_jaccard_index_empty_set(self):
        set1 = set()
        set2 = {1, 2, 3}
        assert jaccard_index(set1, set2) == 0.0
        assert jaccard_index(set2, set1) == 0.0

    def test_jaccard_index_overlap(self):
        set1 = {1, 2, 3, 4, 5}
        set2 = {3, 4, 5, 6}
        set3 = {3, 4, 5, 6, 7}
        set4 = {1}
        assert jaccard_index(set1, set2) == 0.5
        assert jaccard_index(set2, set3) == 4 / 5
        assert jaccard_index(set1, set4) == 1 / 5

    def test_jaccard_index_identical_sets(self):
        set1 = {2, 4, 6, 8}
        set2 = {8, 6, 4, 2}
        assert jaccard_index(set1, set2) == 1.0

    def test_jaccard_index_no_overlap(self):
        set1 = {1, 3, 5, 7}
        set2 = {2, 4, 6, 8}
        assert jaccard_index(set1, set2) == 0.0
