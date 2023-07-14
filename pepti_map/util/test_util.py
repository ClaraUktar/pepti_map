from pepti_map.util.k_mer import split_into_kmer
from pepti_map.util.three_frame_translation import get_three_frame_translations


class TestThreeFrameTranslation:
    def test_three_frame_translation_dna(self):
        test_sequence = "ATGTCCGACGGGACTTGAC"
        expected_result = ["MSDGT*", "CPTGLD", "VRRDL"]
        assert list(get_three_frame_translations(test_sequence)) == expected_result


class TestKmerSplitting:
    def test_splitting_into_kmers(self):
        test_sequence = "WHQVRNWCKHVEIEQCLECV"
        expected_result = [
            "WHQVRNWCKH",
            "HQVRNWCKHV",
            "QVRNWCKHVE",
            "VRNWCKHVEI",
            "RNWCKHVEIE",
            "NWCKHVEIEQ",
            "WCKHVEIEQC",
            "CKHVEIEQCL",
            "KHVEIEQCLE",
            "HVEIEQCLEC",
            "VEIEQCLECV",
        ]
        assert list(split_into_kmer(test_sequence, 10)) == expected_result

    def test_big_step(self):
        test_sequence = "WHQVRNWCKHVEIEQCLECV"
        expected_result = ["WHQVRNWCKH", "NWCKHVEIEQ", "VEIEQCLECV"]
        assert list(split_into_kmer(test_sequence, 10, 5)) == expected_result
