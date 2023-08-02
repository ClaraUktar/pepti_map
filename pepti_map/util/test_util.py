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

    def test_big_step(self):
        test_sequence = "WHQVRNWCKHVEIEQCLECV"
        expected_result = [("WHQVRNWCKH", 0), ("NWCKHVEIEQ", 5), ("VEIEQCLECV", 10)]
        assert list(split_into_kmer(test_sequence, 10, 5)) == expected_result

    def test_split_by_stop_codons(self):
        test_sequence = "ACNVMILCLF*SARFFG*VLP"
        for kmer in split_into_kmer(test_sequence, 6):
            print(kmer)
        expected_result = [
            ("ACNVMI", 0),
            ("CNVMIL", 1),
            ("NVMILC", 2),
            ("VMILCL", 3),
            ("MILCLF", 4),
            ("SARFFG", 11),
        ]
        assert list(split_into_kmer(test_sequence, 6)) == expected_result
