from pepti_map.util.three_frame_translation import get_three_frame_translations


class TestThreeFrameTranslation:
    def test_three_frame_translation_dna(self):
        test_sequence = "ATGTCCGACGGGACTTGAC"
        expected_result = ["MSDGT*", "CPTGLD", "VRRDL"]
        assert get_three_frame_translations(test_sequence) == expected_result
