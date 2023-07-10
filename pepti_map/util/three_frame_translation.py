from typing import Dict
from Bio.Seq import translate


def get_three_frame_translations(sequence: str) -> Dict[int, str]:
    # TODO: Split into k-mers here?
    translated_frames = {}
    seq_length = len(sequence)
    for i in range(3):
        translatable_length = 3 * ((seq_length - i) // 3)
        translated_frames[i] = translate(
            sequence[i : i + translatable_length]  # noqa: E203
        )
    return translated_frames
