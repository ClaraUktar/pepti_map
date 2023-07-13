from typing import Any, Generator
from Bio.Seq import translate


def translate_for_frame(sequence: str, frame: int) -> str:
    def shorten_sequence_to_translatable_len(sequence: str) -> str:
        mod = len(sequence) % 3
        if mod == 0:
            return sequence
        return sequence[:-mod]

    # TODO: Exchange all I for L?
    translation = translate(shorten_sequence_to_translatable_len(sequence[frame:]))
    assert isinstance(translation, str)
    return translation


def get_three_frame_translations(sequence: str) -> Generator[str, Any, Any]:
    # TODO: Split into k-mers here?
    for i in range(3):
        yield translate_for_frame(sequence, i)
    return [translate_for_frame(sequence, i) for i in range(3)]
