from typing import Generator, Tuple
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


def get_three_frame_translations(
    sequence: str,
) -> Generator[Tuple[str, int], None, None]:
    """
    For the given nucleic acid sequence, performs a 3-frame translation
    into amino acid sequences.
    :param str sequence: The nucleic acid sequence to be translated.
    :returns A Generator yielding each translated sequence as a Tuple
    (translated_sequence: str, frame: int), where the frame is either 0, 1, or 2,
    depending on the position of the nucleic acid used as the start of the translation.
    :rtype Generator[Tuple[str, int], None, None]
    """
    for i in range(3):
        yield (translate_for_frame(sequence, i), i)
