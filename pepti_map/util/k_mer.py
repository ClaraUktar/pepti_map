from typing import Generator, Tuple


def split_into_kmer(
    sequence: str, kmer_length: int = 6, step: int = 1
) -> Generator[Tuple[str, int], None, None]:
    for i in range(0, len(sequence) - kmer_length + 1, step):
        yield (sequence[i : i + kmer_length], i)  # noqa: E203
