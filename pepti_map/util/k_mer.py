from typing import Generator, Tuple


def split_into_kmer(
    sequence: str, kmer_length: int = 7
) -> Generator[Tuple[str, int], None, None]:
    for i in range(0, len(sequence) - kmer_length + 1):
        kmer = sequence[i : i + kmer_length]  # noqa: E203
        if "*" in kmer:
            continue
        yield (kmer, i)
