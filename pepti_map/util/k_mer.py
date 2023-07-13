def split_into_kmer(sequence: str, kmer_length: int = 6, step: int = 1):
    for i in range(0, len(sequence) - kmer_length + 1, step):
        yield sequence[i : i + kmer_length]  # noqa: E203
