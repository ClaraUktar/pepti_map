from pathlib import Path
from typing import Iterable, Tuple


class AssemblyHelper:
    @staticmethod
    def write_fasta_with_sequences(
        sequences: Iterable[Tuple[str, str]], filepath: Path
    ) -> None:
        with open(filepath, "wt", encoding="utf-8") as fasta_file:
            for sequence_id, sequence in sequences:
                fasta_file.write(f">{sequence_id}\n")
                fasta_file.write(sequence + "\n")
