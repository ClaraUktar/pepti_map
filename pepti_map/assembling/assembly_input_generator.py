from pathlib import Path
from typing import List, Tuple


class AssemblyInputGenerator:
    @staticmethod
    def write_fasta_with_sequences(
        sequences: List[Tuple[int, str]], filepath: Path
    ) -> None:
        with open(filepath, "wt", encoding="utf-8") as fasta_file:
            for sequence_id, sequence in sequences:
                fasta_file.write(f">{str(sequence_id)}\n")
                fasta_file.write(sequence + "\n")
