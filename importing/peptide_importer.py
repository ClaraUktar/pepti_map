from typing import List


def import_file(file_path: str) -> List[str]:
    peptides: List[str] = []
    with open(file_path) as peptide_file:
        for line in peptide_file:
            peptides.append(line.rstrip())
    return peptides
