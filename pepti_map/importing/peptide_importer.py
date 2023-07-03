from typing import List


# TODO: Add support for custom format
def import_file(file_path: str) -> List[str]:
    peptides: List[str] = []
    with open(file_path) as peptide_file:
        for line in peptide_file:
            peptides.append(line.strip())
    return peptides
