from typing import Dict, Tuple

import pandas as pd


# TODO: Add support for custom format
# TODO: Add support for MZTab
def import_file(file_path: str) -> pd.DataFrame:
    peptides: Dict[str, Tuple[str, int]] = {}
    with open(file_path) as peptide_file:
        for line in peptide_file:
            sequence = line.strip()
            duplicate = peptides.get(sequence)
            if duplicate is not None:
                peptides[sequence] = (duplicate[0], duplicate[1] + 1)
            else:
                peptides[sequence] = (sequence, 1)
    peptide_df = pd.DataFrame(list(peptides.values()), columns=["sequence", "count"])
    print(peptide_df)
    return peptide_df