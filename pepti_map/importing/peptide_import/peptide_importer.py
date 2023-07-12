from typing import Dict, Tuple

import pandas as pd


# TODO: Add support for custom format
# TODO: Add support for MZTab
def import_file(file_path: str) -> pd.DataFrame:
    """
    Reads the file given and transforms it into a pandas DataFrame,
    with columns `ids`, `sequence`, and `count`.

    :param str file_path: The paths to the file that should be imported.
    :returns A pandas DataFrame. The column `ids` contains all ids with duplicate
    sequences. The column `sequence` contains the corresponding sequence.
    The column `count` contains the number of duplicates for the sequence.
    :rtype pandas.DataFrame
    raises FileNotFoundError: Raised if no file could be found for the given path.
    """
    peptides: Dict[str, Tuple[str, int]] = {}
    with open(file_path) as peptide_file:
        for line in peptide_file:
            # TODO: Exchange all I for L?
            sequence = line.strip()
            duplicate = peptides.get(sequence)
            if duplicate is not None:
                peptides[sequence] = (duplicate[0], duplicate[1] + 1)
            else:
                peptides[sequence] = (sequence, 1)
    peptide_df = pd.DataFrame(list(peptides.values()), columns=["sequence", "count"])
    print(peptide_df)
    return peptide_df
