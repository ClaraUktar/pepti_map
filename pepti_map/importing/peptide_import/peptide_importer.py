from typing import Dict, Tuple

import pandas as pd


class PeptideImporter:
    _peptide_dict: Dict[str, Tuple[str, int]] = {}

    def reset(self) -> None:
        self._peptide_dict = {}

    # TODO: Add support for custom format
    # TODO: Add support for MZTab
    # TODO: Use numpy instead?
    def import_file(self, file_path: str) -> pd.DataFrame:
        """
        Reads the file given and transforms it into a pandas DataFrame,
        with columns `sequence` and `count`.

        :param str file_path: The paths to the file that should be imported.
        :returns A pandas DataFrame. It has one row per unique sequence in the
        given file. The column `sequence` contains the sequence.
        The column `count` contains the number of duplicates for the sequence.
        :rtype pandas.DataFrame
        raises FileNotFoundError: Raised if no file could be found for the given path.
        """
        self.reset()

        with open(file_path) as peptide_file:
            for line in peptide_file:
                # TODO: Exchange all I for L?
                sequence = line.strip()
                duplicate = self._peptide_dict.get(sequence)
                if duplicate is not None:
                    self._peptide_dict[sequence] = (duplicate[0], duplicate[1] + 1)
                else:
                    self._peptide_dict[sequence] = (sequence, 1)
        peptide_df = pd.DataFrame(
            list(self._peptide_dict.values()), columns=["sequence", "count"]
        ).astype({"sequence": "string", "count": "uint32"})
        print(peptide_df)
        print(peptide_df.info(verbose=True))
        return peptide_df
