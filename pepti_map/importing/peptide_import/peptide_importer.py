from typing import Dict, List, Tuple

import pandas as pd


class PeptideImporter:
    def __init__(self):
        self._peptide_dict: Dict[str, Tuple[List[int], str, int]] = {}

    def reset(self) -> None:
        self._peptide_dict = {}

    # TODO: Add support for custom format
    # TODO: Add support for MZTab
    # TODO: Use numpy instead?
    def import_file_with_deduplication(self, filepath: str) -> pd.DataFrame:
        """
        Reads the file given and transforms it into a pandas DataFrame,
        with columns `ids`, `sequence` and `count`.

        :param str filepath: The paths to the file that should be imported.
        :returns A pandas DataFrame. It has one row per unique sequence in the
        given file. The column `sequence` contains the sequence.
        The column `ids` contains all ids for this sequence.
        In case no ids are present in the file, the index of each line is used as its
        id. The column `count` contains the number of duplicates for the sequence.
        :rtype pandas.DataFrame
        raises FileNotFoundError: Raised if no file could be found for the given path.
        """
        self.reset()

        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for index, line in enumerate(peptide_file):
                # TODO: Exchange all I for L?
                sequence = line.strip()
                duplicate = self._peptide_dict.get(sequence)
                if duplicate is not None:
                    self._peptide_dict[sequence] = (
                        duplicate[0] + [index],
                        duplicate[1],
                        duplicate[2] + 1,
                    )
                else:
                    self._peptide_dict[sequence] = ([index], sequence, 1)
        peptide_df = pd.DataFrame(
            list(self._peptide_dict.values()), columns=["ids", "sequence", "count"]
        ).astype(dtype={"ids": "object", "sequence": "string", "count": "uint32"})
        print(peptide_df)
        print(peptide_df.info(verbose=True))
        return peptide_df

    # TODO: Still save the peptide data in separate class, similar to the rna data?
    def import_file(self, filepath: str) -> pd.DataFrame:
        """
        TODO
        """
        peptides = []
        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for line in peptide_file:
                peptides.append(line.strip())
        peptide_df = pd.DataFrame(peptides, columns=["sequence"]).astype(
            dtype={"sequence": "string"}
        )
        print(peptide_df)
        print(peptide_df.info(verbose=True))
        return peptide_df
