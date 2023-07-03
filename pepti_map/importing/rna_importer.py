from typing import Dict, Tuple
import pandas as pd


# TODO: Give option to import two files or deal with one pairend end file, see https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#pairedend-fastq
# TODO: Enable reading of zipped file
# e.g. https://stackoverflow.com/questions/19371860/python-open-file-in-zip-without-temporarily-extracting-it
def import_file(file_path: str) -> pd.DataFrame:
    with open(file_path) as rna_file:
        df_data: Dict[str, Tuple[str, str, str, int]] = {}
        line_count_for_current_sequence: int = 0
        id = ""
        sequence = ""
        quality = ""
        found_duplicate = False

        for line in rna_file:
            if line_count_for_current_sequence == 0:
                id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()
                duplicate = df_data.get(sequence)
                if duplicate is not None:
                    found_duplicate = True
                    # TODO: Should use list here instead?
                    # TODO: Deal with quality for duplicates
                    df_data[sequence] = (
                        duplicate[0],
                        duplicate[1],
                        duplicate[2],
                        duplicate[3] + 1,
                    )
            # Information from field 3 is not needed
            elif line_count_for_current_sequence == 3:
                quality = line.strip()

            line_count_for_current_sequence = line_count_for_current_sequence + 1

            # Always read 4 lines per sequence, as per FASTQ format
            if line_count_for_current_sequence == 4:
                if found_duplicate:
                    found_duplicate = False
                else:
                    df_data[sequence] = (id, sequence, quality, 1)
                line_count_for_current_sequence = 0
                id = ""
                sequence = ""
                quality = ""

        rna_df = pd.DataFrame(
            list(df_data.values()), columns=["id", "sequence", "quality", "count"]
        )
        print(rna_df)
        return rna_df
