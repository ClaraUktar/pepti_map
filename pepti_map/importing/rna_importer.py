from typing import Dict, Tuple
import pandas as pd


# TODO: Give option to import two files or deal with one pairend end file, see https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#pairedend-fastq
# TODO: Enable reading of zipped file
# e.g. https://stackoverflow.com/questions/19371860/python-open-file-in-zip-without-temporarily-extracting-it
def import_file(file_path: str) -> pd.DataFrame:
    with open(file_path) as rna_file:
        df_data: Dict[str, Tuple[str, str, int]] = {}
        line_count_for_current_sequence: int = 0
        id = ""
        sequence = ""
        duplicate = None

        for line in rna_file:
            if line_count_for_current_sequence == 0:
                id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()
                duplicate = df_data.get(sequence)
            # Information from field 2 (line 3) is not needed
            # For now skip quality info, getting cutoff value supplied by user
            # elif line_count_for_current_sequence == 3:
            #     quality = line.strip()

            line_count_for_current_sequence = line_count_for_current_sequence + 1

            # Always read 4 lines per sequence, as per FASTQ format
            if line_count_for_current_sequence == 4:
                if duplicate is not None:
                    # TODO: Should use list here instead?
                    # TODO: Deal with quality for duplicates
                    df_data[sequence] = (
                        "".join([duplicate[0], ",", id]),
                        duplicate[1],
                        duplicate[2] + 1,
                    )
                    duplicate = None
                else:
                    df_data[sequence] = (id, sequence, 1)
                line_count_for_current_sequence = 0
                id = ""
                sequence = ""

        rna_df = pd.DataFrame(
            list(df_data.values()), columns=["ids", "sequence_after_cutoff", "count"]
        )
        print(rna_df)
        return rna_df
