import pandas as pd


# TODO: Give option to import multiple files
def import_file(file_path: str) -> pd.DataFrame:
    with open(file_path) as rna_file:
        df_data = []
        line_count_for_current_sequence: int = 0
        id = ""
        sequence = ""
        quality = ""

        for line in rna_file:
            if line_count_for_current_sequence == 0:
                id = line.strip()
            elif line_count_for_current_sequence == 1:
                sequence = line.strip()
            # Information from field 3 is not needed
            elif line_count_for_current_sequence == 3:
                quality = line.strip()

            line_count_for_current_sequence = line_count_for_current_sequence + 1

            # Always read 4 lines per sequence, as per FASTQ format
            if line_count_for_current_sequence == 4:
                df_data.append([id, sequence, quality])
                line_count_for_current_sequence = 0
                id = ""
                sequence = ""
                quality = ""

        rna_df = pd.DataFrame(df_data, columns=["id", "sequence", "quality"])
        return rna_df
