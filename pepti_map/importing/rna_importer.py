from typing import Dict, TextIO, Tuple
import pandas as pd
import gzip


def _convert_file_to_df(rna_file: TextIO, cutoff: int = -1) -> pd.DataFrame:
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
            if cutoff > 0:
                sequence = sequence[0:cutoff]
            duplicate = df_data.get(sequence)
        # Information from field 2 (line 3) is not needed
        # For now skip quality info, getting cutoff value supplied by user
        # elif line_count_for_current_sequence == 3:
        #     quality = line.strip()

        line_count_for_current_sequence = line_count_for_current_sequence + 1

        # Always read 4 lines per sequence, as per FASTQ format
        if line_count_for_current_sequence == 4:
            if duplicate is not None:
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


# TODO: Give option to import two files or deal with one pairend end file,
# see https://www.ncbi.nlm.nih.gov/sra/docs/submitformats/#pairedend-fastq
# TODO: Add parameter documentation
# cutoff is the position of the last bp in the reads
# after which the cutoff should be performed (starting with 1) (should be > 0)
def import_file(file_path: str, cutoff: int = -1) -> pd.DataFrame:
    with gzip.open(file_path, "rt") as rna_file_gzipped:
        try:
            rna_file_gzipped.read(1)
            print("Detected gzip file. Reading in compressed format...")
            return _convert_file_to_df(rna_file_gzipped, cutoff)
        except gzip.BadGzipFile:
            print("File is not a gzip file. Trying to read as uncompressed file...")

    with open(file_path, "rt") as rna_file:
        return _convert_file_to_df(rna_file, cutoff)
