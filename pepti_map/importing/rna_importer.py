from typing import Dict, List, TextIO, Tuple
import pandas as pd
import gzip

from pepti_map.util.complementarity import BASE_COMPLEMENT


def _convert_rna_data_to_df(
    rna_data: TextIO,
    rna_dict: Dict[str, Tuple[str, str, int]] = {},
    cutoff: int = -1,
    is_reverse_complement: bool = False,
) -> Dict[str, Tuple[str, str, int]]:
    line_count_for_current_sequence: int = 0
    id = ""
    sequence = ""
    duplicate = None

    for line in rna_data:
        if line_count_for_current_sequence == 0:
            id = line.strip()
        elif line_count_for_current_sequence == 1:
            sequence = line.strip()

            if cutoff > 0:
                sequence = sequence[0:cutoff]

            if is_reverse_complement:
                # TODO Make sure this is correct (verify with data)
                sequence = "".join(
                    [BASE_COMPLEMENT[base] for base in reversed(sequence)]
                )

            duplicate = rna_dict.get(sequence)
        # Information from field 2 (line 3) is not needed
        # For now skip quality info, getting cutoff value supplied by user
        # elif line_count_for_current_sequence == 3:
        #     quality = line.strip()

        line_count_for_current_sequence = line_count_for_current_sequence + 1

        # Always read 4 lines per sequence, as per FASTQ format
        if line_count_for_current_sequence == 4:
            if duplicate is not None:
                rna_dict[sequence] = (
                    "".join([duplicate[0], ",", id]),
                    duplicate[1],
                    duplicate[2] + 1,
                )
                duplicate = None
            else:
                rna_dict[sequence] = (id, sequence, 1)
            line_count_for_current_sequence = 0
            id = ""
            sequence = ""

    return rna_dict


def _fill_dict_from_file(
    file_path: str,
    rna_dict: Dict[str, Tuple[str, str, int]],
    cutoff: int = -1,
    is_reverse_complement: bool = False,
) -> Dict[str, Tuple[str, str, int]]:
    with gzip.open(file_path, "rt") as rna_data_gzipped:
        try:
            rna_data_gzipped.read(1)
            print("Detected gzip file. Reading in compressed format...")
            return _convert_rna_data_to_df(
                rna_data_gzipped, rna_dict, cutoff, is_reverse_complement
            )
        except gzip.BadGzipFile:
            print("File is not a gzip file. Trying to read as uncompressed file...")

    with open(file_path, "rt") as rna_data:
        return _convert_rna_data_to_df(
            rna_data, rna_dict, cutoff, is_reverse_complement
        )


# TODO: Add parameter documentation
# cutoff is the position of the last bp in the reads
# after which the cutoff should be performed (starting with 1) (should be > 0)
def import_file(file_paths: List[str], cutoff: int = -1) -> pd.DataFrame:
    if len(file_paths) > 2:
        raise ValueError(
            (
                "Only one file (single-read sequencing) "
                "or two files (pairend-end sequencing) expected"
                "for the RNA-seq data"
            )
        )

    rna_dict: Dict[str, Tuple[str, str, int]] = {}
    # TODO: Base on file name etc instead?
    for index, file_path in enumerate(file_paths):
        # TODO: Add complement behavior for second file
        rna_dict = _fill_dict_from_file(file_path, rna_dict, cutoff, index == 1)

    rna_df = pd.DataFrame(
        list(rna_dict.values()), columns=["ids", "sequence_after_cutoff", "count"]
    )
    print(rna_df)
    return rna_df
