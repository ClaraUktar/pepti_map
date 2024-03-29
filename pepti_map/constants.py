from enum import Enum
from pathlib import Path
from dotenv import dotenv_values

temp_dir_path = dotenv_values().get("TEMP_DIR_PATH")
if temp_dir_path is None:
    temp_dir_path = "./temp"

PATH_TO_TEMP_FILES = Path(temp_dir_path)
PATH_TO_LAST_STEP_FILE = PATH_TO_TEMP_FILES / "last_step.txt"
PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE = (
    PATH_TO_TEMP_FILES / "peptide_to_cluster_mapping.txt"
)
PATH_TO_MATCHING_RESULT = PATH_TO_TEMP_FILES / "matching_result.txt"
PATH_TO_PRECOMPUTED_INTERSECTIONS = PATH_TO_TEMP_FILES / "precomputed_intersections.npz"
PATH_TO_MERGED_MATCHES = PATH_TO_TEMP_FILES / "merged_matches.txt"
PATH_TO_MERGED_INDEXES = PATH_TO_TEMP_FILES / "merged_indexes.txt"
PATH_TO_TRINITY_RESULTS_FILEPATHS = PATH_TO_TEMP_FILES / "trinity_results_filepaths.txt"

PEPTIDE_READ_QUANT_FILENAME = "peptide_read_quant.tsv"
OUTPUT_GTF_FILENAME = "pepti_map_output.gtf"
OUTPUT_BED_FILENAME = "pepti_map_output.bed"


class Step(Enum):
    MATCHING = 1
    MERGING = 2
    TRINITY_INPUT = 3
    TRINITY_RUN = 4
    ALIGNMENT = 5
    POGO_INPUT = 6
    POGO_RUN = 7
    POGO_OUTPUT = 8
