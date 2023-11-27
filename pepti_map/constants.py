from enum import Enum
from pathlib import Path

PATH_TO_TEMP_FILES = Path("./temp")
PATH_TO_LAST_STEP_FILE = PATH_TO_TEMP_FILES / "last_step.txt"
PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE = (
    PATH_TO_TEMP_FILES / "peptide_to_cluster_mapping.txt"
)
PATH_TO_MATCHING_RESULT = PATH_TO_TEMP_FILES / "matching_result.txt"
PATH_TO_PRECOMPUTED_INTERSECTIONS = PATH_TO_TEMP_FILES / "precomputed_intersections.npz"
PATH_TO_MERGED_MATCHES = PATH_TO_TEMP_FILES / "merged_matches.txt"
PATH_TO_MERGED_INDEXES = PATH_TO_TEMP_FILES / "merged_indexes.txt"
PATH_TO_TRINITY_RESULTS_FILEPATHS = PATH_TO_TEMP_FILES / "trinity_results_filepaths.txt"

OUTPUT_FILENAME = "pepti_map_output.gff"


class Step(Enum):
    MATCHING = 1
    MERGING = 2
    TRINITY_INPUT = 3
    TRINITY_RUN = 4
    ALIGNMENT = 5
    PEPGENOME_INPUT = 6
    PEPGENOME_RUN = 7
