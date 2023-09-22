from pathlib import Path

PATH_TO_TEMP_FILES = Path("./temp")
PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE = (
    PATH_TO_TEMP_FILES / "peptide_to_cluster_mapping.txt"
)
PATH_TO_MATCHING_RESULT = PATH_TO_TEMP_FILES / "matching_result.txt"
PATH_TO_PRECOMPUTED_INTERSECTIONS = PATH_TO_TEMP_FILES / "precomputed_intersections.npz"
