from pathlib import Path
import shutil
from typing import List

from pepti_map.constants import OUTPUT_BED_FILENAME, OUTPUT_GTF_FILENAME


class PoGoOutputHelper:
    @staticmethod
    def concat_gtf_output_into_single_file(
        paths_to_directories: List[Path], path_to_output_dir: Path
    ) -> None:
        with open(path_to_output_dir / OUTPUT_GTF_FILENAME, "wb") as output_file:
            for path_to_directory in paths_to_directories:
                with open(
                    path_to_directory / "pogo_peptides_in_out.gtf", "rb"
                ) as file_to_concat:
                    shutil.copyfileobj(file_to_concat, output_file)

    @staticmethod
    def concat_bed_output_into_single_file(
        paths_to_directories: List[Path], path_to_output_dir: Path
    ) -> None:
        with open(path_to_output_dir / OUTPUT_BED_FILENAME, "wb") as output_file:
            for path_to_directory in paths_to_directories:
                with open(
                    path_to_directory / "pogo_peptides_in.bed", "rb"
                ) as file_to_concat:
                    shutil.copyfileobj(file_to_concat, output_file)
