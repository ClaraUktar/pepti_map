import logging
import os
from pathlib import Path
import subprocess
from typing import List


class PoGoWrapper:
    def __init__(self):
        path_to_pogo = os.getenv("POGO_PATH")
        if path_to_pogo is None:
            raise AssertionError(
                (
                    "No path found for PoGo in the .env file. "
                    "Specify a path via POGO_PATH"
                )
            )
        self._path_to_pogo = Path(path_to_pogo)

    def run_pogo_for_directory(self, path_to_directory: Path) -> None:
        command_to_run = [
            "./PoGo",
            "-fasta",
            (path_to_directory / "pogo_fasta_in.fa").absolute().as_posix(),
            "-gtf",
            (path_to_directory / "pogo_gtf_in.gtf").absolute().as_posix(),
            "-in",
            (path_to_directory / "pogo_peptides_in.tsv").absolute().as_posix(),
            "-format",
            "bed,gtf",
            "-source",
            "pepti_map",
        ]

        subprocess.run(
            command_to_run,
            cwd=self._path_to_pogo,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    def run_pogo_for_multiple_directories(
        self, paths_to_directories: List[Path]
    ) -> None:
        logging.info("Running PoGo for all generated alignments")
        for path_to_directory in paths_to_directories:
            self.run_pogo_for_directory(path_to_directory)
