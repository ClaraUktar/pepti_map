import logging
import os
from pathlib import Path
import subprocess
from typing import List


class PepGenomeWrapper:
    def __init__(self):
        path_to_pepgenome = os.getenv("PEPGENOME_PATH")
        if path_to_pepgenome is None:
            raise AssertionError(
                (
                    "No path found for PepGenome in the .env file. "
                    "Specify a path via PEPGENOME_PATH"
                )
            )
        self._path_to_pepgenome = Path(path_to_pepgenome)

    def run_pepgenome_for_directory(self, path_to_directory: Path) -> None:
        pepgenome_command = [
            "java",
            "-cp",
            "'./*'",
            "org.bigbio.pgatk.pepgenome.PepGenomeTool",
            "-fasta",
            (path_to_directory / "pepgenome_fasta_in.fa").absolute().as_posix(),
            "-gff",
            (path_to_directory / "pepgenome_gff_in.gff3").absolute().as_posix(),
            "-in",
            (path_to_directory / "pepgenome_peptides_in.tsv").absolute().as_posix(),
            "-format",
            "bed,gtf",
            "-exco",
            "true",
            "-source",
            "pepti_map",
        ]
        path_to_pepgenome_script = path_to_directory / "run-pepgenome.sh"
        with open(path_to_pepgenome_script, "wt", encoding="utf-8") as pepgenome_script:
            pepgenome_script.write(" ".join(pepgenome_command))

        command_to_run = ["sh", path_to_pepgenome_script.absolute().as_posix()]
        subprocess.run(
            command_to_run,
            cwd=self._path_to_pepgenome,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

    def run_pepgenome_for_multiple_directories(
        self, paths_to_directories: List[Path]
    ) -> None:
        logging.info("Running PepGenome for all generated alignments")
        for path_to_directory in paths_to_directories:
            self.run_pepgenome_for_directory(path_to_directory)
