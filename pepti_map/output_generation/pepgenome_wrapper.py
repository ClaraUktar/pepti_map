import logging
from pathlib import Path
import subprocess
from typing import List


class PepGenomeWrapper:
    def __init__(self, path_to_pepgenome: Path):
        self._path_to_pepgenome = path_to_pepgenome

    def run_pepgenome_for_directory(self, path_to_directory: Path) -> None:
        # TODO: Do we need to allow mismatches? (-mm)
        command_to_run = [
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
            "-source",
            "pepti_map",
        ]

        subprocess.run(command_to_run, cwd=self._path_to_pepgenome)

    def run_pepgenome_for_multiple_directories(
        self, paths_to_directories: List[Path]
    ) -> None:
        logging.info("Running PepGenome for all generated alignments")
        for path_to_directory in paths_to_directories:
            self.run_pepgenome_for_directory(path_to_directory)
