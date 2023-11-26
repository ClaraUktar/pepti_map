import logging
import multiprocessing
import os
from pathlib import Path
import subprocess
from typing import List


class PepGenomeWrapper:
    def __init__(self, path_to_pepgenome: Path):
        self._path_to_pepgenome = path_to_pepgenome
        try:
            n_threads = os.getenv("PEPGENOME_N_THREADS")
            assert isinstance(n_threads, str)
            n_threads = int(n_threads)
        except (AssertionError, ValueError):
            n_threads = multiprocessing.cpu_count()

        self._n_threads = n_threads

    def run_pepgenome_for_directory(self, path_to_directory: Path) -> None:
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
        with multiprocessing.Pool(self._n_threads) as pool:
            pool.map(self.run_pepgenome_for_directory, paths_to_directories)
