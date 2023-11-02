import logging
import multiprocessing
import os
from pathlib import Path
from typing import List
import subprocess


class BowtieWrapper:
    def __init__(self):
        try:
            n_threads = os.getenv("BOWTIE_N_THREADS")
            assert isinstance(n_threads, str)
            n_threads = int(n_threads)
        except (AssertionError, ValueError):
            n_threads = multiprocessing.cpu_count()

        self._n_threads = n_threads

    def build_index(
        self, files_to_index: List[Path], index_basename: str = "bowtie_index"
    ) -> None:
        self._index_basename = index_basename
        logging.info(
            (
                f"Generating Bowtie index (basename: {self._index_basename}) "
                f"with {self._n_threads} threads"
            )
        )
        # TODO: Do we want the --noref option?
        build_index_command = [
            "bowtie2-build",
            "--threads",
            str(self._n_threads),
            "--noref",  # do not build the index portions for paired-end alignment
            ",".join(
                [
                    file_to_index.absolute().as_posix()
                    for file_to_index in files_to_index
                ]
            ),
            index_basename,
        ]
        subprocess.run(build_index_command)

    def use_existing_index(self, index_basename: str) -> None:
        self._index_basename = index_basename
        logging.info(
            f"Using existing Bowtie index with basename: {self._index_basename}"
        )

    def produce_alignment(
        self, paths_to_sequences: List[Path], output_filepath: Path
    ) -> None:
        if not self._index_basename:
            logging.error(
                "Cannot generate an alignment without an index generated or set."
            )
            raise ValueError(
                "Cannot generate an alignment without an index generated or set."
            )

        logging.info(f"Generating Bowtie alignment with {self._n_threads} threads")

        alignment_command = [
            "bowtie2",
            "-f",  # reads are given in FASTA format
            # TODO: Should we use "--reorder" option for same order of output as input?
            "--threads",
            str(self._n_threads),
            "-x",
            self._index_basename,
            "-U",  # reads are unpaired
            ",".join(
                [
                    path_to_sequences.absolute().as_posix()
                    for path_to_sequences in paths_to_sequences
                ]
            ),
            "-S",
            output_filepath.absolute().as_posix(),
        ]
        subprocess.run(alignment_command)
