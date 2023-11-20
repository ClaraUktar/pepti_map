import logging
import multiprocessing
import os
from pathlib import Path
from typing import List
import subprocess


class GmapWrapper:
    def __init__(self):
        try:
            n_threads = os.getenv("GMAP_N_THREADS")
            assert isinstance(n_threads, str)
            n_threads = int(n_threads)
        except (AssertionError, ValueError):
            n_threads = multiprocessing.cpu_count()

        self._n_threads = n_threads

    def build_index(
        self,
        files_to_index: List[Path],
        index_basename: str = "gmap_index",
        path_to_index_directory: Path = Path("./"),
    ) -> None:
        self._index_basename = index_basename
        self._path_to_index_directory = path_to_index_directory
        logging.info(
            (
                f"Generating GMAP index with name {self._index_basename} "
                f"in directory {self._path_to_index_directory.absolute().as_posix()}"
            )
        )

        files_are_gzipped = files_to_index[0].name.endswith(".gz")
        build_index_command = ["gmap_build"]
        if files_are_gzipped:
            build_index_command.append("-g")
        build_index_command.extend(
            [
                "-D",
                path_to_index_directory.absolute().as_posix(),
                "-d",
                index_basename,
            ]
        )
        build_index_command.extend(
            [file_to_index.absolute().as_posix() for file_to_index in files_to_index]
        )

        subprocess.run(
            build_index_command,
            stdout=subprocess.DEVNULL,
        )

    def use_existing_index(
        self, index_basename: str, path_to_index_directory: Path
    ) -> None:
        self._index_basename = index_basename
        self._path_to_index_directory = path_to_index_directory
        logging.info(
            (
                f"Using existing GMAP index with name {self._index_basename}, "
                f"located in directory {self._path_to_index_directory.absolute().as_posix()}"
            )
        )

    def produce_alignment(
        self, paths_to_sequences: List[Path], output_filepath: Path
    ) -> None:
        if not self._index_basename or not self._path_to_index_directory:
            error_message = (
                "Cannot generate an alignment without an index generated or set."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        logging.info(f"Generating GMAP alignment with {self._n_threads} threads")

        alignment_command = [
            "gmap",
            "-D",
            self._path_to_index_directory.absolute().as_posix(),
            "-d",
            self._index_basename,
            "-t",
            str(self._n_threads),
            "-f",
            "gff3_gene",
        ]
        alignment_command.extend(
            [
                path_to_sequences.absolute().as_posix()
                for path_to_sequences in paths_to_sequences
            ]
        )
        alignment_command.extend([">", output_filepath.absolute().as_posix()])

        subprocess.run(alignment_command)
