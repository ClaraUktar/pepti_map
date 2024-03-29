import logging
import multiprocessing
from pathlib import Path
from typing import List
import subprocess

from dotenv import dotenv_values


class GmapWrapper:
    def __init__(self):
        env_vars = dotenv_values()
        try:
            n_threads = env_vars.get("GMAP_N_THREADS")
            assert isinstance(n_threads, str)
            n_threads = int(n_threads)
        except (AssertionError, ValueError):
            n_threads = multiprocessing.cpu_count()

        try:
            batch_mode = env_vars.get("GMAP_BATCH_MODE")
            assert isinstance(batch_mode, str)
            batch_mode = int(batch_mode)
        except (AssertionError, ValueError):
            batch_mode = 2

        self._n_threads = n_threads
        self._batch_mode = batch_mode

    def build_index(
        self,
        files_to_index: List[Path],
        path_to_index_directory: Path = Path("./"),
        index_basename: str = "gmap_index",
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
                f"Using existing GMAP index with name {self._index_basename}, located "
                f"in directory {self._path_to_index_directory.absolute().as_posix()}"
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

        files_list_str = ",".join(
            [
                path_to_sequences.absolute().as_posix()
                for path_to_sequences in paths_to_sequences
            ]
        )
        logging.info(
            (
                f"Generating GMAP alignment with {self._n_threads} threads, "
                f"batch mode {str(self._batch_mode)} "
                f"for file(s): {files_list_str}"
            )
        )

        alignment_command = [
            "gmap",
            "-D",
            self._path_to_index_directory.absolute().as_posix(),
            "-d",
            self._index_basename,
            "-B",
            str(self._batch_mode),
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

        with open(output_filepath, "wt", encoding="utf-8") as output_file:
            subprocess.run(
                alignment_command, stderr=subprocess.DEVNULL, stdout=output_file
            )

    def _produce_alignment_for_single_file(self, path_to_sequences: Path) -> None:
        alignment_command = [
            "gmap",
            "-D",
            self._path_to_index_directory.absolute().as_posix(),
            "-d",
            self._index_basename,
            "-B",
            str(self._batch_mode),
            "-t",
            str(self._n_threads),
            "-f",
            "gff3_gene",
            path_to_sequences.absolute().as_posix(),
        ]
        with open(
            path_to_sequences.parent / "alignment_result.gff3", "wt", encoding="utf-8"
        ) as output_file:
            subprocess.run(
                alignment_command, stderr=subprocess.DEVNULL, stdout=output_file
            )

    def produce_alignment_for_multiple_files(
        self, paths_to_sequences: List[Path]
    ) -> None:
        if not self._index_basename or not self._path_to_index_directory:
            error_message = (
                "Cannot generate an alignment without an index generated or set."
            )
            logging.error(error_message)
            raise ValueError(error_message)

        logging.info(
            (
                f"Generating GMAP alignments with {self._n_threads} threads, "
                f"batch mode {str(self._batch_mode)}"
            )
        )
        for path_to_sequences in paths_to_sequences:
            self._produce_alignment_for_single_file(path_to_sequences)
