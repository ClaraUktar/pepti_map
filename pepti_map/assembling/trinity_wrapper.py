import logging
from pathlib import Path
import multiprocessing
import shutil
import subprocess
from typing import List

from dotenv import dotenv_values

from pepti_map.assembling.assembly_helper import AssemblyHelper
from pepti_map.constants import PATH_TO_TRINITY_RESULTS_FILEPATHS


class TrinityWrapper:
    def __init__(self, output_dir: Path, min_contig_length: int = 100):
        self._output_dir = output_dir
        self._command = []
        self._env_vars = dotenv_values()
        use_docker = self._env_vars.get("TRINITY_USE_DOCKER")
        if use_docker is not None and use_docker == "True":
            self._using_docker = True
            self._command.extend(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{output_dir.absolute().as_posix()}:/{output_dir.as_posix()}",
                    "trinityrnaseq/trinityrnaseq",
                    "Trinity",
                ]
            )
        else:
            self._using_docker = False
            trinity_path = self._env_vars.get("TRINITY_PATH")
            if trinity_path is None:
                trinity_path = "Trinity"
            self._command.append(trinity_path)

        max_mem = self._env_vars.get("TRINITY_MAX_MEM")
        if max_mem is None:
            # Enough for ~1M ~76 base Illumina paired reads
            # https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity#typical-trinity-command-line
            max_mem = "1G"
        num_cpus = self._env_vars.get("TRINITY_N_CPUS")
        if num_cpus is None:
            # Trinity default
            num_cpus = "2"
        self._num_cpus = num_cpus
        self._command.extend(
            [
                "--seqType",
                "fa",
                "--max_memory",
                max_mem,
                "--CPU",
                self._num_cpus,
                "--SS_lib_type",
                "F",
                "--full_cleanup",
                "--min_contig_length",
                str(min_contig_length),
            ]
        )

    def _get_results_from_trinity_output_file(
        self, output_dir_for_file: Path
    ) -> List[str]:
        # Delete unneeded Trinity output files
        try:
            (
                output_dir_for_file / "trinity_out_dir.Trinity.fasta.gene_trans_map"
            ).unlink()
        except FileNotFoundError:
            shutil.rmtree(output_dir_for_file / "trinity_out_dir", ignore_errors=True)

        resulting_sequences = []
        try:
            with open(
                output_dir_for_file / "trinity_out_dir.Trinity.fasta",
                "rt",
                encoding="utf-8",
            ) as result_fasta:
                current_sequence = ""
                for line in result_fasta:
                    if line.strip().startswith(">"):
                        if current_sequence != "":
                            resulting_sequences.append(current_sequence)
                            current_sequence = ""
                        continue

                    current_sequence += line.strip()

                if current_sequence != "":
                    resulting_sequences.append(current_sequence)

        except OSError:
            pass

        return resulting_sequences

    def _get_trinity_result_for_file(self, relative_filepath: Path) -> List[str]:
        command_to_run = self._command.copy()
        if self._using_docker:
            command_to_run.extend(
                [
                    "--output",
                    "/"
                    + (
                        self._output_dir / relative_filepath.parent / "trinity_out_dir"
                    ).as_posix(),
                ]
            )
            command_to_run.extend(
                ["--single", "/" + (self._output_dir / relative_filepath).as_posix()]
            )
        else:
            command_to_run.extend(
                [
                    "--output",
                    (
                        self._output_dir / relative_filepath.parent / "trinity_out_dir"
                    ).as_posix(),
                ]
            )
            command_to_run.extend(
                ["--single", (self._output_dir / relative_filepath).as_posix()]
            )

        subprocess.run(
            command_to_run, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return self._get_results_from_trinity_output_file(
            self._output_dir / relative_filepath.parent
        )

    def _write_trinity_result_to_file(self, relative_filepath: Path) -> bool:
        trinity_result = self._get_trinity_result_for_file(relative_filepath)
        if len(trinity_result) == 0:
            return False
        AssemblyHelper.write_fasta_with_sequences(
            zip(
                [
                    relative_filepath.parent.name + "-" + str(read_number)
                    for read_number in range(len(trinity_result))
                ],
                trinity_result,
            ),
            self._output_dir / relative_filepath.parent / "resulting_contigs.fa",
        )
        (
            self._output_dir
            / relative_filepath.parent
            / "trinity_out_dir.Trinity.fasta"
        ).unlink()
        return True

    def get_trinity_result_for_single_file(self, relative_filepath: Path) -> List[str]:
        trinity_results = self._get_trinity_result_for_file(relative_filepath)
        if len(trinity_results) == 0:
            logging.info(
                "No Trinity output sequences generated for set "
                + relative_filepath.parent.name
            )
        return trinity_results

    def get_trinity_result_for_multiple_files(
        self, relative_filepaths: List[Path]
    ) -> List[List[str]]:
        try:
            n_processes = self._env_vars.get("TRINITY_N_PROCESSES")
            assert isinstance(n_processes, str)
            n_processes = int(n_processes)
        except (AssertionError, ValueError):
            n_processes = multiprocessing.cpu_count()
        logging.info(f"Generating Trinity results with {n_processes} processes")
        with multiprocessing.Pool(n_processes) as pool:
            trinity_results = pool.map(
                self._get_trinity_result_for_file, relative_filepaths
            )
        for result_index, result in enumerate(trinity_results):
            if len(result) == 0:
                logging.info(
                    "No Trinity output sequences generated for set "
                    + relative_filepaths[result_index].parent.name
                )

        return trinity_results

    # TODO: Refactor?
    def write_trinity_result_for_multiple_files(
        self, relative_filepaths: List[Path]
    ) -> List[Path]:
        try:
            n_processes = self._env_vars.get("TRINITY_N_PROCESSES")
            assert isinstance(n_processes, str)
            n_processes = int(n_processes)
        except (AssertionError, ValueError):
            n_processes = multiprocessing.cpu_count()
            n_processes = n_processes // int(self._num_cpus)
        logging.info(f"Generating Trinity results with {n_processes} processes")
        with multiprocessing.Pool(n_processes) as pool:
            trinity_results = pool.map(
                self._write_trinity_result_to_file, relative_filepaths
            )
        for result_index, result in enumerate(trinity_results):
            if not result:
                logging.info(
                    "No Trinity output sequences generated for set "
                    + relative_filepaths[result_index].parent.name
                )
        return [
            (self._output_dir / relative_filepath.parent / "resulting_contigs.fa")
            for relative_filepath, has_result in zip(
                relative_filepaths, trinity_results
            )
            if has_result
        ]

    # TODO: Add test
    @staticmethod
    def save_results_filepaths(
        results_filepaths: List[Path],
        path_to_results_filepaths: Path = PATH_TO_TRINITY_RESULTS_FILEPATHS,
    ) -> None:
        with open(
            path_to_results_filepaths, "wt", encoding="utf-8"
        ) as results_filepaths_file:
            results_filepaths_file.writelines(
                [
                    results_filepath.as_posix() + "\n"
                    for results_filepath in results_filepaths
                ]
            )

    # TODO: Add test
    @staticmethod
    def load_results_filepaths(
        path_to_results_filepaths: Path = PATH_TO_TRINITY_RESULTS_FILEPATHS,
    ) -> List[Path]:
        results_filepaths: List[Path] = []
        with open(
            path_to_results_filepaths, "rt", encoding="utf-8"
        ) as results_filepaths_file:
            for line in results_filepaths_file:
                line = line.strip()
                results_filepaths.append(Path(line))
        return results_filepaths
