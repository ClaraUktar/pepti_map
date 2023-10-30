import os
import logging
from pathlib import Path
import shutil
import subprocess
from typing import List


class TrinityWrapper:
    def __init__(self, output_dir: Path, min_contig_length: int = 100):
        self._output_dir = output_dir
        self._command = []
        use_docker = os.getenv("USE_DOCKER")
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
            trinity_path = os.getenv("TRINITY_PATH")
            if trinity_path is None:
                raise AssertionError(
                    (
                        "No path found for Trinity in the .env file. Specify a path "
                        "via TRINITY_PATH or use Docker (USE_DOCKER=True)"
                    )
                )
            self._command.append(trinity_path)

        self._command.extend(
            [
                "--seqType",
                "fa",
                "--max_memory",
                "1G",
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
            logging.info(
                "No Trinity output sequences generated for set "
                + output_dir_for_file.name
            )

        return resulting_sequences

    def get_trinity_result_for_file(self, relative_filepath: Path) -> List[str]:
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
