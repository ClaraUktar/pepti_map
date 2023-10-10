import os
from pathlib import Path
import subprocess
from typing import List


class TrinityWrapper:
    def __init__(self, output_dir: Path):
        self._output_dir = output_dir
        self._command = []
        use_docker = os.getenv("USE_DOCKER")
        if use_docker is not None and bool(use_docker):
            self._using_docker = True
            self._command.extend(
                [
                    "docker",
                    "run",
                    "--rm",
                    "-v",
                    f"{output_dir.absolute().as_posix()}:/temp",
                    "trinityrnaseq/trinityrnaseq",
                    "Trinity",
                ]
            )
            self._command.extend(["--output", "/temp/trinity_out_dir"])
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
                ["--output", (output_dir / "trinity_out_dir").as_posix()]
            )

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
                "100",
                "--no_normalize_reads",
                "--single",
            ]
        )

    def _get_results_from_trinity_output_file(self) -> List[str]:
        return []

    def get_trinity_result_for_file(self, filename: str) -> List[str]:
        command_to_run = self._command.copy()
        if self._using_docker:
            command_to_run.append(f"/temp/{filename}")
        else:
            command_to_run.append((self._output_dir / filename).as_posix())
        subprocess.run(
            command_to_run, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return self._get_results_from_trinity_output_file()
