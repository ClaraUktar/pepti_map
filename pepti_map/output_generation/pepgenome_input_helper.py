from collections import defaultdict
from pathlib import Path
from typing import List


class PepGenomeInputHelper:
    def __init__(
        self, path_to_peptides: Path, path_to_peptide_to_cluster_mapping: Path
    ):
        self._peptides = []
        with open(path_to_peptides, "rt", encoding="utf-8") as peptides_file:
            for line in peptides_file:
                self._peptides.append(
                    line.split("\t")[0].strip()
                )  # ignores group information if present

        # Reverse mapping to get cluster_id -> peptide_id
        self._cluster_to_peptide_mapping: "defaultdict[int, List[int]]" = defaultdict(
            list
        )
        with open(
            path_to_peptide_to_cluster_mapping, "rt", encoding="utf-8"
        ) as peptide_to_cluster_mapping:
            for peptide_id, line in enumerate(peptide_to_cluster_mapping):
                cluster_id = int(line.strip())
                # Exclude peptides that were too short to be included
                if cluster_id == -1:
                    continue
                self._cluster_to_peptide_mapping[cluster_id].append(peptide_id)

    def generate_peptide_input_file(
        self, output_path: Path, merged_indexes: List[int]
    ) -> None:
        peptide_already_written: List[bool] = [
            False for _ in range(len(self._peptides))
        ]
        with open(output_path, "wt", encoding="utf-8") as new_input_file:
            for cluster_id in merged_indexes:
                peptide_ids = self._cluster_to_peptide_mapping[cluster_id]
                for peptide_id in peptide_ids:
                    if peptide_already_written[peptide_id]:
                        continue
                    # Use 1 as default for Sample, PSMs and Quant
                    new_input_file.write(
                        "\t".join(["1", self._peptides[peptide_id], "1", "1"]) + "\n"
                    )
                    peptide_already_written[peptide_id] = True

    def generate_all_peptide_input_files(
        self, output_paths: List[Path], path_to_merged_indexes: Path
    ) -> None:
        # TODO: Refactor to use code from match merger
        merged_indexes: List[List[int]] = []
        with open(
            path_to_merged_indexes, "rt", encoding="utf-8"
        ) as peptide_indexes_file:
            for line in peptide_indexes_file:
                line = line.strip()
                merged_indexes.append(
                    [int(match_elem) for match_elem in line.split(",")]
                )

        for set_index, output_path in enumerate(output_paths):
            self.generate_peptide_input_file(output_path, merged_indexes[set_index])


def generate_gff_input_file(path_to_gff: Path, output_path: Path) -> None:
    # TODO: For each entry in GFF, adapt positions
    pass


def generate_protein_fasta_input_file(
    path_to_contig_sequences: Path, output_path: Path
) -> None:
    # TODO: Create 3-frame translation sequences for each contig
    # Create correct "header" for each sequence in FASTA
    pass
