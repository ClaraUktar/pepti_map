from pathlib import Path
from typing import List


def generate_peptide_input_file(
    path_to_peptides: Path, output_path: Path, peptide_to_cluster_mapping: List[int]
) -> None:
    # TODO: Only use peptides for set
    peptides = []
    with open(path_to_peptides, "rt", encoding="utf-8") as peptides_file:
        for line in peptides_file:
            peptides.append(
                line.split("\t")[0].strip()
            )  # ignores group information if present

    with open(output_path, "wt", encoding="utf-8") as new_input_file:
        for peptide_index, peptide in enumerate(peptides):
            if (
                peptide_to_cluster_mapping[peptide_index] == -1
            ):  # Exclude peptides that were not part of the pipeline
                continue
            # Use 1 as default for Sample, PSMs and Quant
            new_input_file.write("\t".join(["1", peptide, "1", "1"]) + "\n")


def generate_gff_input_file(path_to_gff: Path, output_path: Path) -> None:
    # TODO: For each entry in GFF, adapt positions
    pass


def generate_protein_fasta_input_file(
    path_to_contig_sequences: Path, output_path: Path
) -> None:
    # TODO: Create 3-frame translation sequences for each contig
    # Create correct "header" for each sequence in FASTA
    pass
