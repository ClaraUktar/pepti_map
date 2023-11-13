from pathlib import Path
from typing import List
from pybedtools import BedTool
from pybedtools.featurefuncs import bed2gff


def convert_bam_to_gff(path_to_bam: Path, output_path: Path):
    bam_file = BedTool(path_to_bam)
    bed_file = bam_file.bam_to_bed(cigar=True)
    gff_file = bed_file.each(bed2gff)  # pyright: ignore[reportGeneralTypeIssues]
    gff_file.moveto(output_path)


def _get_set_index_from_attributes(attributes_list: List[str]) -> int:
    attributes_dict = {
        attribute.strip().split("=")[0]: attribute.strip().split("=")[1]
        for attribute in attributes_list
    }
    return int(attributes_dict["Name"].replace('"', "").split("-")[0])


# TODO: Test?
def add_peptides_to_gff(
    path_to_gff: Path,
    peptides: List[str],
    peptide_indexes: List[List[int]],
    peptide_to_cluster_mapping: List[int],
) -> None:
    updated_lines = []
    cluster_to_peptide_mapping = {
        cluster_id: peptide_index
        for peptide_index, cluster_id in enumerate(peptide_to_cluster_mapping)
    }
    with open(path_to_gff, "rt", encoding="utf-8") as gff_file_read:
        for line in gff_file_read:
            line_parts = line.split("\t")
            entry_attributes = line_parts[8]
            attributes_list = entry_attributes.split(";")
            merged_set_index = _get_set_index_from_attributes(attributes_list)
            peptides_for_merged_set = peptide_indexes[merged_set_index]
            original_peptide_ids = [
                cluster_to_peptide_mapping[cluster_id]
                for cluster_id in peptides_for_merged_set
            ]
            original_peptide_sequences = [
                peptides[peptide_id] for peptide_id in original_peptide_ids
            ]
            attributes_list.append(f'peptides="{",".join(original_peptide_sequences)}"')
            line_parts[8] = ";".join(attributes_list)
            updated_lines.append("\t".join(line_parts) + "\n")

    with open(path_to_gff, "wt", encoding="utf-8") as gff_file_write:
        gff_file_write.writelines(updated_lines)
