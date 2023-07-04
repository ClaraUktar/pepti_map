from typing import List
from pepti_map.importing import peptide_importer, rna_importer


def import_rna_and_peptides(
    rna_file_paths: List[str], peptides_file_path: str, rna_reads_cutoff: int = -1
):
    rna_data = rna_importer.import_file(rna_file_paths, rna_reads_cutoff)
    peptides_data = peptide_importer.import_file(peptides_file_path)
