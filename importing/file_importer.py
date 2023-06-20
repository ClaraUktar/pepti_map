import rna_importer, peptide_importer


def import_rna_and_peptides(rna_file_path: str, peptides_file_path: str):
    rna_data = rna_importer.import_file(rna_file_path)
    peptides_data = peptide_importer.import_file(peptides_file_path)
