from typing import List, Tuple

import pandas as pd
from pepti_map.importing.rna_import import rna_importer
from pepti_map.importing.peptide_import import peptide_importer


def import_rna_and_peptides(
    rna_file_paths: List[str], peptides_file_path: str, rna_reads_cutoff: int = -1
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rna_data = rna_importer.import_file(rna_file_paths, rna_reads_cutoff)
    peptides_data = peptide_importer.import_file(peptides_file_path)
    return rna_data, peptides_data
