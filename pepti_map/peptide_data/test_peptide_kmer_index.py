from collections import defaultdict
from pathlib import Path
from pepti_map.importing.peptide_import.testdata_peptide_importer import (
    EXPECTED_RESULT_INDEX,
)
from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex


class TestPeptideKmerIndex:
    def test_dump_and_load_index_file(self, tmp_path):
        kmer_index = PeptideKmerIndex()
        assert kmer_index.kmer_index == {}

        kmer_index.kmer_index = defaultdict(list, EXPECTED_RESULT_INDEX.copy())
        kmer_index.number_of_peptides = 10
        temp_directory: Path = tmp_path / "kmer_index"
        temp_directory.mkdir()
        filepath = temp_directory.as_posix() + "/index.txt"
        kmer_index.dump_index_to_file(filepath)

        new_kmer_index = PeptideKmerIndex.load_index_from_file(filepath)
        assert new_kmer_index.kmer_index == EXPECTED_RESULT_INDEX
        assert new_kmer_index.number_of_peptides == 10
