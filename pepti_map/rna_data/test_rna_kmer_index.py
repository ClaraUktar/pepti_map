from collections import defaultdict
from pathlib import Path
from pepti_map.importing.rna_import.testdata_rna_importer import (
    EXPECTED_RESULT_INDEX_SINGLE_END,
)
from pepti_map.rna_data.rna_kmer_index import RNAKmerIndex

# TODO: Delete


class TestRNAKmerIndex:
    def test_dump_and_load_index_file(self, tmp_path):
        kmer_index = RNAKmerIndex()
        assert kmer_index.kmer_index == {}

        kmer_index.kmer_index = defaultdict(
            list, EXPECTED_RESULT_INDEX_SINGLE_END.copy()
        )
        temp_directory: Path = tmp_path / "kmer_index"
        temp_directory.mkdir()
        filepath = temp_directory.as_posix() + "/index.txt"
        kmer_index.dump_index_to_file(filepath)

        new_kmer_index = RNAKmerIndex.load_index_from_file(filepath)
        assert new_kmer_index.kmer_index == EXPECTED_RESULT_INDEX_SINGLE_END
