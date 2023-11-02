from pepti_map.assembling.assembly_helper import AssemblyHelper


class TestAssemblyInputGenerator:
    def test_write_single_sequence(self, tmp_path):
        tmp_file = tmp_path / "sequence_fasta.fa"
        read_id = "11"
        read_sequence = (
            "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGT"
            "GTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
        )
        AssemblyHelper.write_fasta_with_sequences(
            [
                (
                    read_id,
                    read_sequence,
                ),
            ],
            tmp_file,
        )
        expected_file_content = [f">{read_id}", read_sequence]
        with open(tmp_file, "rt", encoding="utf-8") as fasta_file:
            assert [
                line.strip() for line in fasta_file.readlines()
            ] == expected_file_content

    def test_write_multiple_sequences(self, tmp_path):
        tmp_file = tmp_path / "sequence_fasta.fa"
        read_id_1 = "11"
        read_sequence_1 = (
            "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGT"
            "GTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
        )
        read_id_2 = "302"
        read_sequence_2 = (
            "CTTAAGTCATCGTGAAGGGAGGCATACGTGTTGGACAGAGCCCGAGACGC"
            "GGAATTTGCAAGATTTTATGGAAGAGGGGCGGCTCATTGCAAATCGAGCTC"
        )
        AssemblyHelper.write_fasta_with_sequences(
            [
                (
                    read_id_1,
                    read_sequence_1,
                ),
                (read_id_2, read_sequence_2),
            ],
            tmp_file,
        )
        expected_file_content = [
            f">{read_id_1}",
            read_sequence_1,
            f">{read_id_2}",
            read_sequence_2,
        ]
        with open(tmp_file, "rt", encoding="utf-8") as fasta_file:
            assert [
                line.strip() for line in fasta_file.readlines()
            ] == expected_file_content
