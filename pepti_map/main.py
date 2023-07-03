import importing.file_importer as file_importer


def main():
    file_importer.import_rna_and_peptides(
        "../../../Daten/rna_example.fastq",
        "../../../Daten/ccle_ks_hct116_peptide_list.txt",
        -1,
    )


if __name__ == "__main__":
    main()
