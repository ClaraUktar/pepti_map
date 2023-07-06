import logging
import click
import pepti_map.importing.file_importer as file_importer


def _setup():
    logging.basicConfig(level=logging.DEBUG)


@click.command()
@click.option(
    "-p",
    "--peptide-file",
    required=True,
    type=str,
    help="The path to the peptide file.",
)
@click.option(
    "-r", "--rna-file", required=True, type=str, help="The path to the RNA file."
)
@click.option(
    "-pa",
    "--paired-end-file",
    required=False,
    type=str,
    default="",
    show_default=True,
    help=(
        "The path to the second RNA file in case of paired-end sequencing. "
        'If none is given, the RNA file given with the "-r" option '
        "is assumed to result from single-end sequencing."
    ),
)
@click.option(
    "-c",
    "--cutoff",
    required=False,
    type=int,
    default=-1,
    show_default=True,
    help=(
        "The position of the last base in the reads "
        "after which a cutoff should be performed (starting at 1). "
        "The cutoff is applied to all reads."
        "If the value is equal to or smaller than 0, no cutoff is performed."
    ),
)
def main(
    peptide_file: str,
    rna_file: str,
    paired_end_file: str,
    cutoff: int,
):
    _setup()

    rna_files = [rna_file]
    if paired_end_file != "":
        rna_files.append(paired_end_file)

    file_importer.import_rna_and_peptides(
        rna_files,
        peptide_file,
        cutoff,
    )


if __name__ == "__main__":
    main()
