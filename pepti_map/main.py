import logging
import click
import time
import sys
from memory_profiler import LogFile
from pepti_map.importing.peptide_import.peptide_importer import PeptideImporter
from pepti_map.importing.rna_import.rna_to_index_importer import RNAToIndexImporter


def _setup():
    logging.basicConfig(
        level=logging.DEBUG,
        filename="peptimap.log",
        filemode="a",
        format="%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    sys.stdout = LogFile("peptimap.log")


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
@click.option(
    "-k",
    "--kmer-length",
    required=False,
    type=int,
    default=7,
    show_default=True,
    help=(
        "The k-mer size used during the mapping of peptides to RNA. "
        "As the RNA is 3-frame translated for the mapping, "
        "the k-mer size refers to amino acids."
    ),
)
def main(
    peptide_file: str,
    rna_file: str,
    paired_end_file: str,
    cutoff: int,
    kmer_length: int,
):
    _setup()

    rna_files = [rna_file]
    if paired_end_file != "":
        rna_files.append(paired_end_file)

    start_time = time.time()
    kmer_index = RNAToIndexImporter(kmer_length).import_files_to_index(
        rna_files, cutoff
    )
    end_time = time.time()
    logging.info(f"Time for index creation: {end_time - start_time} seconds")
    # peptides_data = PeptideImporter().import_file(peptide_file)


if __name__ == "__main__":
    main()
