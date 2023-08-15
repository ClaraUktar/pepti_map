import logging
import click
from os.path import isfile
from pepti_map.importing.peptide_import.peptide_importer import PeptideImporter
from pepti_map.importing.peptide_import.peptide_to_index_importer import (
    PeptideToIndexImporter,
)
from pepti_map.importing.rna_import.rna_to_index_importer import RNAToIndexImporter
from pepti_map.matching.peptide_to_rna_matcher import PeptideToRNAMatcher
from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.rna_data.rna_kmer_index import RNAKmerIndex


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
    show_default=False,
    help=(
        "The path to the second RNA file in case of paired-end sequencing. "
        'If none is given, the RNA file given with the "-r" option '
        "is assumed to result from single-end sequencing."
    ),
)
@click.option(
    "-i",
    "--index-file",
    required=False,
    type=str,
    default="",
    show_default=False,
    help=(
        "The path to an index file for an already constructed "
        "k-mer index based on RNA reads. Alternatively, if the file "
        "for the given path does not yet exist, the k-mer index is constructed "
        "based on the given RNA file(s) and saved under this path."
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
    index_file: str,
    cutoff: int,
    kmer_length: int,
):
    _setup()

    # TODO: Also need to read in the RNA file to access original sequences via id
    kmer_index: PeptideKmerIndex
    # TODO: Improve checking whether index file already exists, or use two options
    # Only necessary if peptide file is big
    if index_file != "" and isfile(index_file):
        kmer_index = PeptideKmerIndex.load_index_from_file(index_file)
    else:
        kmer_index = PeptideToIndexImporter(kmer_length).import_file_to_index(
            peptide_file
        )

        if index_file != "":
            kmer_index.dump_index_to_file(index_file)

    # peptides_data = PeptideImporter().import_file(peptide_file)

    # rna_files = [rna_file]
    # if paired_end_file != "":
    #     rna_files.append(paired_end_file)

    # matched_peptides = PeptideToRNAMatcher(
    #     peptides_data
    # ).find_rna_read_matches_for_peptides(kmer_index)
    del kmer_index


if __name__ == "__main__":
    main()
