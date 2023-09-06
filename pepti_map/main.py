import logging
from pathlib import Path

# import time
import click
from pepti_map.importing.peptide_import.peptide_to_index_importer import (
    PeptideToIndexImporter,
)
from pepti_map.importing.rna_import.lazy_rna_reader import LazyRNAReader
from pepti_map.matching.rna_to_peptide_matcher import RNAToPeptideMatcher


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
# @profile
def main(
    peptide_file: str,
    rna_file: str,
    paired_end_file: str,
    cutoff: int,
    kmer_length: int,
):
    _setup()
    # TODO: Delete all classes not needed here!

    # TODO: Also need to read in the RNA file to access original sequences via id

    # TODO: Add progress indications

    # TODO: Recreate clean version of requirements.txt

    # TODO: Add full docstrings for all relevant methods

    # TODO: Use numpy arrays everywhere?

    kmer_index = PeptideToIndexImporter(kmer_length).import_file_to_index(peptide_file)

    rna_files = [Path(rna_file)]
    if paired_end_file != "":
        rna_files.append(Path(paired_end_file))

    # start_matching = time.time()
    number_peptides = kmer_index.number_of_peptides
    matcher = RNAToPeptideMatcher(kmer_index, kmer_index.number_of_peptides)
    for sequence_id, sequence in LazyRNAReader(rna_files, cutoff):
        matcher.add_peptide_matches_for_rna_read(sequence_id, sequence)
    # end_matching = time.time()
    # logging.info(f"Time for matching: {end_matching - start_matching} sec")

    del kmer_index
    # print(matcher.matches)

    number_none_entries = 0
    sum_set_elements = 0
    for match in matcher.matches:
        if match is None:
            number_none_entries += 1
            continue
        sum_set_elements += len(match)
    logging.info(f"Total number of peptides: {number_peptides}")
    logging.info(f'Number of "None" entries: {number_none_entries}')
    logging.info(
        f'Number of non-"None" entries: {number_peptides - number_none_entries}'
    )
    logging.info(
        (
            'Average number of elements in non-"None" entries: '
            f"{sum_set_elements / (number_peptides - number_none_entries)}"
        )
    )


if __name__ == "__main__":
    main()
