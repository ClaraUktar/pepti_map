import logging
from dotenv import load_dotenv
from pathlib import Path
import shutil
from typing import List, Literal, Set, Tuple, Union
import click
import numpy as np
import numpy.typing as npt
from pepti_map.aligning.bowtie_wrapper import BowtieWrapper
from pepti_map.assembling.assembly_helper import AssemblyHelper
from pepti_map.assembling.trinity_wrapper import TrinityWrapper
from pepti_map.constants import (
    OUTPUT_FILENAME,
    PATH_TO_LAST_STEP_FILE,
    PATH_TO_TEMP_FILES,
    Step,
)
from pepti_map.importing.peptide_import.peptide_importer import PeptideImporter
from pepti_map.importing.peptide_import.peptide_to_index_importer import (
    PeptideToIndexImporter,
)
from pepti_map.importing.rna_import.lazy_rna_reader import LazyRNAReader
from pepti_map.importing.rna_import.rna_reads_retriever import RNAReadsRetriever
from pepti_map.matching.match_merger import MatchMerger
from pepti_map.matching.precomputing_rna_to_peptide_matcher import (
    PrecomputingRNAToPeptideMatcher,
)
from pepti_map.matching.rna_to_peptide_matcher import RNAToPeptideMatcher
from pepti_map.util.file_conversion import convert_bam_to_gff


def _setup():
    logging.basicConfig(level=logging.DEBUG)
    PATH_TO_TEMP_FILES.mkdir(exist_ok=True)
    load_dotenv()


def _teardown():
    shutil.rmtree(PATH_TO_TEMP_FILES)


def _write_last_step(step: int):
    with open(PATH_TO_LAST_STEP_FILE, "wt", encoding="utf-8") as last_step_file:
        last_step_file.write(str(step))


def _get_last_step() -> int:
    try:
        with open(PATH_TO_LAST_STEP_FILE, "rt", encoding="utf-8") as last_step_file:
            return int(last_step_file.readline().strip())
    except FileNotFoundError:
        return -1


def load_matches(
    precompute_intersections: bool,
) -> Tuple[List[Union[Set[int], None]], Union[npt.NDArray[np.uint32], None]]:
    matches = RNAToPeptideMatcher.load_matches()
    logging.info("Loaded matches.")
    if precompute_intersections:
        precomputed_intersections = (
            PrecomputingRNAToPeptideMatcher.load_precomputed_intersections()
        )
        logging.info("Loaded precomputed intersections.")
    else:
        precomputed_intersections = None
    return matches, precomputed_intersections


def compute_matches(
    peptide_file: str,
    rna_file: str,
    paired_end_file: str,
    cutoff: int,
    kmer_length: int,
    output_dir: str,
    precompute_intersections: bool,
) -> Tuple[List[Union[Set[int], None]], Union[npt.NDArray[np.uint32], None]]:
    kmer_index, peptide_to_cluster_mapping = PeptideToIndexImporter(
        kmer_length
    ).import_file_to_index(Path(peptide_file))
    logging.info("Imported peptides to index.")

    rna_files = [Path(rna_file)]
    if paired_end_file != "":
        rna_files.append(Path(paired_end_file))

    if precompute_intersections:
        matcher = PrecomputingRNAToPeptideMatcher(
            kmer_index, kmer_index.number_of_peptides, peptide_to_cluster_mapping
        )
        logging.info("Precomputing intersections during matching.")
    else:
        matcher = RNAToPeptideMatcher(
            kmer_index, kmer_index.number_of_peptides, peptide_to_cluster_mapping
        )
    logging.info("Matching RNA-seq reads to peptides...")
    for sequence_id, sequence in LazyRNAReader(rna_files, cutoff):
        matcher.add_peptide_matches_for_rna_read(sequence_id, sequence)
    logging.info("Generated all matches.")
    del kmer_index
    matcher.write_peptide_read_quant_file(
        Path(output_dir), PeptideImporter().import_file(Path(peptide_file))
    )
    matcher.save_matches()
    matches = matcher.get_matches()
    if isinstance(matcher, PrecomputingRNAToPeptideMatcher):
        matcher.save_precomputed_intersections()
        precomputed_intersections = matcher.get_precomputed_intersections()
    else:
        precomputed_intersections = None
    _write_last_step(Step.MATCHING.value)
    logging.info("Saved results of matching step.")
    return matches, precomputed_intersections


def load_merge_results() -> Tuple[List[Set[int]], List[List[int]]]:
    merged_sets, peptide_indexes = MatchMerger.load_merged_result()
    logging.info("Loaded merge results.")
    return merged_sets, peptide_indexes


def compute_merge_results(
    jaccard_index_threshold: float,
    merging_method: Literal["agglomerative-clustering", "full-matrix"],
    matches: List[Union[Set[int], None]],
    precomputed_intersections: Union[npt.NDArray[np.uint32], None],
) -> Tuple[List[Set[int]], List[List[int]]]:
    logging.info(f"Using merging method: {merging_method}.")
    logging.info(
        f"Only merging sets with a Jaccard Index value above {jaccard_index_threshold}."
    )
    merged_sets, peptide_indexes = MatchMerger(
        matches, jaccard_index_threshold, precomputed_intersections
    ).merge_matches(merging_method)
    logging.info("Completed merging.")
    MatchMerger.save_merged_result(merged_sets, peptide_indexes)
    _write_last_step(Step.MERGING.value)
    logging.info("Saved results of merging step.")
    return merged_sets, peptide_indexes


def generate_trinity_input(
    rna_file: str, paired_end_file: str, cutoff: int, merged_sets: List[Set[int]]
) -> List[Path]:
    rna_files = [Path(rna_file)]
    if paired_end_file != "":
        rna_files.append(Path(paired_end_file))
    rna_reads_retriever = RNAReadsRetriever(rna_files, cutoff)

    relative_filepaths: List[Path] = []
    # TODO: Could be parallelized?
    for set_index, merged_set in enumerate(merged_sets):
        read_ids, read_sequences = rna_reads_retriever.get_read_sequences_for_ids(
            list(merged_set)
        )
        (PATH_TO_TEMP_FILES / f"{str(set_index)}").mkdir()
        relative_filepath = Path(f"{str(set_index)}/{str(set_index)}.fa")
        AssemblyHelper.write_fasta_with_sequences(
            (zip([str(read_id) for read_id in read_ids], read_sequences)),
            PATH_TO_TEMP_FILES / relative_filepath,
        )
        relative_filepaths.append(relative_filepath)
        # TODO: Delete input fasta?

    _write_last_step(Step.TRINITY_INPUT.value)
    return relative_filepaths


def generate_relative_filepaths_for_trinity(merged_sets_len: int) -> List[Path]:
    return [
        Path(f"{str(set_index)}/{str(set_index)}.fa")
        for set_index in range(merged_sets_len)
    ]


def generate_trinity_results(
    min_contig_length: int, relative_filepaths: List[Path]
) -> List[Path]:
    trinity_results_paths = TrinityWrapper(
        PATH_TO_TEMP_FILES, min_contig_length
    ).write_trinity_result_for_multiple_files(relative_filepaths)
    TrinityWrapper.save_results_filepaths(trinity_results_paths)
    _write_last_step(Step.TRINITY_RUN.value)
    return trinity_results_paths


def load_trinity_results_paths() -> List[Path]:
    return TrinityWrapper.load_results_filepaths()


def align_reads_to_genome(
    trinity_results_paths: List[Path],
    genome: Union[str, None],
    bowtie_index: Union[str, None],
) -> None:
    bowtie_wrapper = BowtieWrapper()
    # TODO: How to automatically use previously generated index?
    if bowtie_index is not None and bowtie_index != "":
        bowtie_wrapper.use_existing_index(bowtie_index)
    elif genome is not None and genome != "":
        bowtie_wrapper.build_index(
            [Path(index_file) for index_file in genome.split(",")]
        )
    else:
        missing_option_message = (
            "One of '-g/--genome', '-x/--bowtie-index' must be set."
        )
        logging.error(missing_option_message)
        raise ValueError(missing_option_message)

    bowtie_wrapper.produce_alignment(
        trinity_results_paths, PATH_TO_TEMP_FILES / "alignment_result.sam"
    )
    _write_last_step(Step.ALIGNMENT.value)


@click.command()
@click.option(
    "-p",
    "--peptide-file",
    required=True,
    type=str,
    help="The path to the peptide file.",
)
@click.option(
    "-r",
    "--rna-file",
    required=True,
    type=str,
    help=(
        "The path to the RNA-seq file. In case of paired-end sequencing, "
        "this file is expected to be in forward orientation."
    ),
)
@click.option(
    "-pa",
    "--paired-end-file",
    required=False,
    type=str,
    default="",
    show_default=False,
    help=(
        "The path to the second RNA-seq file in case of paired-end sequencing. "
        "This file is expected to be in reverse orientation. "
        'If none is given, the RNA-seq file given with the "-r" option '
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
        "The cutoff is applied to all reads. "
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
        "The k-mer size used during the mapping of peptides to RNA-seq reads. "
        "As the RNA-seq reads are 3-frame translated for the mapping, "
        "the k-mer size refers to amino acids."
    ),
)
@click.option(
    "-o",
    "--output-dir",
    required=False,
    type=str,
    default="./",
    help="The path to the output directory for all generated files.",
)
@click.option(
    "-pi",
    "--precompute-intersections",
    is_flag=True,
    help=(
        "If used, the intersection sizes for the Jaccard Index "
        "calculation are precomputed during the matching phase."
    ),
)
@click.option(
    "-j",
    "--jaccard-index-threshold",
    required=False,
    type=float,
    default=0.5,
    help=(
        "Sets of matched RNA-seq reads per peptide will only be merged "
        "together if their Jaccard Index has a value above the given threshold."
    ),
)
@click.option(
    "-m",
    "--merging-method",
    required=False,
    type=click.Choice(["agglomerative-clustering", "full-matrix"]),
    default="full-matrix",
    help="Which merging method to use for sets of matched RNA-seq reads.",
)
@click.option(
    "-cl",
    "--min-contig-length",
    required=False,
    type=int,
    default=100,
    help="Sets the '--min_contig_length' option for Trinity during assembly.",
)
@click.option(
    "-g",
    "--genome",
    required=False,
    type=str,
    help=(
        "The path to the genome file(s) to align to. "
        "In case of multiple files, the paths must be separated by comma."
    ),
)
@click.option(
    "-x",
    "--bowtie-index",
    required=False,
    type=str,
    help=(
        "The basename of an existing Bowtie index that should be used "
        "instead of building a new one. If this option is set, "
        "the '-g/--genome' option is ignored."
    ),
)
def main(
    peptide_file: str,
    rna_file: str,
    paired_end_file: str,
    cutoff: int,
    kmer_length: int,
    output_dir: str,
    precompute_intersections: bool,
    jaccard_index_threshold: float,
    merging_method: Literal["agglomerative-clustering", "full-matrix"],
    min_contig_length: int,
    genome: Union[str, None],
    bowtie_index: Union[str, None],
):
    _setup()
    # TODO: Delete all classes not needed here!

    # TODO: Also need to read in the RNA file to access original sequences via id

    # TODO: Add progress indications

    # TODO: Recreate clean version of requirements.txt

    # TODO: Add full docstrings for all relevant methods

    # TODO: Use numpy arrays everywhere?

    last_step = _get_last_step()
    # TODO: Make more beautiful?
    matches: List[Union[Set[int], None]] = []
    precomputed_intersections: Union[npt.NDArray[np.uint32], None] = None
    if last_step < Step.MATCHING.value:
        logging.info("Computing matches from the given files.")
        matches, precomputed_intersections = compute_matches(
            peptide_file,
            rna_file,
            paired_end_file,
            cutoff,
            kmer_length,
            output_dir,
            precompute_intersections,
        )
    elif last_step == Step.MATCHING.value:
        logging.info("Using already computed matches from last run.")
        matches, precomputed_intersections = load_matches(precompute_intersections)

    if last_step < Step.MERGING.value:
        logging.info("Merging sets of matched RNA-seq reads.")
        merged_sets, peptide_indexes = compute_merge_results(
            jaccard_index_threshold, merging_method, matches, precomputed_intersections
        )
    # TODO: What if merged sets not needed?
    else:
        logging.info("Using already computed merge results from last run.")
        merged_sets, peptide_indexes = load_merge_results()

    if last_step < Step.TRINITY_INPUT.value:
        logging.info("Generating input files for Trinity.")
        relative_filepaths = generate_trinity_input(
            rna_file, paired_end_file, cutoff, merged_sets
        )
    # TODO: What if Trinity step not needed?
    else:
        logging.info("Using already generated Trinity input files.")
        relative_filepaths = generate_relative_filepaths_for_trinity(
            len(peptide_indexes)
        )

    if last_step < Step.TRINITY_RUN.value:
        logging.info("Assembling sets of RNA-seq reads using Trinity.")
        trinity_results_paths = generate_trinity_results(
            min_contig_length, relative_filepaths
        )
    else:
        logging.info("Using already generated Trinity output files.")
        trinity_results_paths = load_trinity_results_paths()

    print(trinity_results_paths)

    if last_step < Step.ALIGNMENT.value:
        logging.info("Aligning assembled RNA-seq reads to the genome.")
        align_reads_to_genome(trinity_results_paths, genome, bowtie_index)

    convert_bam_to_gff(
        PATH_TO_TEMP_FILES / "alignment_result.sam", Path(output_dir) / OUTPUT_FILENAME
    )

    _teardown()


if __name__ == "__main__":
    main()
