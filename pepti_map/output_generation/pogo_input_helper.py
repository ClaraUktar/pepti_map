from collections import defaultdict
from functools import partial
import logging
import multiprocessing
from dotenv import dotenv_values
import gffutils
from pathlib import Path
from typing import List, TextIO, Tuple, Union

from pepti_map.util.three_frame_translation import get_three_frame_translations


class PoGoInputHelper:
    def __init__(
        self, path_to_peptides: Path, path_to_peptide_to_cluster_mapping: Path
    ):
        self._peptides = []
        with open(path_to_peptides, "rt", encoding="utf-8") as peptides_file:
            for line in peptides_file:
                self._peptides.append(
                    line.split("\t")[0].strip()
                )  # ignores group information if present

        # Reverse mapping to get cluster_id -> peptide_id
        self._cluster_to_peptide_mapping: "defaultdict[int, List[int]]" = defaultdict(
            list
        )
        with open(
            path_to_peptide_to_cluster_mapping, "rt", encoding="utf-8"
        ) as peptide_to_cluster_mapping:
            for peptide_id, line in enumerate(peptide_to_cluster_mapping):
                cluster_id = int(line.strip())
                # Exclude peptides that were too short to be included
                if cluster_id == -1:
                    continue
                self._cluster_to_peptide_mapping[cluster_id].append(peptide_id)

    def generate_peptide_input_file(
        self, output_directory: Path, merged_indexes: List[int]
    ) -> None:
        peptide_already_written: List[bool] = [
            False for _ in range(len(self._peptides))
        ]
        with open(
            output_directory / "pogo_peptides_in.tsv", "wt", encoding="utf-8"
        ) as new_input_file:
            for cluster_id in merged_indexes:
                peptide_ids = self._cluster_to_peptide_mapping[cluster_id]
                for peptide_id in peptide_ids:
                    if peptide_already_written[peptide_id]:
                        continue
                    # Use 1 as default for Sample, PSMs and Quant
                    new_input_file.write(
                        "\t".join(["1", self._peptides[peptide_id], "1", "1"]) + "\n"
                    )
                    peptide_already_written[peptide_id] = True

    def generate_all_peptide_input_files(
        self, output_directories: List[Path], path_to_merged_indexes: Path
    ) -> None:
        # TODO: Refactor to use code from match merger
        merged_indexes: List[List[int]] = []
        with open(
            path_to_merged_indexes, "rt", encoding="utf-8"
        ) as peptide_indexes_file:
            for line in peptide_indexes_file:
                line = line.strip()
                merged_indexes.append(
                    [int(match_elem) for match_elem in line.split(",")]
                )

        for output_directory in output_directories:
            # TODO: Refactor?
            set_index = int(output_directory.name)
            self.generate_peptide_input_file(
                output_directory, merged_indexes[set_index]
            )

    @staticmethod
    def _get_exon_feature_id(exon: gffutils.Feature) -> int:
        return int(exon.attributes["ID"][0].split(".")[-1].replace("exon", ""))

    @staticmethod
    def _write_feature_in_gtf_format(
        output_gtf: TextIO,
        feature: gffutils.Feature,
        gene_id: str,
        transcript_id: Union[str, None] = None,
    ) -> None:
        featuretype = feature.featuretype
        if featuretype == "mRNA":
            featuretype = "transcript"
        attributes = f"gene_id \"{gene_id.replace('.', '_')}\";"
        if transcript_id is not None:
            attributes = (
                attributes + f"transcript_id \"{transcript_id.replace('.', '_')}\";"
            )

        output_gtf.write(
            "\t".join(
                [
                    feature.chrom,  # pyright: ignore[reportGeneralTypeIssues]
                    feature.source,
                    featuretype,
                    str(feature.start),
                    str(feature.end),
                    ".",
                    feature.strand,
                    ".",
                    attributes,
                ]
            )
            + "\n"
        )

    @classmethod
    def _write_new_feature_coordinates(
        cls,
        output_gtf: TextIO,
        gene: gffutils.Feature,
        mrna: gffutils.Feature,
        exons: List[gffutils.Feature],
        strand: str,
        direction: str,
        contig: str,
    ) -> str:
        # The exons need to be sorted to follow the same order as in the original GFF.
        # This is necessary because gffutils does not necessarily return the children
        # of a feature in order when calling children(). The parameter order_by
        # cannot be used here because it does not allow us to select for an
        # exon identifier or specify the same order as in the original.
        # TODO: Is there an easier option for ordering
        # that does not involve string splitting?
        exons.sort(key=cls._get_exon_feature_id)

        # If direction is antisense, change alignment to be on the complementary strand
        if direction == "-":
            if strand == "+":
                strand = "-"
            else:
                strand = "+"

            gene.strand = strand
            mrna.strand = strand
            for exon in exons:
                exon.strand = strand

            exons.reverse()

        # TODO: Is assumption that exons are already in correct order
        # if dir=indeterminate true?
        # Potentially only labeled "indeterminate" if only one exon?

        exon_starts: List[int] = []
        exon_ends: List[int] = []
        for exon in exons:
            if direction == ".":
                _, exon_end, exon_start, _ = exon.attributes["Target"][0].split(" ")
            else:
                _, exon_start, exon_end, _ = exon.attributes["Target"][0].split(" ")
            exon_starts.append(int(exon_start))
            exon_ends.append(int(exon_end))

        cut_contig_parts: List[str] = []
        for i in range(len(exon_starts)):
            current_start = exon_starts[i]
            current_end = exon_ends[i]
            current_part = contig[current_start - 1 : current_end]  # noqa: E203
            cut_contig_parts.append(current_part)
        start_end_cut_contig = "".join(cut_contig_parts)

        if len(start_end_cut_contig) < (0.7 * len(contig)):
            # TODO: Filter out alignment and report
            pass

        # TODO: Not needed?
        # exon_ids = [exon.attributes["ID"][0] for exon in exons]
        # cds_ids = [exon_id.replace("exon", "cds") for exon_id in exon_ids]
        mrna_id = mrna.attributes["ID"][0]
        exon_start_coords = [exon.start for exon in exons]
        exon_end_coords = [exon.end for exon in exons]

        cls._write_feature_in_gtf_format(
            output_gtf,
            gene,
            gene.attributes["ID"][0],
        )
        for frame in range(3):
            mrna.attributes["ID"] = mrna_id + "." + str(frame)
            cls._write_feature_in_gtf_format(
                output_gtf, mrna, gene.attributes["ID"][0], mrna.attributes["ID"][0]
            )
            for exon_idx, exon in enumerate(exons):
                # TODO: Not needed?
                # exon.attributes["ID"] = exon_ids[exon_index] + "." + str(frame)
                # exon.attributes["Parent"] = mrna.attributes["ID"]
                exon.featuretype = "exon"
                exon.start = exon_start_coords[exon_idx]
                exon.end = exon_end_coords[exon_idx]
                cls._write_feature_in_gtf_format(
                    output_gtf,
                    exon,
                    gene.attributes["ID"][0],
                    mrna.attributes["ID"][0],
                )
            for exon_idx, exon in enumerate(exons):
                # TODO: Is there a better solution,
                # e.g. copying and modifying the feature?
                exon.featuretype = "CDS"
                # strand = +, dir = sense
                # -> add frame to start of first CDS
                # strand = +, dir = antisense
                # -> subtract frame from end of first CDS (is first after reversing)
                # strand = -, dir = sense
                # -> subtract frame from end of first CDS
                # strand = -, dir = antisense
                # -> add frame to start of first CDS (is first after reversing)
                # --> differentiation between +/- strand should suffice after reversing
                if strand == "+":
                    if exon_idx == 0:
                        exon.start = exon.start + frame  # pyright: ignore
                    if exon_idx == (len(exons) - 1):
                        exon.end = exon.end - (  # pyright: ignore
                            (len(start_end_cut_contig) - frame) % 3
                        )
                else:
                    if exon_idx == 0:
                        exon.end = exon.end - frame  # pyright: ignore
                    if exon_idx == (len(exons) - 1):
                        exon.start = exon.start + (  # pyright: ignore
                            (len(start_end_cut_contig) - frame) % 3
                        )

                cls._write_feature_in_gtf_format(
                    output_gtf,
                    exon,
                    gene.attributes["ID"][0],
                    mrna.attributes["ID"][0],
                )

        return start_end_cut_contig

    @classmethod
    def generate_gtf_input_file(
        cls,
        path_to_gff: Path,
        output_directory: Path,
        contig_sequences: List[Tuple[str, str]],
        no_indels=False,
    ) -> Tuple[List[List[int]], List[List[str]]]:
        # Track contig alignment ids to write protein FASTA with matching ids
        alignment_ids_per_contig: List[List[int]] = [
            [] for _ in range(len(contig_sequences))
        ]
        # Per original contig, there can be several new contigs
        # based on different cutoffs
        new_contig_sequences: List[List[str]] = [
            [] for _ in range(len(contig_sequences))
        ]

        gffutils_db = gffutils.create_db(
            path_to_gff.absolute().as_posix(),
            (path_to_gff.parent / "gffutils_db.sqlite").absolute().as_posix(),
        )

        with open(
            output_directory / "pogo_gtf_in.gtf", "wt", encoding="utf-8"
        ) as output_gtf:
            for gene_feature in gffutils_db.features_of_type("gene"):
                gene_children = list(gffutils_db.children(gene_feature))

                mrna = [
                    gene_child
                    for gene_child in gene_children
                    if gene_child.featuretype == "mRNA"
                ][
                    0
                ]  # There can be only one mRNA per gene

                # If no indels allowed: Check if feature contains indels
                # -> If so, gene feature is skipped
                if no_indels:
                    mrna_indels = int(mrna.attributes["indels"][0])
                    if mrna_indels != 0:
                        continue

                exons = [
                    gene_child
                    for gene_child in gene_children
                    if gene_child.featuretype == "exon"
                ]

                first_exon = exons[0]
                strand = first_exon.strand
                target: str = first_exon.attributes["Target"][0]
                contig_id, _, _, direction = target.split(" ")
                contig_id = int(contig_id.split("-")[-1])
                contig = contig_sequences[contig_id]

                new_contig = cls._write_new_feature_coordinates(
                    output_gtf,
                    gene_feature,
                    mrna,
                    exons,
                    strand,
                    direction,
                    contig[1],
                )

                mrna_id = mrna.attributes["ID"][0]
                contig_idx = int(mrna_id.split(".")[0].split("-")[-1])
                path_number = int(mrna_id.split(".")[1].replace("mrna", ""))
                new_contig_sequences[contig_idx].append(new_contig)
                alignment_ids_per_contig[contig_idx].append(path_number)

        return (alignment_ids_per_contig, new_contig_sequences)

    # TODO: Remove unneeded arguments
    @staticmethod
    def generate_protein_fasta_input_file(
        contig_ids: List[str],
        contig_sequences: List[List[str]],
        output_directory: Path,
        alignment_ids_per_contig: List[List[int]],
    ) -> None:
        # TODO: Adapt to new separation of ids and seqs
        with open(
            output_directory / "pogo_fasta_in.fa", "wt", encoding="utf-8"
        ) as output_file:
            for contig_id, contig_cut_sequences, alignment_ids in zip(
                contig_ids, contig_sequences, alignment_ids_per_contig
            ):
                for alignment_id, contig_sequence in zip(
                    alignment_ids, contig_cut_sequences
                ):
                    for translation, frame in get_three_frame_translations(
                        contig_sequence, False
                    ):
                        gene_id = f"{contig_id}_path{str(alignment_id)}"
                        transcript_id = (
                            f"{contig_id}_mrna{str(alignment_id)}_{str(frame)}"
                        )
                        output_file.write(
                            (
                                f">{contig_id} gene:{gene_id} "
                                f"transcript:{transcript_id}\n"
                            )
                        )
                        output_file.write(translation + "\n")

    @staticmethod
    def _get_contig_sequences(path_to_contig_sequences: Path) -> List[Tuple[str, str]]:
        contig_sequences: List[Tuple[str, str]] = []
        with open(
            path_to_contig_sequences, "rt", encoding="utf-8"
        ) as contig_sequences_file:
            current_id = ""
            for line_index, line in enumerate(contig_sequences_file):
                if line_index % 2 == 1:
                    contig_sequences.append((current_id, line.strip()))
                else:
                    current_id = line.strip().replace(">", "")
        return contig_sequences

    @classmethod
    def generate_gtf_and_protein_files_for_directory(
        cls, path_to_directory: Path, no_indels=False
    ) -> None:
        contig_sequences = cls._get_contig_sequences(
            path_to_directory / "resulting_contigs.fa"
        )
        alignment_ids_per_contig, updated_contig_sequences = (
            cls.generate_gtf_input_file(
                path_to_directory / "alignment_result.gff3",
                path_to_directory,
                contig_sequences,
                no_indels,
            )
        )
        cls.generate_protein_fasta_input_file(
            [contig_sequence[0] for contig_sequence in contig_sequences],
            updated_contig_sequences,
            path_to_directory,
            alignment_ids_per_contig,
        )

    @classmethod
    def generate_gtf_and_protein_files_for_multiple_directories(
        cls, paths_to_directories: List[Path], no_indels=False
    ) -> None:
        # TODO: Refactor code duplication
        try:
            n_processes = dotenv_values().get("IO_N_PROCESSES")
            assert isinstance(n_processes, str)
            n_processes = int(n_processes)
        except (AssertionError, ValueError):
            n_processes = multiprocessing.cpu_count()
        logging.info(
            f"Generating GTF and FASTA files for PoGo with {n_processes} processes"
        )
        with multiprocessing.Pool(n_processes) as pool:
            pool.map(
                partial(
                    cls.generate_gtf_and_protein_files_for_directory,
                    no_indels=no_indels,
                ),
                paths_to_directories,
            )
