from collections import defaultdict
import logging
import multiprocessing
import os
import gffutils
from pathlib import Path
from typing import List, TextIO, Tuple

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
    def _write_new_feature_coordinates(
        output_gtf: TextIO,
        gene: gffutils.Feature,
        mrna: gffutils.Feature,
        exons: List[gffutils.Feature],
        cds: List[gffutils.Feature],
        start_exon_index: int,
        end_exon_index: int,
        strand: str,
        direction: str,
        contig_length: int,
    ) -> None:
        start_exon = exons[start_exon_index]
        end_exon = exons[end_exon_index]
        if (
            not start_exon.start
            or not start_exon.end
            or not end_exon.start
            or not end_exon.end
        ):
            raise ValueError("No start and end coordinates given.")

        if direction == ".":
            contig_id, start_exon_contig_end, contig_start, _ = start_exon.attributes[
                "Target"
            ][0].split(" ")
            _, contig_end, end_exon_contig_start, _ = end_exon.attributes["Target"][
                0
            ].split(" ")
            contig_start = int(contig_start)
            contig_end = int(contig_end)

            if start_exon == end_exon:
                start_exon.attributes["Target"] = " ".join(
                    [contig_id, str(contig_length), "1", direction]
                )
            else:
                start_exon.attributes["Target"] = " ".join(
                    [contig_id, start_exon_contig_end, "1", direction]
                )
                end_exon.attributes["Target"] = " ".join(
                    [contig_id, str(contig_length), end_exon_contig_start, direction]
                )
        else:
            contig_id, contig_start, start_exon_contig_end, _ = start_exon.attributes[
                "Target"
            ][0].split(" ")
            _, end_exon_contig_start, contig_end, _ = end_exon.attributes["Target"][
                0
            ].split(" ")
            contig_start = int(contig_start)
            contig_end = int(contig_end)

            if start_exon == end_exon:
                start_exon.attributes["Target"] = " ".join(
                    [contig_id, "1", str(contig_length), direction]
                )
            else:
                start_exon.attributes["Target"] = " ".join(
                    [contig_id, "1", start_exon_contig_end, direction]
                )
                end_exon.attributes["Target"] = " ".join(
                    [contig_id, end_exon_contig_start, str(contig_length), direction]
                )

        exon_ids = [exon.attributes["ID"][0] for exon in exons]
        cds_ids = [cds_entry.attributes["ID"][0] for cds_entry in cds]
        mrna_id = mrna.attributes["ID"][0]

        if (
            (direction == "+" and strand == "+")
            or (direction == "-" and strand == "-")
            or (direction == "." and strand == "+")
        ):
            new_start = start_exon.start - (contig_start - 1)
            start_exon.start = new_start
            mrna.start = new_start
            gene.start = new_start
            new_end = end_exon.end + (contig_length - contig_end)
            end_exon.end = new_end
            mrna.end = new_end
            gene.end = new_end

            output_gtf.write(str(gene) + "\n")
            for frame in range(3):
                mrna.attributes["ID"] = mrna_id + "." + str(frame)
                output_gtf.write(str(mrna) + "\n")
                for exon_index, exon in enumerate(exons):
                    exon.attributes["ID"] = exon_ids[exon_index] + "." + str(frame)
                    exon.attributes["Parent"] = mrna.attributes["ID"]
                    output_gtf.write(str(exon) + "\n")
                for cds_index, cds_entry in enumerate(cds):
                    cds_entry.start = exons[cds_index].start
                    cds_entry.end = exons[cds_index].end
                    if cds_index == start_exon_index:
                        cds_entry.start = cds_entry.start + frame
                    if cds_index == end_exon_index:
                        cds_entry.end = cds_entry.end - ((contig_length - frame) % 3)
                    cds_entry.attributes["ID"] = cds_ids[cds_index] + "." + str(frame)
                    cds_entry.attributes["Parent"] = mrna.attributes["ID"]
                    output_gtf.write(str(cds_entry) + "\n")

        elif (
            (direction == "+" and strand == "-")
            or (direction == "-" and strand == "+")
            or (direction == "." and strand == "-")
        ):
            new_end = start_exon.end + (contig_start - 1)
            start_exon.end = new_end
            mrna.end = new_end
            gene.end = new_end
            new_start = end_exon.start - (contig_length - contig_end)
            end_exon.start = new_start
            mrna.start = new_start
            gene.start = new_start

            output_gtf.write(str(gene) + "\n")
            for frame in range(3):
                mrna.attributes["ID"] = mrna_id + "." + str(frame)
                output_gtf.write(str(mrna) + "\n")
                for exon_index, exon in enumerate(exons):
                    exon.attributes["ID"] = exon_ids[exon_index] + "." + str(frame)
                    exon.attributes["Parent"] = mrna.attributes["ID"]
                    output_gtf.write(str(exon) + "\n")
                for cds_index, cds_entry in enumerate(cds):
                    cds_entry.start = exons[cds_index].start
                    cds_entry.end = exons[cds_index].end
                    if cds_index == start_exon_index:
                        cds_entry.end = cds_entry.end - frame
                    if cds_index == end_exon_index:
                        cds_entry.start = cds_entry.start + (
                            (contig_length - frame) % 3
                        )
                    cds_entry.attributes["ID"] = cds_ids[cds_index] + "." + str(frame)
                    cds_entry.attributes["Parent"] = mrna.attributes["ID"]
                    output_gtf.write(str(cds_entry) + "\n")
        else:
            raise ValueError("Strand must be one of '+', '-'.")

    @classmethod
    def generate_gtf_input_file(
        cls,
        path_to_gff: Path,
        output_directory: Path,
        sequence_lengths_per_contig: List[int],
    ) -> List[int]:
        # Track number of transcripts to write protein FASTA with matching ids
        number_of_transcripts_per_contig: List[int] = [
            0 for _ in range(len(sequence_lengths_per_contig))
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
                exons = [
                    gene_child
                    for gene_child in gene_children
                    if gene_child.featuretype == "exon"
                ]

                first_exon = exons[0]
                strand = first_exon.strand
                target: str = first_exon.attributes["Target"][0]
                contig_id, _, _, direction = target.split(" ")
                contig_length = sequence_lengths_per_contig[int(contig_id[-1])]
                if direction == ".":  # indeterminate
                    start_exon_index = exons.index(
                        min(
                            exons,
                            key=lambda exon: int(
                                exon.attributes["Target"][0].split(" ")[2]
                            ),
                        )
                    )
                    end_exon_index = exons.index(
                        max(
                            exons,
                            key=lambda exon: int(
                                exon.attributes["Target"][0].split(" ")[1]
                            ),
                        )
                    )
                else:  # sense or antisense
                    start_exon_index = exons.index(
                        min(
                            exons,
                            key=lambda exon: int(
                                exon.attributes["Target"][0].split(" ")[1]
                            ),
                        )
                    )
                    end_exon_index = exons.index(
                        max(
                            exons,
                            key=lambda exon: int(
                                exon.attributes["Target"][0].split(" ")[2]
                            ),
                        )
                    )
                mrna = [
                    gene_child
                    for gene_child in gene_children
                    if gene_child.featuretype == "mRNA"
                ][
                    0
                ]  # There can be only one mRNA per gene
                mrna_id = mrna.attributes["ID"][0]
                number_of_transcripts_per_contig[int(mrna_id.split(".")[0][-1])] += 1

                cls._write_new_feature_coordinates(
                    output_gtf,
                    gene_feature,
                    mrna,
                    exons,
                    [
                        gene_child
                        for gene_child in gene_children
                        if gene_child.featuretype == "CDS"
                    ],
                    start_exon_index,
                    end_exon_index,
                    strand,
                    direction,
                    contig_length,
                )

        return number_of_transcripts_per_contig

    @staticmethod
    def generate_protein_fasta_input_file(
        contig_sequences: List[Tuple[str, str]],
        output_directory: Path,
        number_of_transcripts_per_contig: List[int],
    ) -> None:
        with open(
            output_directory / "pogo_fasta_in.fa", "wt", encoding="utf-8"
        ) as output_file:
            for contig_index, (contig_id, contig_sequence) in enumerate(
                contig_sequences
            ):
                for translation, frame in get_three_frame_translations(
                    contig_sequence, False
                ):
                    for transcript_index in range(
                        number_of_transcripts_per_contig[contig_index]
                    ):
                        gene_id = f"{contig_id}-path{str(transcript_index + 1)}"
                        transcript_id = (
                            f"{contig_id}-mrna{str(transcript_index + 1)}-{str(frame)}"
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
    def generate_gff_and_protein_files_for_directory(
        cls, path_to_directory: Path
    ) -> None:
        contig_sequences = cls._get_contig_sequences(
            path_to_directory / "resulting_contigs.fa"
        )
        number_of_transcripts_per_contig = cls.generate_gtf_input_file(
            path_to_directory / "alignment_result.gff3",
            path_to_directory,
            [len(contig_sequence[1]) for contig_sequence in contig_sequences],
        )
        cls.generate_protein_fasta_input_file(
            contig_sequences,
            path_to_directory,
            number_of_transcripts_per_contig,
        )

    @classmethod
    def generate_gff_and_protein_files_for_multiple_directories(
        cls,
        paths_to_directories: List[Path],
    ) -> None:
        # TODO: Refactor code duplication
        try:
            n_processes = os.getenv("IO_N_PROCESSES")
            assert isinstance(n_processes, str)
            n_processes = int(n_processes)
        except (AssertionError, ValueError):
            n_processes = multiprocessing.cpu_count()
        logging.info(
            f"Generating GTF and FASTA files for PoGo with {n_processes} processes"
        )
        with multiprocessing.Pool(n_processes) as pool:
            pool.map(
                cls.generate_gff_and_protein_files_for_directory, paths_to_directories
            )
