from collections import defaultdict
import logging
import multiprocessing
import os
import gffutils
from pathlib import Path
from typing import List, Tuple

from pepti_map.util.three_frame_translation import get_three_frame_translations


class PepGenomeInputHelper:
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
        self, output_path: Path, merged_indexes: List[int]
    ) -> None:
        peptide_already_written: List[bool] = [
            False for _ in range(len(self._peptides))
        ]
        with open(output_path, "wt", encoding="utf-8") as new_input_file:
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
        self, output_paths: List[Path], path_to_merged_indexes: Path
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

        for set_index, output_path in enumerate(output_paths):
            self.generate_peptide_input_file(output_path, merged_indexes[set_index])

    @staticmethod
    def _calculate_new_feature_coordinates(
        contig_start: int,
        contig_end: int,
        feature_start: int,
        feature_end: int,
        contig_length: int,
        strand: str,
        direction: str,
    ) -> Tuple[int, int]:
        # TODO: Add algorithm
        pass

    @classmethod
    def generate_gff_input_file(
        cls,
        path_to_gff: Path,
        output_path: Path,
        sequence_lengths_per_contig: List[int],
    ) -> List[int]:
        # TODO: Return number of transcripts per contig (?)
        # TODO: For each entry in GFF, adapt positions
        gffutils_db = gffutils.create_db(
            path_to_gff.absolute().as_posix(),
            (path_to_gff.parent / "gffutils_db.sqlite").absolute().as_posix(),
        )
        # TODO: Add counts
        number_of_transcripts_per_contig: List[int] = [
            0 for _ in range(len(sequence_lengths_per_contig))
        ]
        for gene_feature in gffutils_db.features_of_type("gene"):
            gene_children = list(gffutils_db.children(gene_feature))
            # Delete unneeded CDS features
            for gene_child in gene_children:
                if gene_child.featuretype == "CDS":
                    gffutils_db.delete(gene_child.id, False)
            exons = [
                gene_child
                for gene_child in gene_children
                if gene_child.featuretype == "exon"
            ]
            if len(exons) == 1:
                only_exon = exons[0]
                strand = only_exon.strand
                target: str = only_exon.attributes["Target"][0]
                contig_id, contig_start, contig_end, direction = target.split(" ")
                contig_length = sequence_lengths_per_contig[int(contig_id[-1])]
                contig_start = int(contig_start)
                contig_end = int(contig_end)
                if contig_start < contig_end:
                    only_exon.attributes["Target"] = " ".join(
                        [contig_id, "1", str(contig_length), direction]
                    )
                else:
                    only_exon.attributes["Target"] = " ".join(
                        [contig_id, str(contig_length), "1", direction]
                    )
                new_start, new_end = cls._calculate_new_feature_coordinates(
                    contig_start,
                    contig_end,
                    only_exon.start,  # pyright: ignore[reportGeneralTypeIssues]
                    only_exon.end,  # pyright: ignore[reportGeneralTypeIssues]
                    contig_length,
                    strand,
                    direction,
                )
                only_exon.start = new_start
                only_exon.end = new_end
            else:
                # TODO: Case of multiple exons
                pass

            mrna = [
                gene_child
                for gene_child in gene_children
                if gene_child.featuretype == "mRNA"
            ][
                0
            ]  # There can be only one mRNA per gene
            mrna.start = new_start
            mrna.end = new_end
            gene_feature.start = new_start
            gene_feature.end = new_end

        with open(output_path, "wt", encoding="utf-8") as output_gff:
            for feature in gffutils_db.all_features():
                output_gff.write(str(feature) + "\n")

        # Idea:
        # 1. Iterate over all "gene" features
        # 2. For each, get all children
        # 2a. Delete all CDS features (?)
        # 3. Look at exon attributes first to determine new coordinates
        # 4. Adapt coordinates of exon + start/end in attrs
        # - need to be able to deal with multiple exons
        # 5. Adapt coordinates of mRNA and gene
        # 6. Iterate over whole DB to generate new output file

    @staticmethod
    def generate_protein_fasta_input_file(
        path_to_contig_sequences: Path,
        output_path: Path,
        number_of_transcripts_per_contig: List[int],
    ) -> None:
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

        # TODO: Check if correct format
        with open(output_path, "wt", encoding="utf-8") as output_file:
            for contig_index, (contig_id, contig_sequence) in enumerate(
                contig_sequences
            ):
                for translation, frame in get_three_frame_translations(
                    contig_sequence, False
                ):
                    for transcript_index in range(
                        number_of_transcripts_per_contig[contig_index]
                    ):
                        gene_id = f"{contig_id}.path{transcript_index + 1}"
                        transcript_id = f"{contig_id}.mrna{transcript_index + 1}"
                        output_file.write(
                            (
                                f">{contig_id} geneID={gene_id} "
                                f"transcriptID={transcript_id} offset={str(frame)}\n"
                            )
                        )
                        output_file.write(translation + "\n")

    @classmethod
    def generate_gff_and_protein_files_for_directory(
        cls, path_to_directory: Path
    ) -> None:
        number_of_transcripts_per_contig = cls.generate_gff_input_file(
            path_to_directory / "alignment_result.gff3",
            path_to_directory / "pepgenome_gff_in.gff3",
        )
        cls.generate_protein_fasta_input_file(
            path_to_directory / "resulting_contigs.fa",
            path_to_directory / "pepgenome_fasta_in.fa",
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
            f"Generating GFF and FASTA files for PepGenome with {n_processes} processes"
        )
        with multiprocessing.Pool(n_processes) as pool:
            pool.map(
                cls.generate_gff_and_protein_files_for_directory, paths_to_directories
            )
