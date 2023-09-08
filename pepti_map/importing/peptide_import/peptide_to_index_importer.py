from pathlib import Path
from typing import Dict, List
from hashlib import sha256

from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.util.k_mer import split_into_kmer

PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE = Path("./temp/peptide_to_cluster_mapping.txt")


class PeptideToIndexImporter:
    def __init__(self, kmer_length: int = 7):
        self.kmer_length: int = kmer_length
        self.kmer_index: PeptideKmerIndex = PeptideKmerIndex(kmer_length)

    def reset(self) -> None:
        self.kmer_index.clear()

    @staticmethod
    def _file_contains_protein_groups(filepath: Path) -> bool:
        contains_protein_groups = False
        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            first_line = peptide_file.readline()
            if "\t" in first_line:
                contains_protein_groups = True
        return contains_protein_groups

    @staticmethod
    def _write_peptide_to_cluster_mapping_file(
        peptide_to_cluster_mapping: List[int],
    ) -> None:
        PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE.parent.mkdir(exist_ok=True)
        with open(
            PATH_PEPTIDE_TO_CLUSTER_MAPPING_FILE, "wt", encoding="utf-8"
        ) as peptide_to_cluster_file:
            peptide_to_cluster_file.writelines(
                [str(cluster_id) + "\n" for cluster_id in peptide_to_cluster_mapping]
            )

    def _process_simple_peptide_file(
        self, filepath: Path, replace_isoleucine=True
    ) -> None:
        peptide_to_cluster_mapping: List[int] = []
        number_of_peptides = 0
        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for line in peptide_file:
                sequence = line.strip()
                if len(sequence) < self.kmer_length:
                    peptide_to_cluster_mapping.append(-1)
                    continue

                if replace_isoleucine:
                    sequence = sequence.replace("I", "L")
                for kmer in split_into_kmer(sequence, self.kmer_length):
                    self.kmer_index.appendToEntryForKmer(kmer[0], number_of_peptides)
                peptide_to_cluster_mapping.append(number_of_peptides)
                number_of_peptides += 1

        # TODO: Do this in a prettier way
        self.kmer_index.number_of_peptides = number_of_peptides

        PeptideToIndexImporter._write_peptide_to_cluster_mapping_file(
            peptide_to_cluster_mapping
        )

    def _process_peptide_file_with_protein_groups(
        self, filepath: Path, replace_isoleucine=True
    ) -> None:
        protein_group_cluster_index: Dict[str, int] = {}
        peptide_to_cluster_mapping: List[int] = []
        number_of_clusters = 0

        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for line in peptide_file:
                peptide, protein_group = line.split("\t")
                peptide = peptide.strip()
                # TODO: Should even too short peptides be taken into account later on
                # if their protein group matches other peptides that are long enough?
                if len(peptide) < self.kmer_length:
                    peptide_to_cluster_mapping.append(-1)
                    continue

                protein_group = sha256(
                    protein_group.strip().encode("utf-8")
                ).hexdigest()
                if protein_group in protein_group_cluster_index:
                    cluster_id = protein_group_cluster_index[protein_group]
                else:
                    cluster_id = number_of_clusters
                    number_of_clusters += 1
                    protein_group_cluster_index[protein_group] = cluster_id
                peptide_to_cluster_mapping.append(cluster_id)

                if replace_isoleucine:
                    peptide = peptide.replace("I", "L")
                for kmer in split_into_kmer(peptide, self.kmer_length):
                    self.kmer_index.appendToEntryForKmer(
                        kmer[0], cluster_id, check_for_duplicates=True
                    )

        # TODO: Do this in a prettier way
        self.kmer_index.number_of_peptides = number_of_clusters

        PeptideToIndexImporter._write_peptide_to_cluster_mapping_file(
            peptide_to_cluster_mapping
        )

    def import_file_to_index(
        self, filepath: Path, replace_isoleucine=True
    ) -> PeptideKmerIndex:
        # TODO: Should reset not happen automatically
        # to enable import of multiple files
        self.reset()

        if PeptideToIndexImporter._file_contains_protein_groups(filepath):
            self._process_peptide_file_with_protein_groups(filepath, replace_isoleucine)
        else:
            self._process_simple_peptide_file(filepath, replace_isoleucine)

        return self.kmer_index
