from typing import Dict, List
from hashlib import sha256

from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.util.k_mer import split_into_kmer


class PeptideToIndexImporter:
    def __init__(self, kmer_length: int = 7):
        self.kmer_length: int = kmer_length
        self.kmer_index: PeptideKmerIndex = PeptideKmerIndex(kmer_length)

    def reset(self) -> None:
        self.kmer_index.clear()

    @staticmethod
    def _file_contains_protein_groups(filepath: str) -> bool:
        contains_protein_groups = False
        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            first_line = peptide_file.readline()
            if "\t" in first_line:
                contains_protein_groups = True
        return contains_protein_groups

    def _process_simple_peptide_file(
        self, filepath: str, replace_isoleucine=True
    ) -> None:
        # TODO: We could throw out the too short peptides here as well
        # Then need to adapt count of peptides
        # and how they are treated during matching!
        # e.g. write mapping id -> internal index to file?
        number_of_peptides = 0
        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for index, line in enumerate(peptide_file):
                sequence = line.strip()
                if replace_isoleucine:
                    sequence = sequence.replace("I", "L")
                for kmer in split_into_kmer(sequence, self.kmer_length):
                    self.kmer_index.appendToEntryForKmer(kmer[0], index)
                if sequence != "":
                    number_of_peptides += 1

        # TODO: Do this in a prettier way
        self.kmer_index.number_of_peptides = number_of_peptides

    def _process_peptide_file_with_protein_groups(
        self, filepath: str, replace_isoleucine=True
    ) -> None:
        protein_group_cluster_index: Dict[str, int] = {}
        cluster_mappings: List[List[int]] = []

        with open(filepath, "rt", encoding="utf-8") as peptide_file:
            for index, line in enumerate(peptide_file):
                peptide, protein_group = line.split("\t")
                peptide = peptide.strip()
                # TODO: Should even too short peptides be taken into account later on
                # if their protein group matches other peptides that are long enough?
                if len(peptide) < self.kmer_length:
                    continue
                protein_group = sha256(
                    protein_group.strip().encode("utf-8")
                ).hexdigest()
                if protein_group in protein_group_cluster_index:
                    cluster_id = protein_group_cluster_index[protein_group]
                else:
                    cluster_id = len(cluster_mappings)
                    cluster_mappings.append([])
                    protein_group_cluster_index[protein_group] = cluster_id
                cluster_mappings[cluster_id].append(index)

                if replace_isoleucine:
                    peptide = peptide.replace("I", "L")
                for kmer in split_into_kmer(peptide, self.kmer_length):
                    self.kmer_index.appendToEntryForKmer(
                        kmer[0], cluster_id, check_for_duplicates=True
                    )

        # TODO: Do this in a prettier way
        self.kmer_index.number_of_peptides = len(cluster_mappings)

        # TODO: Store cluster_mappings/write to file for later use!!!

    def import_file_to_index(
        self, filepath: str, replace_isoleucine=True
    ) -> PeptideKmerIndex:
        # TODO: Should reset not happen automatically
        # to enable import of multiple files
        self.reset()

        if self._file_contains_protein_groups(filepath):
            self._process_peptide_file_with_protein_groups(filepath, replace_isoleucine)
        else:
            self._process_simple_peptide_file(filepath, replace_isoleucine)

        return self.kmer_index
