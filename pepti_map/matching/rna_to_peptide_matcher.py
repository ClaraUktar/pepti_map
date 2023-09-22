import csv
import logging
from pathlib import Path
from typing import Generator, List, Set, Union
from pepti_map.constants import PATH_TO_MATCHING_RESULT

from pepti_map.peptide_data.peptide_kmer_index import PeptideKmerIndex
from pepti_map.util.k_mer import split_into_kmer
from pepti_map.util.three_frame_translation import get_three_frame_translations

PEPTIDE_READ_QUANT_FILENAME = "peptide_read_quant.tsv"


class RNAToPeptideMatcher:
    def __init__(
        self,
        kmer_index: PeptideKmerIndex,
        number_of_clusters: int,
        peptide_to_cluster_mapping: List[int],
    ):
        # TODO: We probably want to delete the kmer index after the matching
        self._kmer_index: PeptideKmerIndex = kmer_index
        self._matches: List[Union[Set[int], None]] = [
            None for _ in range(0, number_of_clusters)
        ]
        self._peptide_to_cluster_mapping = peptide_to_cluster_mapping
        # TODO: Use numpy array instead?
        self._matches_per_peptide: List[int] = [
            0 for _ in range(0, len(peptide_to_cluster_mapping))
        ]

    def _process_rna_read_to_kmers(self, sequence: str) -> Generator[str, None, None]:
        # TODO: Exchange all T for U? (inplace?)
        # TODO: Construct reverse complement here or during file reading?
        for translation in get_three_frame_translations(sequence):
            for kmer in split_into_kmer(translation[0], self._kmer_index.kmer_length):
                yield kmer[0]

    def add_peptide_matches_for_rna_read(
        self, rna_read_id: int, rna_read_sequence: str
    ) -> None:
        matched_peptides = set()
        for kmer in self._process_rna_read_to_kmers(rna_read_sequence):
            peptide_matches = self._kmer_index.getEntryForKmer(kmer)
            for match in peptide_matches:
                cluster_match = self._peptide_to_cluster_mapping[match]
                if self._matches[cluster_match] is None:
                    self._matches[cluster_match] = set()
                self._matches[
                    cluster_match
                ].add(  # pyright: ignore[reportOptionalMemberAccess]
                    rna_read_id
                )
                matched_peptides.add(match)
        for match in matched_peptides:
            self._matches_per_peptide[match] += 1

    def get_matches(self) -> List[Union[Set[int], None]]:
        return self._matches

    def _get_file_entry_for_match(self, match_index: int, is_last: bool) -> str:
        match = self._matches[match_index]
        entry = ""
        if match is not None:
            entry = ",".join([str(match_elem) for match_elem in match])
        if not is_last:
            entry += "\n"
        return entry

    def save_matches(self) -> None:
        with open(
            PATH_TO_MATCHING_RESULT, "wt", encoding="utf-8"
        ) as matching_result_file:
            matching_result_file.writelines(
                [
                    self._get_file_entry_for_match(
                        match_index, match_index == len(self._matches) - 1
                    )
                    for match_index in range(0, len(self._matches))
                ]
            )

    def write_peptide_read_quant_file(
        self, dirpath: Path, peptide_sequences: List[str]
    ) -> None:
        try:
            assert len(peptide_sequences) == len(self._peptide_to_cluster_mapping)
        except AssertionError as assertion_error:
            logging.error(
                (
                    "Expected the number of given peptide sequences to be the same "
                    "as used for the mapping, but was different."
                )
            )
            raise assertion_error
        with open(
            dirpath / PEPTIDE_READ_QUANT_FILENAME, "wt", encoding="utf-8"
        ) as peptide_quant_file:
            writer = csv.writer(peptide_quant_file, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["peptide_sequence", "group_id", "n_reads_peptide", "n_reads_group"]
            )
            for peptide_index, peptide_sequence in enumerate(peptide_sequences):
                n_cluster_matches = (
                    len(
                        self._matches[
                            self._peptide_to_cluster_mapping[peptide_index]
                        ]  # pyright: ignore[reportGeneralTypeIssues]
                    )
                    if (
                        self._peptide_to_cluster_mapping[peptide_index] >= 0
                        and self._matches[
                            self._peptide_to_cluster_mapping[peptide_index]
                        ]
                        is not None
                    )
                    else 0
                )
                writer.writerow(
                    [
                        peptide_sequence,
                        self._peptide_to_cluster_mapping[peptide_index],
                        self._matches_per_peptide[peptide_index],
                        n_cluster_matches,
                    ]
                )
