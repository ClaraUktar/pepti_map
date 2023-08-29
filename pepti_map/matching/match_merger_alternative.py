from typing import List, Set, Tuple, Union

from pepti_map.util.jaccard_index import jaccard_index


class MatchMerger:
    def __init__(
        self, matches: List[Union[Set[int], None]], jaccard_index_threshold: float = 0.7
    ):
        # TODO: Want to delete matches after merging
        self.matches = matches
        self.jaccard_index_threshold = jaccard_index_threshold
        self._merge_indications: List[List[int]] = [
            [] for _ in range(0, len(self.matches))
        ]
        self.peptide_mappings: List[List[int]] = []
        self.merged_matches: List[Set[int]] = []

    def _process_match_entry(
        self, match_entry: Union[Set[int], None], entry_index: int
    ) -> None:
        if match_entry is None:
            return
        for set_index in range(entry_index + 1, len(self.matches)):
            set_to_compare = self.matches[set_index]
            if set_to_compare is None:
                continue

            if (
                jaccard_index(match_entry, set_to_compare)
                >= self.jaccard_index_threshold
            ):
                self._merge_indications[set_index].append(entry_index)

    def _construct_final_merge_result(self) -> None:
        current_set_index = 0
        offset = 0
        while len(self.matches) > 0:
            current_set = self.matches.pop(0)
            if current_set is None:
                continue

            current_peptide_mappings = [current_set_index]

            for merge_indication in self._merge_indications[current_set_index]:
                set_to_merge = self.merged_matches.pop(merge_indication - offset)
                peptide_mappings_to_merge = self.peptide_mappings.pop(
                    merge_indication - offset
                )
                offset += 1
                current_set.update(set_to_merge)
                current_peptide_mappings.extend(peptide_mappings_to_merge)
                current_peptide_mappings.sort()

            self.merged_matches.append(current_set)
            self.peptide_mappings.append(current_peptide_mappings)

            current_set_index += 1

    def merge_matches(self) -> Tuple[List[Set[int]], List[List[int]]]:
        entry_index = 0
        for entry_index, match_entry in enumerate(self.matches):
            self._process_match_entry(match_entry, entry_index)

        self._construct_final_merge_result()

        return (self.merged_matches, self.peptide_mappings)
