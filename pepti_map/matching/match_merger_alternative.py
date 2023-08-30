from typing import Dict, List, Set, Tuple, Union

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
        self._updated_indexes: Dict[int, int] = {}

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
                self._merge_indications[entry_index].append(set_index)
                self._merge_indications[set_index].append(entry_index)

    def _update_set_with_merges(
        self, set_to_update: set, set_index: int, previous_indexes: List[int]
    ) -> Tuple[Set, List[int]]:
        merge_indications = self._merge_indications[set_index]
        previous_indexes.append(set_index)
        self.matches[set_index] = None
        if len(merge_indications) == 0:
            return (set_to_update, previous_indexes)
        for index_to_merge in merge_indications:
            set_to_merge = self.matches[index_to_merge]
            if set_to_merge is None:
                continue
            updated_set, previous_indexes = self._update_set_with_merges(
                set_to_merge, index_to_merge, previous_indexes
            )
            set_to_update.update(updated_set)
        return (set_to_update, previous_indexes)

    def _construct_final_merge_result(self) -> None:
        for current_index, current_set in enumerate(self.matches):
            if current_set is None:
                continue

            merged_set, merged_indexes = self._update_set_with_merges(
                current_set, current_index, []
            )
            merged_indexes.sort()
            self.merged_matches.append(merged_set)
            self.peptide_mappings.append(merged_indexes)

    def merge_matches(self) -> Tuple[List[Set[int]], List[List[int]]]:
        for entry_index, match_entry in enumerate(self.matches):
            self._process_match_entry(match_entry, entry_index)

        self._construct_final_merge_result()

        return (self.merged_matches, self.peptide_mappings)
