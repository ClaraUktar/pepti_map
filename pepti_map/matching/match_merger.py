from typing import List, Set, Union

from pepti_map.util.jaccard_index import jaccard_index


class MatchMerger:
    def __init__(
        self, matches: List[Union[Set[int], None]], jaccard_index_threshold: float = 0.8
    ):
        # TODO: Want to delete matches after merging
        self.matches = matches
        self.jaccard_index_threshold = jaccard_index_threshold
        # TODO: Want to pop this one as well I think
        # (otherwise no benefit over computing all merge-matches at once?)
        self._merge_indications: List[List[int]] = [
            [] for _ in range(0, len(self.matches))
        ]
        self._peptide_mappings: List[List[int]] = [
            [] for _ in range(0, len(self.matches))
        ]
        self.merged_matches: List[Set[int]] = []

    def _process_match_entry(
        self, match_entry: Union[Set[int], None], entry_index: int
    ):
        if match_entry is None:
            return
        extended_set = match_entry.copy()
        merged_sets = []
        for set_index in range(entry_index + 1, len(self.matches)):
            current_set = self.matches[set_index]
            if current_set is None:
                continue
            if jaccard_index(match_entry, current_set) < self.jaccard_index_threshold:
                continue
            extended_set.update(current_set)
            # Set the index of the newly built merged set as reference
            # self._merge_indications[set_index].append(len(self.merged_matches))
            merged_sets.append(set_index)
        if len(self._merge_indications[entry_index]) == 0:
            self.merged_matches.append(extended_set)
            extended_set_index = len(self.merged_matches) - 1
            self._peptide_mappings[entry_index].append(extended_set_index)
            for merged_set_index in merged_sets:
                self._peptide_mappings[merged_set_index].append(extended_set_index)
        else:
            # TODO
            pass

    # TODO: Add return type
    def merge_matches(self):
        # TODO: Pop entries instead to reduce mem consumption?
        for entry_index, match_entry in enumerate(self.matches):
            self._process_match_entry(match_entry, entry_index)
