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
        sets_to_merge = []
        for set_index, set_to_compare in enumerate(self.matches):
            if set_to_compare is None:
                continue
            if (
                jaccard_index(match_entry, set_to_compare)
                < self.jaccard_index_threshold
            ):
                continue
            sets_to_merge.append(set_index)

        entry_merge_indications = self._merge_indications.pop(0)
        if len(entry_merge_indications) == 0:
            self.merged_matches.append(match_entry)
            self.peptide_mappings.append([entry_index])
            for set_to_merge in sets_to_merge:
                self._merge_indications[set_to_merge].append(
                    len(self.merged_matches) - 1
                )
        else:
            for merge_indication in entry_merge_indications:
                self.merged_matches[merge_indication].update(match_entry)
                self.peptide_mappings[merge_indication].append(entry_index)
                for set_to_merge in sets_to_merge:
                    self._merge_indications[set_to_merge].append(merge_indication)

    def merge_matches(self) -> Tuple[List[Set[int]], List[List[int]]]:
        entry_index = 0
        while len(self.matches) > 0:
            self._process_match_entry(self.matches.pop(0), entry_index)
            entry_index += 1

        return (self.merged_matches, self.peptide_mappings)
