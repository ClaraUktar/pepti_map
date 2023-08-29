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
        self._moved_merge_indications: Dict[int, int] = {}
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
        elif len(entry_merge_indications) == 1:
            merge_indication = entry_merge_indications[0]
            merge_indication = self._moved_merge_indications.get(
                merge_indication, merge_indication
            )
            self.merged_matches[merge_indication].update(match_entry)
            self.peptide_mappings[merge_indication].append(entry_index)
            for set_to_merge in sets_to_merge:
                self._merge_indications[set_to_merge].append(merge_indication)
        else:
            # TODO: refactor with same part above
            self.merged_matches.append(match_entry)
            self.peptide_mappings.append([entry_index])
            for set_to_merge in sets_to_merge:
                self._merge_indications[set_to_merge].append(
                    len(self.merged_matches) - 1
                )
            # TODO: Do these merge indications also need to be looked up
            # in the moved merge indications?
            for merge_indication in entry_merge_indications:
                self._moved_merge_indications[merge_indication] = (
                    len(self.merged_matches) - 1
                )

    def _construct_final_merge_result(self) -> None:
        set_index = 0
        offset = 0
        while set_index < len(self.merged_matches):
            original_set_index = set_index - offset
            if original_set_index not in self._moved_merge_indications:
                set_index += 1
                continue

            set_to_merge = self.merged_matches.pop(set_index)
            peptide_mapping_to_merge = self.peptide_mappings.pop(set_index)
            offset += 1
            merge_into_index = (
                self._moved_merge_indications[original_set_index] - offset
            )
            self.merged_matches[merge_into_index].update(set_to_merge)
            self.peptide_mappings[merge_into_index].extend(peptide_mapping_to_merge)
            self.peptide_mappings[merge_into_index].sort()
            set_index += 1

    def merge_matches(self) -> Tuple[List[Set[int]], List[List[int]]]:
        entry_index = 0
        while len(self.matches) > 0:
            self._process_match_entry(self.matches.pop(0), entry_index)
            entry_index += 1

        self._construct_final_merge_result()

        return (self.merged_matches, self.peptide_mappings)
