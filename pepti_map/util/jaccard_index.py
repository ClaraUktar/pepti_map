from typing import Set


def jaccard_index(set1: Set, set2: Set) -> float:
    if len(set1) == 0 or len(set2) == 0:
        return 0.0
    return len(set1.intersection(set2)) / len(set1.union(set2))
