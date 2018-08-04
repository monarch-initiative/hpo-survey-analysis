from typing import Set, Union, Sequence, List
from functools import reduce
from itertools import chain


# Union types
Num = Union[int, float]


def jaccard(set1: Set, set2: Set) -> float:
    return len(set1.intersection(set2))/len(set1.union(set2))


def geometric_mean(values: Sequence[Num]) -> float:
    return reduce(lambda x,y: x*y, values)**(1/(len(values)))


def mean(values: Sequence[Num]) -> float:
    return reduce(lambda x,y: x+y, values)/(len(values))


def flip_matrix(matrix: List[List]) -> List[List]:
    """
    swap rows and columns in a list of lists
    """
    flipped_matrix = [[0 for column_size in range(len(matrix))] for row_size in range(len(matrix[0]))]
    for row_index, row in enumerate(matrix):
        for col_index, cell in enumerate(row):
            flipped_matrix[col_index][row_index] = cell
    return flipped_matrix


def max_score(matrix: Sequence[Sequence[float]]) -> float:
    return max(list(chain.from_iterable(matrix)))


def avg_score(matrix: Sequence[Sequence[float]]) -> float:
    return mean([max(row) for row in matrix])


def max_percentage_score(
        query_matrix: Sequence[Sequence[float]],
        optimal_matrix: Sequence[Sequence[float]]) -> float:
    return max_score(query_matrix) / max_score(optimal_matrix)


def avg_percentage_score(
        query_matrix: Sequence[Sequence[float]],
        optimal_matrix: Sequence[Sequence[float]]) -> float:
    return avg_score(query_matrix) / avg_score(optimal_matrix)