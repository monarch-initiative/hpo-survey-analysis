from typing import Union, Sequence, List
from phenom.math import math_utils
from itertools import chain


# Union types
Num = Union[int, float]


def flip_matrix(matrix: Sequence[Sequence]) -> List[List]:
    """
    swap rows and columns in a list of lists
    """
    column_len = len(matrix)
    row_len = len(matrix[0])  # assumes all rows same length
    flipped_matrix = [[0 for x in range(column_len)] for y in range(row_len)]
    for row_index, row in enumerate(matrix):
        for col_index, cell in enumerate(row):
            flipped_matrix[col_index][row_index] = cell
    return flipped_matrix


def sym_bma_score(matrix: Sequence[Sequence[Num]]) -> float:
    """
    symmetric best max average score
    """
    forwards = [max(row) for row in matrix]
    flipped = flip_matrix(matrix)
    backwards = [max(row) for row in flipped]
    combined_maxes = forwards + backwards
    return mean(combined_maxes)


def max_score(matrix: Sequence[Sequence[Num]]) -> float:
    return max(list(chain.from_iterable(matrix)))


def bma_score(matrix: Sequence[Sequence[Num]]) -> float:
    """
    best max average score
    """
    return math_utils.mean([max(row) for row in matrix])

def best_min_avg(matrix: Sequence[Sequence[Num]]) -> float:
    """
    best min average score
    """
    return math_utils.mean([min(row) for row in matrix])


def avg_score(matrix: Sequence[Sequence[Num]]) -> float:
    """
    average of every value in the matrix
    """
    return math_utils.mean(list(chain.from_iterable(matrix)))


def max_percentage_score(
        query_matrix: Sequence[Sequence[Num]],
        optimal_matrix: Sequence[Sequence[float]]) -> float:
    return max_score(query_matrix) / max_score(optimal_matrix)


def bma_percentage_score(
        query_matrix: Sequence[Sequence[Num]],
        optimal_matrix: Sequence[Sequence[Num]]) -> float:
    return sym_bma_score(query_matrix) / sym_bma_score(optimal_matrix)


def avg_percentage_score(
        query_matrix: Sequence[Sequence[Num]],
        optimal_matrix: Sequence[Sequence[Num]]) -> float:
    return avg_score(query_matrix) / avg_score(optimal_matrix)
