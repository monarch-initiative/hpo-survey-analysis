from typing import Set, Union, Sequence, List
from functools import reduce
from itertools import chain
import math


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
    column_len = len(matrix)
    row_len = len(matrix[0])  # assumes all rows same length
    flipped_matrix = [[0 for x in range(column_len)] for y in range(row_len)]
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


# Should this be in a separate enrichment module?


def hypergeometric_test(matrix: List[List[int]]) -> float:
    """
    Hypergeometric test
    https://en.wikipedia.org/wiki/Hypergeometric_distribution
    :param matrix: 2x2 matrix
    :return: float, probability value between 0-1
    """
    a = matrix[0][0]
    b = matrix[0][1]
    c = matrix[1][0]
    d = matrix[1][1]
    numerator = math.factorial(a + b) * math.factorial(c + d) \
                * math.factorial(a + c) * math.factorial(b + d)
    denominator = math.factorial(a) * math.factorial(b) \
                  * math.factorial(c) * math.factorial(d) \
                  * math.factorial(a + b + c + d)
    return numerator/denominator


def fisher_exact(matrix, direction="greater"):
    """
    References: https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    This is also a good bio related description:
    http://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/#setup
    :param matrix:
    :param direction: str greater or lesser
    :return:
    """
    p_value =  hypergeometric_test(matrix)
    if direction == "greater":
        while matrix[0][1] > 0 and matrix[1][0] > 0:
            matrix[0][0] += 1
            matrix[0][1] -= 1
            matrix[1][0] -= 1
            matrix[1][1] += 1
            p_value += hypergeometric_test(matrix)

    elif direction == "lesser":
        while matrix[0][0] > 0 and matrix[1][1] > 0:
            matrix[0][0] -= 1
            matrix[0][1] += 1
            matrix[1][0] += 1
            matrix[1][1] -= 1
            p_value += hypergeometric_test(matrix)
    else:
        raise ValueError("only accepts greater or lesser")
    return p_value
