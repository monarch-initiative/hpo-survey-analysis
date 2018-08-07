from typing import List
import math


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
    :param matrix: 2x2 matrix
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