from typing import List
from copy import deepcopy
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


def fisher_exact(matrix, direction = "two-sided"):
    """
    References: https://en.wikipedia.org/wiki/Fisher%27s_exact_test
    This is also a good bio related description:
    http://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/#setup

    Slower but more accurate than scipy fisher exact,
    passes https://github.com/scipy/scipy/issues/4130

    :param matrix: 2x2 matrix
    :param direction: str two-tailed, greater or less
    :return:
    """
    if direction not in ["two-sided", "greater", "less"]:
        raise ValueError("only accepts two-sided, greater, less")

    # Probability of associations being random
    hyper_prob =  hypergeometric_test(matrix)
    probabilities = [hyper_prob]

    # odds ratio
    try:
        odds_ratio = (matrix[0][0] * matrix[1][1]) / \
                     (matrix[0][1] * matrix[1][0])
    except ZeroDivisionError:
        odds_ratio = 'inf'

    if direction == "two-sided" or direction == "greater":
        matrix_tmp = deepcopy(matrix)
        while matrix_tmp[0][1] > 0 and matrix_tmp[1][0] > 0:
            matrix_tmp[0][0] += 1
            matrix_tmp[0][1] -= 1
            matrix_tmp[1][0] -= 1
            matrix_tmp[1][1] += 1
            probabilities.append(hypergeometric_test(matrix_tmp))

    if direction == "two-sided" or direction == "less":
        matrix_tmp = deepcopy(matrix)
        while matrix_tmp[0][0] > 0 and matrix_tmp[1][1] > 0:
            matrix_tmp[0][0] -= 1
            matrix_tmp[0][1] += 1
            matrix_tmp[1][0] += 1
            matrix_tmp[1][1] -= 1
            probabilities.append(hypergeometric_test(matrix_tmp))

    if direction == "less" or direction == "greater":
        p_value = sum(probabilities)
    else:
        p_value = sum([prob for prob in probabilities if prob <= hyper_prob])

    return odds_ratio, p_value
