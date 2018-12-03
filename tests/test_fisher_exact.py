import pytest


# Copied from https://github.com/brentp/fishers_exact_test/blob/398c8c/tests.py
# Computed by ``fisher.test`` in R 3.2.2 and printed with
# ``sprintf(".16f")``.

test_data = [
    (
        # Input contingency table
        [[100, 2], [1000, 5]],
        # left_tail, right_tail, two sided
        (0.1300759363430016, 0.9797904453147230, 0.1300759363430016)),
    (
        [[2, 100], [5, 1000]],
        (0.9797904453147230, 0.1300759363430016, 0.1300759363430016)),
    (
        [[2, 7], [8, 2]],
        (0.0185217259520665, 0.9990149169715733, 0.0230141375652212)),
    (
        [[5, 1], [10, 10]],
        (0.9782608695652173, 0.1652173913043478, 0.1973244147157191)),
    (
        [[5, 15], [20, 20]],
        (0.0562577507439996, 0.9849086665340765, 0.0958044001247763)),
    (
        [[5, 16], [20, 25]],
        (0.0891382278309642, 0.9723490195633506, 0.1725864953812995)),
    (
        [[10, 5], [10, 1]],
        (0.1652173913043479, 0.9782608695652174, 0.1973244147157192)),
    (
        [[10, 5], [10, 0]],
        (0.0565217391304348, 1.0000000000000000, 0.0612648221343874)),
    (
        [[5, 0], [1, 4]],
        (1.0000000000000000, 0.0238095238095238, 0.0476190476190476)),
    (
        [[0, 5], [1, 4]],
        (0.5000000000000000, 1.0000000000000000, 1.0000000000000000)),
    (
        [[5, 1], [0, 4]],
        (1.0000000000000000, 0.0238095238095238, 0.0476190476190476)),
    (
        [[0, 1], [3, 2]],
        (0.4999999999999999, 1.0000000000000000, 1.0000000000000000))
]


def compare_with_r(table, expected, fisher_exact):
    epsilon = 1e-10
    _, left_tail = fisher_exact(table, "less")
    _, right_tail = fisher_exact(table, "greater")
    _, two_sided = fisher_exact(table, "two-sided")
    assert abs(left_tail - expected[0]) < epsilon
    assert abs(right_tail - expected[1]) < epsilon
    assert abs(two_sided - expected[2]) < epsilon


def compare_with_scipy_gh_issue(fisher_exact):
    """
    https://github.com/scipy/scipy/issues/4130
    """
    table = [[345, 455], [260, 345]]
    _, two_sided = fisher_exact(table, "two-sided")
    expected = 0.9567
    assert round(two_sided, 4) == expected


@pytest.mark.parametrize("table,expected", test_data)
def test_phenom_against_r(table, expected):
    from phenom.math.enrichment import fisher_exact
    compare_with_r(table, expected, fisher_exact)


@pytest.mark.parametrize("table,expected", test_data)
def test_scipy_against_r(table, expected):
    from scipy.stats import fisher_exact
    compare_with_r(table, expected, fisher_exact)


def test_phenom_against_scipy_issue():
    from phenom.math.enrichment import fisher_exact
    compare_with_scipy_gh_issue(fisher_exact)


@pytest.mark.skip(reason="https://github.com/scipy/scipy/issues/4130")
def test_scipy_against_scipy_issue():
    from scipy.stats import fisher_exact
    compare_with_scipy_gh_issue(fisher_exact)
