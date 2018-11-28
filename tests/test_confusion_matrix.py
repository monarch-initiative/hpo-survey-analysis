import pytest
from phenom.utils.simulate import rerank_by_average

test_data = [
    (
        [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 6, 6, 7, 8],
        [3, 3, 3, 3, 3, 4, 5, 6, 7, 9, 9, 9, 10, 11]
    ),
    (
        [1, 2, 2, 2, 2, 3, 3, 4, 5, 6, 7, 8, 8, 9, 9, 9, 9],
        [1, 4, 4, 4, 4, 6, 6, 7, 8, 9, 10, 12, 12, 14, 14, 14, 14]
    ),
    (
        [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10],
        [1, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 11, 12, 12, 14, 14, 14, 15]
    ),
    (
        [1, 1, 2, 3],
        [2, 2, 3, 4]
    )
]


@pytest.mark.parametrize("input_ranks, expected_ranks", test_data)
def test_reranking(input_ranks, expected_ranks):
    """
    Test reranking by adjusting ties
    function: phenom.utils.simulate.rerank_by_average
    """
    rankings = rerank_by_average(input_ranks)
    assert rankings == expected_ranks
