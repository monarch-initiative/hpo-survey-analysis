import pytest
from phenom.utils.simulate import rerank_by_average

test_data = [
    (
        [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 6, 6, 7, 8],
        [4, 4, 4, 4, 4, 4, 5, 6, 7, 8, 9, 9, 9, 10]
        #[3, 3, 3, 3, 3, 4, 5, 4, 5, 6, 6, 6, 7, 8]
    ),
    (
        [1, 2, 2, 2, 2, 3, 3, 4, 5, 6, 7, 8, 8, 9, 9, 9, 9],
        [1, 2, 4, 4, 4, 4, 4, 4, 5, 6, 7, 8, 8, 8, 9, 9, 9]
    ),
    (
        [1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 8, 8, 9, 9, 9, 10],
        [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 5, 5, 5]
    ),
]


@pytest.mark.parametrize("input_ranks, expected_ranks", test_data)
def test_reranking(input_ranks, expected_ranks):
    """
    Test reranking by adjusting ties
    function: phenom.utils.simulate.rerank_by_average
    """
    mock_owlsim = [{"rank": rank} for rank in input_ranks]
    rankings = rerank_by_average(mock_owlsim)
    print(rankings)
    assert rankings == expected_ranks
