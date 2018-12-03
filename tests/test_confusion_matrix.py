import pytest
import os
import json
from phenom.model.synthetic import SyntheticProfile
from phenom.utils.simulate import rerank_by_average, create_confusion_matrix_per_threshold
from unittest.mock import MagicMock, patch

# Test files and data
mock_owlsim_output = open(os.path.join(os.path.dirname(__file__),
                                       'resources/owlsim3-bayes.json'), 'r')
owlsim_output = json.load(mock_owlsim_output)
mock_owlsim_output.close()

test_rank_data = [
    (
        # Input rankings
        [1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 6, 6, 7, 8],
        # Expected output rankings
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


@pytest.mark.parametrize("input_ranks, expected_ranks", test_rank_data)
def test_reranking(input_ranks, expected_ranks):
    """
    Test reranking by adjusting ties
    function: phenom.utils.simulate.rerank_by_average
    """
    rankings = rerank_by_average(input_ranks)
    assert rankings == expected_ranks


@patch('phenom.utils.simulate.owlsim_classify', MagicMock(return_value=owlsim_output))
def test_create_confusion_by_rank():
    """
    Test output confusion matrix when using rank as threshold
    """
    classes_to_eval = 10
    thresholds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    threshold_type = 'rank'
    fake_profiles = [SyntheticProfile('1', ['2'], "MONDO:0007924")]
    owlsim_url = ''

    expected_table = {
        # true_pos, false_pos, false_neg, true_neg
        1:  [0, 1, 1, 8],
        2:  [0, 1, 1, 8],
        3:  [1, 3, 0, 6],
        4:  [1, 4, 0, 5],
        5:  [1, 5, 0, 4],
        6:  [1, 6, 0, 3],
        7:  [1, 7, 0, 2],
        8:  [1, 8, 0, 1],
        9:  [1, 9, 0, 0],
        10: [1, 9, 0, 0],
        11: [1, 9, 0, 0],
        12: [1, 9, 0, 0]
    }

    confusion_by_rank = create_confusion_matrix_per_threshold(
        fake_profiles, owlsim_url, classes_to_eval, thresholds, threshold_type)

    assert confusion_by_rank == expected_table


@patch('phenom.utils.simulate.owlsim_classify', MagicMock(return_value=owlsim_output))
def test_create_confusion_by_prob():
    """
    Test output confusion matrix when using post probability as threshold
    """
    classes_to_eval = 10
    thresholds = [1.0, .005, 1.e-030, 1.e-060, 1.e-090, 1.e-120, 1.e-150,
                  1.e-180, 1.e-210, 1.e-240, 1.e-270, 1.e-300, 0.0]
    threshold_type = 'probability'
    fake_profiles = [SyntheticProfile('1',['2'],"MONDO:0007924")]
    owlsim_url = ''

    expected_table = {
        # true_pos, false_pos, false_neg, true_neg
        1.0:    [0, 0, 1, 9],
        0.005:  [0, 1, 1, 8],
        1e-30:  [1, 7, 0, 2],
        1e-60:  [1, 7, 0, 2],
        1e-90:  [1, 7, 0, 2],
        1e-120: [1, 8, 0, 1],
        1e-150: [1, 8, 0, 1],
        1e-180: [1, 8, 0, 1],
        1e-210: [1, 8, 0, 1],
        1e-240: [1, 8, 0, 1],
        1e-270: [1, 8, 0, 1],
        1e-300: [1, 8, 0, 1],
        0.0:    [1, 9, 0, 0]
    }

    confusion_by_rank = create_confusion_matrix_per_threshold(
        fake_profiles, owlsim_url, classes_to_eval, thresholds, threshold_type)

    assert confusion_by_rank == expected_table
