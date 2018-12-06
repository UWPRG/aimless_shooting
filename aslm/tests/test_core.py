import os
import os.path as op

import numpy as np
from numpy.testing import assert_almost_equal

from ..core import AimlessShooting, ShootingPoint, find_and_replace


def test_read_cv_values():
    test_file_loc = op.join(op.dirname(op.abspath(__file__)),
                            'test_data', 'COLVAR2')
    sp = ShootingPoint(name='test', input_file='.')

    sp._read_cv_values(test_file_loc)
    test_values = sp.cv_values
    true_values = np.array([1.000000, 2.000000, 3.000000])
    for test, true in zip(test_values, true_values):
        assert_almost_equal(test, true)
    return


def test_find_and_replace():
    test_line = "abcABC123456"
    test_sub = {"ABC": "EDF", "456": "789"}
    test_result = find_and_replace(test_line, test_sub)
    assert test_result == "abcEDF123789"
    return


def test_check_if_commited():
    sp = ShootingPoint(name='test', input_file='.')
    test_file_loc = op.join(op.dirname(op.abspath(__file__)),
                            'test_data', 'commit.log')
    test_result = sp.check_if_committed(test_file_loc)
    assert test_result == 2, "Not getting correct basin."

    test_file_loc = op.join(op.dirname(op.abspath(__file__)),
                            'test_data', 'no_commit.log')
    test_result = sp.check_if_committed(test_file_loc)
    assert not test_result, "Not reporting 'None' when it does not commit."
    return


def test_log():
    sp = ShootingPoint(name='test')
    sp.name = 'test'
    sp.cv_values = [1, 2, 3]
    sp.result = 'accepted'
    sp.log(filename="test.log")
    with open("test.log", 'r') as f:
        line = f.readline()
    os.remove("test.log")
    line = line.split(' ')
    line[-1] = line[-1].rstrip()
    assert line[0] == "test"
    for test, true in zip(line[1:4], sp.cv_values):
        assert float(test) == float(true)
    # assert line[1:6] == test_CVs
    assert line[-1] == "accepted"


def test_generate_guesses():
    assert 1 == 0, "Write unit test."
    return
