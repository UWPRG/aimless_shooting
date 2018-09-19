import os
import os.path as op

from ..utils import log_run, get_colvar_header, log_header


def test_log_run():
    test_name = "test1"
    test_CVs = [1, 2, 3, 4, 5]
    test_result = "accepted"
    log_run(test_name, test_CVs, test_result, filename="test.log")
    with open("test.log", 'r') as f:
        line = f.readline()
    os.remove("test.log")
    line = line.split()
    line[-1] = line[-1].rstrip()
    assert line[0] == "test1"
    for test, true in zip(line[1:6], test_CVs):
        assert float(test) == float(true)
    # assert line[1:6] == test_CVs
    assert test_result == "accepted"
    return


def test_get_colvar_header():
    test_file_loc = op.join(op.dirname(op.abspath(__file__)),
                            'test_data', 'COLVAR')

    test_header = get_colvar_header(test_file_loc)
    assert test_header == ['d1', 'd2', 'd3']
    return


def test_log_header():
    cv_header = ['d1', 'd2', 'd3']
    log_header(cv_header, filename='test.log')

    with open('test.log', 'r') as f:
        line = f.readline()
    os.remove('test.log')

    line = line.rstrip()
    true_line = "NAME d1 d2 d3 RESULT"
    assert line == true_line
    return
