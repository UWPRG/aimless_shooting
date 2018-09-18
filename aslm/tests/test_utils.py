import os

from ..utils import get_colvar_header


def test_get_colvar_header():
    test_file_loc = os.path.join(os.path.abspath(__file__),
                                 'test_data', 'COLVAR')

    test_header = get_colvar_header(test_file_loc)
    assert test_header == ['d1', 'd2', 'd3']
    return
