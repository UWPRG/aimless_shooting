import os
import os.path as op

import numpy as np
from numpy.testing import assert_almost_equal

from ..core import ShootingPoint


def test_read_cv_values():
    test_file_loc = op.join(op.dirname(op.abspath(__file__)),
                            'test_data', 'COLVAR')
    sp = ShootingPoint()

    sp._read_cv_values(test_file_loc)
    test_values = sp.cv_values
    true_values = np.array([1.000000, 2.000000, 3.000000])
    for test, true in zip(test_values, true_values):
        assert_almost_equal(test, true)
    return
