#!/usr/bin/env python

#stdlib imports
import os.path
import tempfile
import shutil

#third party imports
from configobj import ConfigObj
import numpy as np

#neic imports
from shakemap.utils.config import get_config_paths
from shakelib.utils.containers import OutputContainer
from shakemap.coremods.contour import _get_default_intervals

def test_intervals():
    fgrid = np.linspace(0.007,11.7,50)
    intervals = _get_default_intervals(fgrid)
    np.testing.assert_almost_equal(intervals[0],0.01)
    np.testing.assert_almost_equal(intervals[-1],10)

    fgrid = np.linspace(0.1,70,50)
    intervals = _get_default_intervals(fgrid)
    np.testing.assert_almost_equal(intervals[0],0.3)
    np.testing.assert_almost_equal(intervals[-1],30)

    fgrid = np.array([3.5,5.6,8.2])
    intervals = _get_default_intervals(fgrid,interval_type='linear')
    np.testing.assert_almost_equal(intervals[0],4.0)
    np.testing.assert_almost_equal(intervals[-1],8.0)


if __name__ == '__main__':
    test_intervals()
