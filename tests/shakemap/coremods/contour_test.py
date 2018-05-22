#!/usr/bin/env python

import numpy as np
from shakelib.plotting.contour import getContourLevels


def test_intervals():
    fgrid = np.linspace(0.007, 11.7, 50)
    intervals = getContourLevels(np.min(fgrid), np.max(fgrid))
    np.testing.assert_almost_equal(intervals[0], 0.01)
    np.testing.assert_almost_equal(intervals[-1], 10)

    fgrid = np.linspace(0.1, 70, 50)
    intervals = getContourLevels(np.min(fgrid), np.max(fgrid))
    np.testing.assert_almost_equal(intervals[0], 0.2)
    np.testing.assert_almost_equal(intervals[-1], 50)

    fgrid = np.array([3.5, 5.6, 8.2])
    intervals = getContourLevels(np.min(fgrid), np.max(fgrid), itype='linear')
    np.testing.assert_almost_equal(intervals[0], 3.5)
    np.testing.assert_almost_equal(intervals[-1], 8.0)


if __name__ == '__main__':
    test_intervals()
