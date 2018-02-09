#!/usr/bin/env python

import os
import pytest

import numpy as np

from shakemap.utils.generic_amp import get_generic_amp_factors


do_test = True

class Dummy(object):
    pass

def test_generic_amp():
    #
    # Using the LA basin test data set, make a set of lons and lats
    # Bounds of data set: -119.284/-117.284/32.9/34.6
    # But we want to go somewhat beyond the bounds in order to test
    # that functionality
    #
    lons = np.linspace(-121.0, -118.5, 25)
    lons = np.append(lons, np.linspace(-118.5, -116.0, 25))
    lats = np.linspace(33.0, 35.8, 25)
    lats = np.append(lats, np.linspace(33.0, 35.8, 25))

    sx = Dummy()
    sx.lons = lons
    sx.lats = lats

    gaf = get_generic_amp_factors(sx, 'PGA')
    gaf_target = np.array(
      [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
    if do_test is True:
        assert np.allclose(gaf, gaf_target)
    else:
        print('PGA:')
        print(repr(gaf))

    gaf = get_generic_amp_factors(sx, 'SA(0.3)')
    gaf_target = np.array(
      [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
    if do_test is True:
        assert np.allclose(gaf, gaf_target)
    else:
        print('SA(0.3):')
        print(repr(gaf))

    gaf = get_generic_amp_factors(sx, 'SA(2.0)')
    gaf_target = np.array(
      [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
    if do_test is True:
        assert np.allclose(gaf, gaf_target)
    else:
        print('SA(2.0):')
        print(repr(gaf))

    gaf = get_generic_amp_factors(sx, 'SA(4.0)')
    gaf_target = np.array(
      [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
    if do_test is True:
        assert np.allclose(gaf, gaf_target)
    else:
        print('SA(4.0):')
        print(repr(gaf))

    gaf = get_generic_amp_factors(sx, 'SA(0.1)')
    gaf_target = np.array(
      [ 1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
    if do_test is True:
        assert np.allclose(gaf, gaf_target)
    else:
        print('SA(0.1):')
        print(repr(gaf))


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_generic_amp()

