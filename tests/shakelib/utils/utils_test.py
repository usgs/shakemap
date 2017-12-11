#!/usr/bin/env python
import io
import os.path
import pytest

from shakelib.rupture.origin import Origin
from shakelib.rupture.factory import get_rupture
from shakelib.utils.utils import get_extent, is_stable

homedir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(homedir, 'utils_data')


def test_get_extent_small_point():
    #
    # Do a point rupture
    # Small magnitude
    #
    eventfile = os.path.join(datadir, 'event_wenchuan_small.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(rupture)
    assert W == 102.33333333333333
    assert E == 104.41666666666667
    assert S == 30.083333333333332
    assert N == 31.883333333333333


def test_get_extent_small_complex():
    #
    # Do a complex rupture
    #
    eventfile = os.path.join(datadir, 'event_wenchuan.xml')
    origin = Origin.fromFile(eventfile)
    faultfile = os.path.join(datadir, 'Hartzell11_fault.txt')
    rupture = get_rupture(origin, faultfile)
    W, E, S, N = get_extent(rupture)
    assert W == 100.03333333333333
    assert E == 108.93333333333334
    assert S == 27.9
    assert N == 35.55


def test_get_extent_bad_usage():
    #
    # Test bad usage
    #
    with pytest.raises(TypeError):
        W, E, S, N = get_extent()
    with pytest.raises(TypeError):
        W, E, S, N = get_extent(34)


def test_get_extent_aspect():
    #
    # Test ruptures with weird aspect ratios
    #
    eventfile = os.path.join(datadir, 'event_wenchuan.xml')
    origin = Origin.fromFile(eventfile)
    #
    # Long horizontal rupture
    #
    rrep = io.StringIO(
            """
            30.0 100.0 0.0
            30.0 110.0 0.0
            30.0 110.0 5.0
            30.0 100.0 5.0
            30.0 100.0 0.0
            """
           )
    rupture = get_rupture(origin, rrep)
    W, E, S, N = get_extent(rupture)
    assert W == 97.18333333333334
    assert E == 113.51666666666667
    assert S == 25.616666666666667
    assert N == 34.06666666666666
    #
    # Long vertical rupture
    #
    rrep = io.StringIO(
            """
            24.0 100.0 0.0
            36.0 100.0 0.0
            36.0 100.0 5.0
            24.0 100.0 5.0
            24.0 100.0 0.0
            """
           )
    rupture = get_rupture(origin, rrep)
    W, E, S, N = get_extent(rupture)
    assert W == 92.53333333333333
    assert E == 108.89999999999999
    assert S == 21.03333333333333
    assert N == 38.46666666666667


def test_get_extent_stable_small():
    #
    # Do an event in a stable region
    # Small magnitude
    #
    eventfile = os.path.join(datadir, 'event_oklahoma.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(rupture)
    assert W == -98.5
    assert E == -96.28333333333333
    assert S == 34.766666666666666
    assert N == 36.583333333333336


def test_get_extent_stable_large():
    #
    # Do an event in a stable region
    # Large magnitude
    #
    eventfile = os.path.join(datadir, 'event_oklahoma_large.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(rupture)
    assert W == -106.14999999999999
    assert E == -86.76666666666667
    assert S == 27.55
    assert N == 43.03333333333333


def test_is_stable():

    assert not is_stable(0, 0)
    assert is_stable(-87.0, 36.0)
    assert not is_stable(-112.0, 39.0)
    assert is_stable(135.0, -30.0)
    assert not is_stable(180.0, 0.0)
    assert is_stable(-60.0, -15.0)


if __name__ == '__main__':
    test_get_extent_small_point()
    test_get_extent_small_complex()
    test_get_extent_bad_usage()
    test_get_extent_aspect()
    test_get_extent_stable_small()
    test_get_extent_stable_large()
    test_is_stable()
