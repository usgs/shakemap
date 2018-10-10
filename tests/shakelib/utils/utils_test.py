#!/usr/bin/env python

# stdlib imports
import io
import os.path

# third party imports
import numpy as np
import pytest

# local imports
from shakelib.utils.exception import ShakeLibException
from shakelib.utils.utils import get_extent, is_stable
from shakelib.rupture.factory import get_rupture
from shakelib.rupture.origin import Origin


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
    np.testing.assert_allclose(W, 93.66666666666667)
    np.testing.assert_allclose(E, 115.08333333333333)
    np.testing.assert_allclose(S, 21.55)
    np.testing.assert_allclose(N, 39.55)


def test_get_extent_small_complex():
    #
    # Do a complex rupture
    #
    eventfile = os.path.join(datadir, 'event_wenchuan.xml')
    origin = Origin.fromFile(eventfile)
    faultfile = os.path.join(datadir, 'Hartzell11_fault.txt')
    rupture = get_rupture(origin, faultfile)
    W, E, S, N = get_extent(rupture)
    np.testing.assert_allclose(W, 100.05)
    np.testing.assert_allclose(E, 108.91666666666667)
    np.testing.assert_allclose(S, 27.916666666666668)
    np.testing.assert_allclose(N, 35.53333333333333)


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
            100.0 30.0  0.0
            110.0 30.0 0.0
            110.0 30.0 5.0
            100.0 30.0 5.0
            100.0 30.0 0.0
            """
    )
    rupture = get_rupture(origin, rrep)
    W, E, S, N = get_extent(rupture)
    np.testing.assert_allclose(W, 97.3)
    np.testing.assert_allclose(E, 113.66666666666667)
    np.testing.assert_allclose(S, 24.0)
    np.testing.assert_allclose(N, 35.68333333333333)
    #
    # Long vertical rupture
    #
    rrep = io.StringIO(
        """
            100.0 24.0 0.0
            100.0 36.0 0.0
            100.0 36.0 5.0
            100.0 24.0 5.0
            100.0 24.0 0.0
            """
    )
    rupture = get_rupture(origin, rrep)
    W, E, S, N = get_extent(rupture)
    np.testing.assert_allclose(W, 92.23333333333333)
    np.testing.assert_allclose(E, 109.26666666666667)
    np.testing.assert_allclose(S, 21.016666666666666)
    np.testing.assert_allclose(N, 38.43333333333333)


def test_get_extent_stable_small():
    #
    # Do an event in a stable region
    # Small magnitude
    #
    eventfile = os.path.join(datadir, 'event_oklahoma.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(rupture)
    np.testing.assert_allclose(W, -107.45)
    np.testing.assert_allclose(E, -84.8)
    np.testing.assert_allclose(S, 26.166666666666668)
    np.testing.assert_allclose(N, 44.15)


def test_get_extent_stable_large():
    #
    # Do an event in a stable region
    # Large magnitude
    #
    eventfile = os.path.join(datadir, 'event_oklahoma_large.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(rupture)
    np.testing.assert_allclose(W, -100.61666666666666)
    np.testing.assert_allclose(E, -93.96666666666667)
    np.testing.assert_allclose(S, 32.93333333333333)
    np.testing.assert_allclose(N, 38.333333333333336)


def test_is_stable():

    assert is_stable(0, 0)  # Gulf of Guinea
    assert is_stable(-87.0, 36.0)  # Nashville
    assert not is_stable(-112.0, 39.0)  # Wasatch
    assert is_stable(135.0, -30.0)  # South Australia
    assert is_stable(180.0, 0.0)  # South Pacific (middle of ocean)
    assert is_stable(-60.0, -15.0)  # Brazil, near border of Bolivia


def test_extent_config():
    eventfile = os.path.join(datadir, 'event_oklahoma_large.xml')
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)

    config = {'extent': {
        'magnitude_spans': {
            'span1': [0, 6, 4, 3],
            'span2': [6, 10, 6, 4]
        }
    }}
    extent = get_extent(rupture, config)
    cmp_extent = (-99.416, -95.416, 32.679, 38.679)
    np.testing.assert_almost_equal(cmp_extent, extent)

    config = {
        'extent': {
            'bounds': {
                'extent': [-100, 32, -95, 37]
            }
        }
    }
    extent = get_extent(rupture, config)
    cmp_extent = [-100, -95, 32, 37]
    np.testing.assert_almost_equal(extent, cmp_extent)


def test_exception():
    exception = ShakeLibException('Test exception')
    # Check __str__ override is correct
    assert str(exception) == "'Test exception'"


if __name__ == '__main__':
    test_get_extent_small_point()
    test_get_extent_small_complex()
    test_get_extent_bad_usage()
    test_get_extent_aspect()
    test_get_extent_stable_small()
    test_get_extent_stable_large()
    test_is_stable()
    test_extent_config()
    test_exception()
