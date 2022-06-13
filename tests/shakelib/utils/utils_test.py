#!/usr/bin/env python

# stdlib imports
import io
import os.path

# third party imports
import numpy as np
import pytest
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from impactutils.rupture.factory import get_rupture
from impactutils.rupture.origin import Origin

# local imports
from shakelib.utils.exception import ShakeLibException
from shakelib.utils.utils import get_extent, is_stable, replace_dyfi
from shakelib.station import StationList
from shakelib.virtualipe import VirtualIPE
from shakelib.gmice.wgrw12 import WGRW12
from shakelib.multigmpe import MultiGMPE

homedir = os.path.dirname(os.path.abspath(__file__))
datadir = os.path.join(homedir, "utils_data")


gmpe_cy14 = ChiouYoungs2014()
gmpe = MultiGMPE.__from_list__([gmpe_cy14], [1.0])
gmice = WGRW12()
ipe = VirtualIPE.__fromFuncs__(gmpe, gmice)


def test_replace_dyfi():
    stationfile = os.path.join(datadir, "nepal_dat.xml")
    dyfifile = os.path.join(datadir, "nepal_dyfi_dat.xml")

    original_stations = StationList.loadFromFiles([stationfile])
    dyfi_stations = StationList.loadFromFiles([dyfifile])
    dcursor = dyfi_stations.cursor
    dcursor.execute("SELECT count(*) from station")
    ndyfi1 = dcursor.fetchone()[0]

    original_cursor = original_stations.cursor
    original_cursor.execute("SELECT count(*) from station")

    original_cursor.execute("SELECT count(*) from station WHERE instrumented=0")
    noriginal_mmi = original_cursor.fetchone()[0]

    original_cursor.execute(
        'SELECT count(*) from station WHERE instrumented=0 and network="DYFI"'
    )
    noriginal_dyfi = original_cursor.fetchone()[0]

    noriginal_observed = noriginal_mmi - noriginal_dyfi

    stations = replace_dyfi(stationfile, dyfifile)

    scursor = stations.cursor
    scursor.execute(
        "SELECT count(*) from station where " 'instrumented=0 and network != "DYFI"'
    )
    nobserved = scursor.fetchone()[0]

    assert nobserved == noriginal_observed
    scursor.execute(
        'SELECT count(*) from station where instrumented=0 and network="DYFI"'
    )
    ndyfi = scursor.fetchone()[0]
    assert ndyfi == ndyfi1


def test_get_extent_small_point():
    #
    # Do a point rupture
    # Small magnitude
    #
    config = {"extent": {"mmi": {"threshold": 4.0, "mindist": 100, "maxdist": 1000}}}
    eventfile = os.path.join(datadir, "event_wenchuan_small.xml")
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, 102.21666666666667)
    np.testing.assert_allclose(E, 104.53333333333333)
    np.testing.assert_allclose(S, 29.983333333333334)
    np.testing.assert_allclose(N, 31.983333333333334)


def test_get_extent_small_complex():
    #
    # Do a complex rupture
    #
    config = {"extent": {"mmi": {"threshold": 4.0, "mindist": 100, "maxdist": 1000}}}
    eventfile = os.path.join(datadir, "event_wenchuan.xml")
    origin = Origin.fromFile(eventfile)
    faultfile = os.path.join(datadir, "Hartzell11_fault.txt")
    rupture = get_rupture(origin, faultfile)
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, 96.349999999999994)
    np.testing.assert_allclose(E, 113.61666666666666)
    np.testing.assert_allclose(S, 24.216666666666665)
    np.testing.assert_allclose(N, 38.799999999999997)


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
    config = {"extent": {"mmi": {"threshold": 4.0, "mindist": 100, "maxdist": 1000}}}
    eventfile = os.path.join(datadir, "event_wenchuan.xml")
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
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, 93.700000000000003)
    np.testing.assert_allclose(E, 118.48333333333333)
    np.testing.assert_allclose(S, 20.800000000000001)
    np.testing.assert_allclose(N, 38.233333333333334)
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
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, 89.333333333333329)
    np.testing.assert_allclose(E, 113.68333333333334)
    np.testing.assert_allclose(S, 17.233333333333334)
    np.testing.assert_allclose(N, 41.666666666666664)


def test_get_extent_stable_small():
    #
    # Do an event in a stable region
    # Small magnitude
    #
    config = {"extent": {"mmi": {"threshold": 4.0, "mindist": 100, "maxdist": 1000}}}
    eventfile = os.path.join(datadir, "event_oklahoma.xml")
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, -98.5)
    np.testing.assert_allclose(E, -96.28333333333333)
    np.testing.assert_allclose(S, 34.766666666666666)
    np.testing.assert_allclose(N, 36.583333333333336)


def test_get_extent_stable_large():
    #
    # Do an event in a stable region
    # Large magnitude
    #
    config = {"extent": {"mmi": {"threshold": 4.0, "mindist": 100, "maxdist": 1000}}}
    eventfile = os.path.join(datadir, "event_oklahoma_large.xml")
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)
    W, E, S, N = get_extent(config=config, ipe=ipe, rupture=rupture)
    np.testing.assert_allclose(W, -104.5)
    np.testing.assert_allclose(E, -89.11666666666666)
    np.testing.assert_allclose(S, 29.25)
    np.testing.assert_allclose(N, 41.633333333333333)


def test_is_stable():

    assert is_stable(0, 0)  # Gulf of Guinea
    assert is_stable(-87.0, 36.0)  # Nashville
    assert not is_stable(-112.0, 39.0)  # Wasatch
    assert is_stable(135.0, -30.0)  # South Australia
    assert is_stable(180.0, 0.0)  # South Pacific (middle of ocean)
    assert is_stable(-60.0, -15.0)  # Brazil, near border of Bolivia


def test_extent_config():
    eventfile = os.path.join(datadir, "event_oklahoma_large.xml")
    origin = Origin.fromFile(eventfile)
    rupture = get_rupture(origin)

    config = {
        "extent": {"magnitude_spans": {"span1": [0, 6, 4, 3], "span2": [6, 10, 6, 4]}}
    }
    extent = get_extent(config=config, ipe=ipe, rupture=rupture)
    cmp_extent = (-99.4166667, -95.4083333, 32.675, 38.6833333)
    np.testing.assert_almost_equal(cmp_extent, extent)

    # Test for relative_offset
    config = {
        "extent": {
            "magnitude_spans": {"span1": [0, 6, 4, 3], "span2": [6, 10, 6, 4]},
            "relative_offset": [0.25, 0.5],
        }
    }
    extent = get_extent(config=config, ipe=ipe, rupture=rupture)
    cmp_extent = (-98.4166667, -94.4083333, 35.675, 41.6833333)
    np.testing.assert_almost_equal(cmp_extent, extent)

    config = {"extent": {"bounds": {"extent": [-100, 32, -95, 37]}}}
    extent = get_extent(config=config, ipe=ipe, rupture=rupture)
    cmp_extent = [-100, -95, 32, 37]
    np.testing.assert_almost_equal(extent, cmp_extent)


def test_exception():
    exception = ShakeLibException("Test exception")
    # Check __str__ override is correct
    assert str(exception) == "'Test exception'"


if __name__ == "__main__":
    test_replace_dyfi()
    test_get_extent_small_point()
    test_get_extent_small_complex()
    test_get_extent_bad_usage()
    test_get_extent_aspect()
    test_get_extent_stable_small()
    test_get_extent_stable_large()
    test_is_stable()
    test_extent_config()
    test_exception()
