#!/usr/bin/env python

# stdlib imports
import os
import os.path
import io
import copy
import tempfile
import shutil
import time

# third party
import numpy as np
import pytest
from openquake.hazardlib.geo.geodetic import azimuth
from mapio.geodict import GeoDict
import matplotlib.pyplot as plt
from obspy.core.event import Catalog, FocalMechanism, Event
from obspy.core.event.source import (NodalPlane, NodalPlanes,
                                     Axis, PrincipalAxes)

from shakelib.rupture.origin import Origin, read_moment_quakeml
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.factory import get_rupture, text_to_json
from shakelib.rupture.factory import rupture_from_dict
from shakelib.rupture.factory import validate_json

from shakelib.rupture.utils import (get_local_unit_slip_vector,
                                    get_local_unit_slip_vector_DS,
                                    get_local_unit_slip_vector_SS)
from shakelib.rupture.utils import (get_quad_slip,
                                    get_quad_strike_vector,
                                    get_quad_down_dip_vector)
from shakelib.rupture.utils import rake_to_mech
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture import constants
from impactutils.time.ancient_time import HistoricTime

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?

do_tests = True


def test_no_origin():
    # No argument
    with pytest.raises(Exception):
        get_rupture()

    # Wrong type
    with pytest.raises(Exception):
        get_rupture(7.3)


def test_text_to_json():
    sm_data = '''#Oglesby, D. D., D. S. Dreger, R. A. Harris, \
N. Ratchkovski, and R. Hansen (2004). Inverse kinematic and forward dynamic \
models of the 2002 Denali fault earthquake, Alaska, Bull. Seism. Soc. Am. 94, \
S214-S233.
>
-147.807 63.434        0.000
-147.210 63.472        0.000
-147.267 63.650        22.294
-147.864 63.613        22.294
-147.807 63.434        0.000
>
-146.951 63.551        0.000
-147.551 63.518        0.000
-147.551 63.518        30.000
-146.951 63.551        30.000
-146.951 63.551        0.000
>
-145.968 63.453        0.000
-146.952 63.547        0.000
-146.952 63.547        30.000
-145.968 63.453        30.000
-145.968 63.453        0.000
>
-143.586 62.872        0.000
-145.996 63.427        0.000
-145.996 63.427        30.000
-143.586 62.872        30.000
-143.586 62.872        0.000
>
-142.500 62.114        0.000
-143.669 62.831        0.000
-143.669 62.831        30.000
-142.500 62.114        30.000
-142.500 62.114        0.000'''
    stringio = io.StringIO(sm_data)
    jdict = text_to_json(stringio)
    refcmp = '''Oglesby, D. D., D. S. Dreger, R. A. Harris, \
N. Ratchkovski, and R. Hansen (2004). Inverse kinematic and forward dynamic \
models of the 2002 Denali fault earthquake, Alaska, Bull. Seism. Soc. Am. 94, \
S214-S233.'''
#    assert jdict['metadata']['reference'] == refcmp.replace('\n', ' ')
    assert jdict['metadata']['reference'] == refcmp
    cstart = [[-147.807, 63.434, 0.0],
              [-147.21,  63.472, 0.0],
              [-147.267, 63.65,  22.294],
              [-147.864, 63.613, 22.294],
              [-147.807, 63.434, 0.0]]
    cend = [[-142.5,   62.114, 0.0],
            [-143.669, 62.831, 0.0],
            [-143.669, 62.831, 30.0],
            [-142.5,   62.114, 30.0],
            [-142.5,   62.114, 0.0]]

    assert len(jdict['features'][0]['geometry']['coordinates'][0]) == 5
    assert jdict['features'][0]['geometry']['coordinates'][0][0] == cstart
    assert jdict['features'][0]['geometry']['coordinates'][0][-1] == cend

    sm_data2 = '''#Source: NEIC (2004) based on aftershock distribution
2.98 94.45 10.0
5.0  94.0 10.0
6.0  93.5 10.0
7.0  93.2 10.0
8.0  93.1 10.0
10.0 92.5 10.0
12.0 92.5 10.0
12.0 93.5 40.0
10.0 93.5 40.0
8.0  94.2 40.0
7.0  94.2 40.0
6.0  95.0 40.0
5.0  95.2 40.0
3.30 95.78 40.0
2.98 94.45 10.0'''
    stringio = io.StringIO(sm_data2)
    jdict = text_to_json(stringio, new_format=False)
    coords = [[94.45, 2.98, 10.0],
              [94.0,  5.0,  10.0],
              [93.5,  6.0,  10.0],
              [93.2,  7.0,  10.0],
              [93.1,  8.0,  10.0],
              [92.5, 10.0,  10.0],
              [92.5, 12.0,  10.0],
              [93.5, 12.0,  40.0],
              [93.5, 10.0,  40.0],
              [94.2,  8.0,  40.0],
              [94.2,  7.0,  40.0],
              [95.0,  6.0,  40.0],
              [95.2,  5.0,  40.0],
              [95.78, 3.30, 40.0],
              [94.45, 2.98, 10.0]]
    assert len(jdict['features'][0]['geometry']['coordinates']) == 1
    assert jdict['features'][0]['geometry']['coordinates'][0][0] == coords


def test_rupture_from_dict():

    # Grab an EdgeRupture
    origin = Origin({'id': 'test', 'lat': 0, 'lon': 0, 'depth': 5.0,
                     'mag': 7.0, 'netid': 'us', 'network': '',
                     'locstring': '', 'time':
                     HistoricTime.utcfromtimestamp(time.time())})

    file = os.path.join(homedir, 'rupture_data/cascadia.json')
    rup_original = get_rupture(origin, file)
    d = rup_original._geojson
    rup_from_dict = rupture_from_dict(d)
    assert rup_from_dict._mesh_dx == 0.5

    # Specify mesh_dx
    rup_original = get_rupture(origin, file, mesh_dx=1.0)
    d = rup_original._geojson
    rup_from_dict = rupture_from_dict(d)
    assert rup_from_dict._mesh_dx == 1.0

    # Quad rupture
    file = os.path.join(homedir, 'rupture_data/izmit.json')
    rup_original = get_rupture(origin, file)
    d = rup_original._geojson
    rup_from_dict = rupture_from_dict(d)
    assert rup_from_dict.getArea() == rup_original.getArea()
    # Note, there's a bit of an inconsistency highlighted here because
    # magnitude has key 'magnitude' in the izmit file, but 'mag' in
    # the origin and both get retained.

    # Point rupture from json
    file = os.path.join(homedir, 'rupture_data/point.json')
    rup = get_rupture(origin, file)
    assert rup.lats == 0
    assert rup.lons == 0

    # Point rupture
    origin = Origin({
        'id': 'test',
        'lon': -122.5, 'lat': 37.3,
        'depth': 5.0, 'mag': 7.0, 'netid': 'us',
        'network': '', 'locstring': '',
        'time': HistoricTime.utcfromtimestamp(time.time())
    })

    rup_original = get_rupture(origin)
    d = rup_original._geojson
    rup_from_dict = rupture_from_dict(d)
    assert rup_from_dict.lats == 37.3
    assert rup_from_dict.lons == -122.5

    assert rup_original.getLength() is None
    assert rup_original.getWidth() == constants.DEFAULT_WIDTH
    assert rup_original.getArea() is None
    assert rup_original.getStrike() == constants.DEFAULT_STRIKE
    assert rup_original.getDip() == constants.DEFAULT_DIP
    assert rup_original.getDepthToTop() == constants.DEFAULT_ZTOR
    assert rup_original.getQuadrilaterals() is None
    assert rup_original.depths == 5.0
    # No mech, no tectonic region
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [42.38297204887224, 45.17551530849775])
    else:
        print(rjb[0], rrup[0])
    # Various combinations of mech and tectonic region...
    rup_original._origin._tectonic_region = 'Active Shallow Crust'
    rup_original._origin.mech = 'ALL'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [42.38297204887224, 45.17551530849775])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'RS'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [39.82080242008069, 43.790368434017324])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'NM'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [40.926575292271664, 44.88593633468176])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'SS'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [46.25040107007744, 47.54090263729997])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin._tectonic_region = 'Stable Shallow Crust'
    rup_original._origin.mech = 'ALL'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [43.563645263662636, 46.42217285643245])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'RS'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [42.58702826839786, 46.39811791109895])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'NM'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [43.31021890717387, 47.15925487179388])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin.mech = 'SS'
    rjb, _ = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                     np.array([0.0]))
    rrup, _ = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [47.40232617967255, 49.45750891908971])
    else:
        print(rjb[0], rrup[0])
    rup_original._origin._tectonic_region = 'Somewhere Else'
    rup_original._origin.mech = 'ALL'
    rjb, var = rup_original.computeRjb(np.array([-122.0]), np.array([37.0]),
                                       np.array([0.0]))
    rrup, var = rup_original.computeRrup(np.array([-122.0]), np.array([37.0]),
                                         np.array([0.0]))
    if do_tests is True:
        np.testing.assert_allclose([rjb[0], rrup[0]],
                                   [42.38297204887224, 45.17551530849775])
    else:
        print(rjb[0], rrup[0])

    # This is just zeroes now, so there's not much to check
    gc2 = rup_original.computeGC2(np.array([-122.0]), np.array([37.0]),
                                  np.array([0.0]))
    assert gc2['rx'][0] == 0


def test_EdgeRupture():

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})

    file = os.path.join(homedir, 'rupture_data/cascadia.json')
    rup = get_rupture(origin, file)
    np.testing.assert_allclose(rup.getArea(), 105635.92827547337)

    # Force read Northridge as EdgeRupture
    file = os.path.join(homedir, 'rupture_data/northridge_fault.txt')
    d = text_to_json(file, new_format=True)
    rupt = EdgeRupture(d, origin)
    strike = rupt.getStrike()
    np.testing.assert_allclose(strike, 121.97, atol=0.01)
    dip = rupt.getDip()
    np.testing.assert_allclose(dip, 40.12, atol=0.01)
    L = rupt.getLength()
    np.testing.assert_allclose(L, 17.99, atol=0.01)
    W = rupt.getWidth()
    np.testing.assert_allclose(W, 23.92, atol=0.01)
    ztor = rupt.getDepthToTop()
    np.testing.assert_allclose(ztor, 5, atol=0.01)

    # And again for the same vertices but reversed order
    file = os.path.join(homedir, 'rupture_data/northridge_fixed_fault.txt')
    d = text_to_json(file, new_format=True)
    rupt = EdgeRupture(d, origin)
    strike = rupt.getStrike()
    np.testing.assert_allclose(strike, 121.97, atol=0.01)
    dip = rupt.getDip()
    np.testing.assert_allclose(dip, 40.12, atol=0.01)
    L = rupt.getLength()
    np.testing.assert_allclose(L, 17.99, atol=0.01)
    W = rupt.getWidth()
    np.testing.assert_allclose(W, 23.92, atol=0.01)
    ztor = rupt.getDepthToTop()
    np.testing.assert_allclose(ztor, 5, atol=0.01)

    # Test for fromArrays method
    toplats = np.array([37.0, 38.0])
    toplons = np.array([-120.0, -120.0])
    topdeps = np.array([0.0, 0.0])
    botlats = copy.copy(toplats)
    botlons = copy.copy(toplons)
    botdeps = np.array([10.0, 10.0])
    erup = EdgeRupture.fromArrays(toplons, toplats, topdeps, botlons, botlats,
                                  botdeps, origin)
    # Error: array lengths differ
    with pytest.raises(ShakeLibException) as e:
        qrup = QuadRupture.fromVertices(
            [toplons[0]], [toplats[0]], [topdeps[0]],
            [toplons[1]], [toplats[1]], [topdeps[1]],
            [botlons[1]], [botlats[1]], [botdeps[1]],
            [botlons[0]], [botlats[0]], [botdeps[0]][:-1],
            origin)
    print(str(e))

    # Error: group index too long
    with pytest.raises(ShakeLibException) as e:
        qrup = QuadRupture.fromVertices(
            [toplons[0]], [toplats[0]], [topdeps[0]],
            [toplons[1]], [toplats[1]], [topdeps[1]],
            [botlons[1]], [botlats[1]], [botdeps[1]],
            [botlons[0]], [botlats[0]], [botdeps[0]],
            origin, group_index=[0, 0, 0, 0, 0, 0])
    print(str(e))

    qrup = QuadRupture.fromVertices(
        [toplons[0]], [toplats[0]], [topdeps[0]],
        [toplons[1]], [toplats[1]], [topdeps[1]],
        [botlons[1]], [botlats[1]], [botdeps[1]],
        [botlons[0]], [botlats[0]], [botdeps[0]],
        origin)
    np.testing.assert_allclose(erup.getArea(), 1108.9414759967776)
    np.testing.assert_allclose(erup.getDepthToTop(), 0)
    np.testing.assert_allclose(erup.getLength(), 111.19492664455889)
    np.testing.assert_allclose(
        erup.lats, np.array([37.,  38.,  38.,  37.,  37.,  np.nan]))
    np.testing.assert_allclose(
        erup.lons, np.array([-120., -120., -120., -120., -120.,  np.nan]))
    np.testing.assert_allclose(
        erup.depths, np.array([0.,   0.,  10.,  10.,   0.,  np.nan]))
    np.testing.assert_allclose(
        erup._getGroupIndex(), np.array([0.,   0.]))
    quads = erup.getQuadrilaterals()
    np.testing.assert_allclose(quads[0][0].x, -120.0)

    # Need to also test the distances with EdgeRupture
    lons = np.linspace(-120.1, -121.0, 10)
    lats = np.linspace(37.0, 38, 10)
    deps = np.zeros_like(lons)
    rrup1, _ = qrup.computeRrup(lons, lats, deps)
    rrup2, _ = erup.computeRrup(lons, lats, deps)
    np.testing.assert_allclose(rrup1, rrup2, atol=2e-2)
    rjb1, _ = qrup.computeRjb(lons, lats, deps)
    rjb2, _ = erup.computeRjb(lons, lats, deps)
    np.testing.assert_allclose(rjb1, rjb2, atol=2e-2)
    gc2 = erup.computeGC2(lons, lats, deps)
    targetRy0 = np.array(
        [0., 0.,  0., 0.,  0.,
         0., 0.,  0., 0.,  0.67335931])
    targetRx = np.array(
        [-8.88024949, -17.73390996, -26.56167797, -35.3634266,
         -44.13902929, -52.88835984, -61.61129242, -70.30770154,
         -78.97746209, -87.6204493])
    np.testing.assert_allclose(gc2['ry0'], targetRy0)
    np.testing.assert_allclose(gc2['rx'], targetRx)


def test_QuadRupture():

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})

    # First with json file
    file = os.path.join(homedir, 'rupture_data/izmit.json')
    rupj = get_rupture(origin, file)
    # Then with text file:
    file = os.path.join(homedir, 'rupture_data/Barkaetal02_fault.txt')
    rupt = get_rupture(origin, file)

    np.testing.assert_allclose(rupj.lats, rupt.lats, atol=1e-5)
    np.testing.assert_allclose(rupj.lons, rupt.lons, atol=1e-5)
    np.testing.assert_allclose(rupj._depth, rupt._depth, atol=1e-5)
    np.testing.assert_allclose(rupt.getArea(), 2391.2822653900268, atol=1e-5)

    target = np.array(
        [29.51528,  29.3376,  29.3376,  29.51528005,
         29.51528,       np.nan,  29.87519,  29.61152,
         29.61152,  29.87519021,  29.87519,       np.nan,
         30.11126,  29.88662,  30.11126,  30.11126,
         29.88662,  30.11126,  30.11126,       np.nan,
         30.4654,  30.30494,  30.4654,  30.4654,
         30.30494,  30.4654,  30.4654,       np.nan,
         30.63731,  30.57658,  30.57658,  30.63731011,
         30.63731,       np.nan,  30.93655,  30.729,
         30.729,  30.93655103,  30.93655,       np.nan,
         31.01799,  30.94688,  30.94688,  31.0179905,
         31.01799,       np.nan]
    )
    np.testing.assert_allclose(rupj.lons, target, atol=1e-5)
    target = np.array(
        [40.72733,  40.70985,  40.71185,  40.72932969,
         40.72733,       np.nan,  40.74903,  40.70513,
         40.70713,  40.75102924,  40.74903,       np.nan,
         40.72336,  40.72582,  40.72336,  40.72536,
         40.72782,  40.72536004,  40.72336,       np.nan,
         40.71081,  40.7121,  40.71081,  40.71281,
         40.7141,  40.71281002,  40.71081,       np.nan,
         40.70068,  40.71621,  40.71821,  40.70268025,
         40.70068,       np.nan,  40.79654,  40.69947,
         40.70147,  40.79853872,  40.79654,       np.nan,
         40.84501,  40.80199,  40.80399,  40.84700952,
         40.84501,       np.nan]
    )
    np.testing.assert_allclose(rupj.lats, target, atol=1e-5)
    target = np.array(
        [-0.00000000e+00,  -0.00000000e+00,   2.00000000e+01,
         1.99999325e+01,  -0.00000000e+00,           np.nan,
         -9.31322575e-13,  -0.00000000e+00,   2.00000000e+01,
         1.99998304e+01,  -9.31322575e-13,           np.nan,
         9.31322575e-13,  -0.00000000e+00,   9.31322575e-13,
         2.00000000e+01,   2.00000000e+01,   2.00000095e+01,
         9.31322575e-13,           np.nan,  -0.00000000e+00,
         -0.00000000e+00,  -0.00000000e+00,   2.00000000e+01,
         2.00000000e+01,   2.00000050e+01,  -0.00000000e+00,
         np.nan,  -0.00000000e+00,  -0.00000000e+00,
         2.00000000e+01,   2.00000600e+01,  -0.00000000e+00,
         np.nan,  -0.00000000e+00,  -0.00000000e+00,
         2.00000000e+01,   1.99996249e+01,  -0.00000000e+00,
         np.nan,  -0.00000000e+00,  -0.00000000e+00,
         2.00000000e+01,   1.99998338e+01,  -0.00000000e+00,
         np.nan])
    np.testing.assert_allclose(rupj.depths, target, atol=1e-5)


def test_rupture_depth(interactive=False):
    DIP = 17.0
    WIDTH = 20.0
    GRIDRES = 0.1

    names = ['single', 'double', 'triple',
             'concave', 'concave_simple', 'ANrvSA']
    means = [3.1554422780092461, 2.9224454569459781,
             3.0381968625073563, 2.0522694624400271,
             2.4805390352818755, 2.8740121776209673]
    stds = [2.1895293825074575, 2.0506459673526174,
            2.0244588429154402, 2.0112565876976416,
            2.1599789955270019, 1.6156220309120068]
    xp0list = [np.array([118.3]),
               np.array([10.1, 10.1]),
               np.array([10.1, 10.1, 10.3]),
               np.array([10.9, 10.5, 10.9]),
               np.array([10.9, 10.6]),
               np.array([-76.483, -76.626, -76.757, -76.99, -77.024, -76.925,
                         -76.65, -76.321, -75.997, -75.958])]
    xp1list = [np.array([118.3]),
               np.array([10.1, 10.3]),
               np.array([10.1, 10.3, 10.1]),
               np.array([10.5, 10.9, 11.3]),
               np.array([10.6, 10.9]),
               np.array([-76.626, -76.757, -76.99, -77.024, -76.925, -76.65,
                         -76.321, -75.997, -75.958, -76.006])]
    yp0list = [np.array([34.2]),
               np.array([34.2, 34.5]),
               np.array([34.2, 34.5, 34.8]),
               np.array([34.2, 34.5, 34.8]),
               np.array([35.1, 35.2]),
               np.array([-52.068, -51.377, -50.729, -49.845, -49.192, -48.507,
                         -47.875, -47.478, -47.08, -46.422])]
    yp1list = [np.array([34.5]),
               np.array([34.5, 34.8]),
               np.array([34.5, 34.8, 35.1]),
               np.array([34.5, 34.8, 34.6]),
               np.array([35.2, 35.4]),
               np.array([-51.377, -50.729, -49.845, -49.192, -48.507, -47.875,
                         -47.478, -47.08, -46.422, -45.659])]

    for i in range(0, len(xp0list)):
        xp0 = xp0list[i]
        xp1 = xp1list[i]
        yp0 = yp0list[i]
        yp1 = yp1list[i]
        name = names[i]
        mean_value = means[i]
        std_value = stds[i]

        zp = np.zeros(xp0.shape)
        strike = azimuth(xp0[0], yp0[0], xp1[-1], yp1[-1])
        widths = np.ones(xp0.shape) * WIDTH
        dips = np.ones(xp0.shape) * DIP
        strike = [strike]

    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})

    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths, dips, origin, strike=strike)

    # make a grid of points over both quads, ask for depths
    ymin = np.nanmin(rupture.lats)
    ymax = np.nanmax(rupture.lats)
    xmin = np.nanmin(rupture.lons)
    xmax = np.nanmax(rupture.lons)

    xmin = np.floor(xmin * (1 / GRIDRES)) / (1 / GRIDRES)
    xmax = np.ceil(xmax * (1 / GRIDRES)) / (1 / GRIDRES)
    ymin = np.floor(ymin * (1 / GRIDRES)) / (1 / GRIDRES)
    ymax = np.ceil(ymax * (1 / GRIDRES)) / (1 / GRIDRES)
    geodict = GeoDict.createDictFromBox(
        xmin, xmax, ymin, ymax, GRIDRES, GRIDRES)
    nx = geodict.nx
    ny = geodict.ny
    depths = np.zeros((ny, nx))
    for row in range(0, ny):
        for col in range(0, nx):
            lat, lon = geodict.getLatLon(row, col)
            depth = rupture.getDepthAtPoint(lat, lon)
            depths[row, col] = depth

    np.testing.assert_almost_equal(np.nanmean(depths), mean_value)
    np.testing.assert_almost_equal(np.nanstd(depths), std_value)

    if interactive:
        fig, axes = plt.subplots(nrows=2, ncols=1)
        ax1, ax2 = axes
        xdata = np.append(xp0, xp1[-1])
        ydata = np.append(yp0, yp1[-1])
        plt.sca(ax1)
        plt.plot(xdata, ydata, 'b')
        plt.sca(ax2)
        im = plt.imshow(depths, cmap='viridis_r')  # noqa
        ch = plt.colorbar()  # noqa
        fname = os.path.join(os.path.expanduser('~'),
                             'quad_%s_test.png' % name)
        print('Saving image for %s quad test... %s' % (name, fname))
        plt.savefig(fname)
        plt.close()


def test_slip():

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})

# Make a rupture
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])
    rup = QuadRupture.fromTrace(lon0, lat0, lon1, lat1, z, W, dip, origin)

    slp = get_quad_slip(rup.getQuadrilaterals()[0], 30).getArray()
    slpd = np.array([0.80816457,  0.25350787,  0.53160491])
    np.testing.assert_allclose(slp, slpd)

    slp = get_quad_strike_vector(rup.getQuadrilaterals()[0]).getArray()
    slpd = np.array([0.58311969, 0.27569625, 0.76417472])
    np.testing.assert_allclose(slp, slpd)

    slp = get_quad_down_dip_vector(rup.getQuadrilaterals()[0]).getArray()
    slpd = np.array([0.81219873, -0.17763484, -0.55567895])
    np.testing.assert_allclose(slp, slpd)

    slp = get_local_unit_slip_vector(22, 30, 86).getArray()
    slpd = np.array([0.82714003,  0.38830563,  0.49878203])
    np.testing.assert_allclose(slp, slpd)

    slp = get_local_unit_slip_vector_DS(22, 30, -86).getArray()
    slpd = np.array([-0.80100879, -0.32362856, -0.49878203])
    np.testing.assert_allclose(slp, slpd)

    slp = get_local_unit_slip_vector_SS(22, 80, 5).getArray()
    slpd = np.array([0.3731811, 0.92365564, 0.])
    np.testing.assert_allclose(slp, slpd)

    mech = rake_to_mech(-160)
    assert mech == 'SS'
    mech = rake_to_mech(0)
    assert mech == 'SS'
    mech = rake_to_mech(160)
    assert mech == 'SS'
    mech = rake_to_mech(-80)
    assert mech == 'NM'
    mech = rake_to_mech(80)
    assert mech == 'RS'


def test_northridge():
    rupture_text = """# Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut \
(1996). The Slip History of the 1994 Northridge, California, Earthquake \
Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. \
Seism. Soc. Am. 86, S49-S70.
    -118.421 34.315  5.000
    -118.587 34.401  5.000
    -118.693 34.261 20.427
    -118.527 34.175 20.427
    -118.421 34.315 5.000
    """  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})

    cbuf = io.StringIO(rupture_text)
    rupture = get_rupture(origin, cbuf)
    strike = rupture.getStrike()
    np.testing.assert_allclose(strike, 122.06, atol=0.01)
    dip = rupture.getDip()
    np.testing.assert_allclose(dip, 40.21, atol=0.01)
    L = rupture.getLength()
    np.testing.assert_allclose(L, 17.99, atol=0.01)
    W = rupture.getWidth()
    np.testing.assert_allclose(W, 23.94, atol=0.01)
    nq = rupture.getNumQuads()
    np.testing.assert_allclose(nq, 1)
    ng = rupture.getNumGroups()
    np.testing.assert_allclose(ng, 1)
    nd = rupture.getDeps()
    np.testing.assert_allclose(nd, [5.0, 5.0, 20.427, 20.427, np.nan])
    sind = rupture._getGroupIndex()
    np.testing.assert_allclose(sind, [0])
    ztor = rupture.getDepthToTop()
    np.testing.assert_allclose(ztor, 5, atol=0.01)
    itl = rupture.getIndividualTopLengths()
    np.testing.assert_allclose(itl, 17.99, atol=0.01)
    iw = rupture.getIndividualWidths()
    np.testing.assert_allclose(iw, 23.94, atol=0.01)
    lats = rupture.lats
    lats_d = np.array([34.401, 34.315, 34.175, 34.261, 34.401, np.nan])
    np.testing.assert_allclose(lats, lats_d, atol=0.01)
    lons = rupture.lons
    lons_d = np.array(
        [-118.587, -118.421, -118.527, -118.693, -118.587, np.nan])
    np.testing.assert_allclose(lons, lons_d, atol=0.01)
    ln, lt, de = rupture.getRuptureAsArrays()
    np.testing.assert_allclose(ln, np.array([-118.421, -118.587, -118.693,
                                             -118.527, np.nan]), atol=0.01)
    np.testing.assert_allclose(lt, np.array([34.315, 34.401, 34.261,
                                             34.175, np.nan]), atol=0.01)
    np.testing.assert_allclose(de, [5.0, 5.0, 20.427, 20.427, np.nan])
    mesh = rupture.getRuptureAsMesh()
    np.testing.assert_allclose(mesh.lons, [-118.421, -118.587, -118.693,
                                           -118.527, np.nan])
    np.testing.assert_allclose(mesh.lats, [34.315, 34.401, 34.261,
                                           34.175, np.nan])
    np.testing.assert_allclose(mesh.depths, [5., 5., 20.427, 20.427, np.nan])


def test_parse_complicated_rupture():
    rupture_text = """# SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, \
Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, \
T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, \
A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, \
and W. Page (2002). The Surface Rupture and Slip Distribution of the \
17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. \
Seism. Soc. Am. 92, 43-60.
    29.33760 40.70985 0
    29.51528 40.72733 0
    29.51528 40.72933 20
    29.33760 40.71185 20
    29.33760 40.70985 0
    >
    29.61152 40.70513 0
    29.87519 40.74903 0
    29.87519 40.75103 20
    29.61152 40.70713 20
    29.61152 40.70513 0
    >
    29.88662 40.72582 0
    30.11126 40.72336 0
    30.19265 40.73432 0
    30.19265 40.73632 20
    30.11126 40.72536 20
    29.88662 40.72782 20
    29.88662 40.72582 0
    >
    30.30494 40.71210 0
    30.46540 40.71081 0
    30.56511 40.70739 0
    30.56511 40.70939 20
    30.46540 40.71281 20
    30.30494 40.71410 20
    30.30494 40.71210 0
    >
    30.57658 40.71621 0
    30.63731 40.70068 0
    30.63731 40.70268 20
    30.57658 40.71821 20
    30.57658 40.71621 0
    >
    30.72900 40.69947 0
    30.93655 40.79654 0
    30.93655 40.79854 20
    30.72900 40.70147 20
    30.72900 40.69947 0
    >
    30.94688 40.80199 0
    31.01799 40.84501 0
    31.01799 40.84701 20
    30.94688 40.80399 20
    30.94688 40.80199 0"""  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})
    cbuf = io.StringIO(rupture_text)
    rupture = get_rupture(origin, cbuf)
    strike = rupture.getStrike()
    np.testing.assert_allclose(strike, -100.46, atol=0.01)
    dip = rupture.getDip()
    np.testing.assert_allclose(dip, 89.40, atol=0.01)
    L = rupture.getLength()
    np.testing.assert_allclose(L, 119.56, atol=0.01)
    W = rupture.getWidth()
    np.testing.assert_allclose(W, 20.0, atol=0.01)
    nq = rupture.getNumQuads()
    np.testing.assert_allclose(nq, 9)
    ng = rupture.getNumGroups()
    np.testing.assert_allclose(ng, 7)
    sind = rupture._getGroupIndex()
    np.testing.assert_allclose(sind, [0, 1, 2, 2, 3, 3, 4, 5, 6])
    ztor = rupture.getDepthToTop()
    np.testing.assert_allclose(ztor, 0, atol=0.01)
    itl = rupture.getIndividualTopLengths()
    itl_d = np.array([15.13750778,  22.80237887,  18.98053425,   6.98263853,
                      13.55978731,   8.43444811,   5.41399812,  20.57788056,
                      7.66869463])
    np.testing.assert_allclose(itl, itl_d, atol=0.01)
    iw = rupture.getIndividualWidths()
    iw_d = np.array([20.00122876,  20.00122608,  20.00120173,  20.00121028,
                     20.00121513,  20.00121568,  20.00107293,  20.00105498,
                     20.00083348])
    np.testing.assert_allclose(iw, iw_d, atol=0.01)
    lats = rupture.lats
    lats_d = np.array([40.72733,  40.70985,  40.71185,  40.72932969,
                       40.72733,          np.nan,  40.74903,  40.70513,
                       40.70713,  40.75102924,  40.74903,          np.nan,
                       40.72336,  40.72582,  40.72336,  40.72536,
                       40.72782,  40.72536004,  40.72336,          np.nan,
                       40.71081,  40.7121,  40.71081,  40.71281,
                       40.7141,  40.71281002,  40.71081,          np.nan,
                       40.70068,  40.71621,  40.71821,  40.70268025,
                       40.70068,          np.nan,  40.79654,  40.69947,
                       40.70147,  40.79853872,  40.79654,          np.nan,
                       40.84501,  40.80199,  40.80399,  40.84700952,
                       40.84501,          np.nan])
    np.testing.assert_allclose(lats, lats_d, atol=0.001)
    lons = rupture.lons
    lons_d = np.array([29.51528,  29.3376,  29.3376,  29.51528005,
                       29.51528,          np.nan,  29.87519,  29.61152,
                       29.61152,  29.87519021,  29.87519,          np.nan,
                       30.11126,  29.88662,  30.11126,  30.11126,
                       29.88662,  30.11126,  30.11126,          np.nan,
                       30.4654,  30.30494,  30.4654,  30.4654,
                       30.30494,  30.4654,  30.4654,          np.nan,
                       30.63731,  30.57658,  30.57658,  30.63731011,
                       30.63731,          np.nan,  30.93655,  30.729,
                       30.729,  30.93655103,  30.93655,          np.nan,
                       31.01799,  30.94688,  30.94688,  31.0179905,
                       31.01799,          np.nan])

    np.testing.assert_allclose(lons, lons_d, atol=0.001)


def test_incorrect():
    # Number of points in polyon is even
    rupture_text = """# Source: Ji, C., D. V. Helmberger, D. J. Wald, and \
K.-F. Ma (2003). Slip history and dynamic implications of the 1999 Chi-Chi, \
Taiwan, earthquake, J. Geophys. Res. 108, 2412, doi:10.1029/2002JB001764.
    120.72300 24.27980 	0
    121.00000 24.05000	17
    121.09300 24.07190	17
    121.04300 24.33120	17
    121.04300 24.33120	17
    120.72300 24.27980	0
    >
    120.72300 24.27980	0
    120.68000 23.70000	0
    120.97200 23.60400	17
    121.00000 24.05000	17
    120.72300 24.27980	0
    >
    120.97200 23.60400	17
    120.68000 23.70000	0
    120.58600 23.58850	0
    120.78900 23.40240	17
    120.97200 23.60400	17"""  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf)

    # Top points must be first
    rupture_text = """# Test
    120.72300 24.27980 	0
    121.00000 24.05000	17
    121.09300 24.07190	17
    121.04300 24.33120	17
    120.72300 24.27980	0"""  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf)

    # Wrong order of lat/lon
    rupture_text = """# Test
    -118.421 34.315  5.000
    -118.587 34.401  5.000
    -118.693 34.261 20.427
    -118.527 34.175 20.427
    -118.421 34.315 5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=False)

    # Wrong order of lat/lon
    rupture_text = """# Test
    34.315 -118.421  5.000
    34.401 -118.587  5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421  5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=True)

    # Unclosed segments
    rupture_text = """# Test
    34.315 -118.421  5.000
    34.401 -118.587  5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.6    5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=False)

    # incorrect delimiter
    rupture_text = """#Test
    34.315;-118.421;5.000
    34.401;-118.587;5.000
    34.261;-118.693;20.427
    34.175;-118.527;20.427
    34.315;-118.421;5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=False)

    # incorrect delimiter, new format
    rupture_text = """#Test
    34.315;-118.421;5.000
    34.401;-118.587;5.000
    34.261;-118.693;20.427
    34.175;-118.527;20.427
    34.315;-118.421;5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=True)

    # Not 3 columns
    rupture_text = """#Test
    34.315 -118.421;5.000
    34.401 -118.587;5.000
    34.261 -118.693;20.427
    34.175 -118.527;20.427
    34.315 -118.421;5.000
    """  # noqa
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf, new_format=False)

    # Json incorrect
    test = {
        "metadata": {
            "id": "test",
            "mag": 7.0,
            "lon": 0,
            "mech": "ALL",
            "depth": 5.0,
            "time": "2018-07-02T22:50:03Z",
            "netid": "us",
            "rake": 0.0,
            "lat": 0,
            "network": "",
            "locstring": "",
            "reference": "Test"
        },
        "features": [{
            "type": "Feature",
            "geometry": {
                "coordinates": [[
                    [[-118.421, 34.315, 5.0],
                     [-118.587, 34.401, 5.0],
                     [-118.693, 34.261, 20.427],
                     [-118.527, 34.175, 20.427],
                     [-118.421, 34.315, 5.0]]]],
                "type": "MultiPolygon"
            },
            "properties":{
                "rupture type": "rupture extent"
            }
        }],
        "type": "FeatureCollection"
    }

    # incorrect type
    test_incorrect = copy.deepcopy(test)
    test_incorrect['type'] = 'Feature'
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))

    # Incorrect number of features
    test_incorrect = copy.deepcopy(test)
    test_incorrect['features'].append(['wrong'])
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))

    # no reference
    test_incorrect = copy.deepcopy(test)
    test_incorrect['metadata'].pop('reference', None)
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))

    # incorrect feature type
    test_incorrect = copy.deepcopy(test)
    test_incorrect['features'][0]['type'] = 'fred'
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))

    # incorrect feature geometry type
    test_incorrect = copy.deepcopy(test)
    test_incorrect['features'][0]['geometry']['type'] = 'fred'
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))

    # no coordinates
    test_incorrect = copy.deepcopy(test)
    test_incorrect['features'][0]['geometry'].pop('coordinates', None)
    with pytest.raises(Exception) as e:
        validate_json(test_incorrect)
    print(str(e))


def test_fromTrace():
    xp0 = [0.0]
    xp1 = [0.0]
    yp0 = [0.0]
    yp1 = [0.05]
    zp = [0.0]
    widths = [10.0]
    dips = [45.0]

    # Rupture requires an origin even when not used:
    origin = Origin({
        'id': 'test',
        'lon': -121.81529, 'lat': 37.73707,
        'depth': 5.0, 'mag': 7.0, 'netid': 'us',
        'network': '', 'locstring': '',
        'time': HistoricTime.utcfromtimestamp(time.time())
    })

    # Error: unequal array lengths
    with pytest.raises(ShakeLibException) as e:
        rupture = QuadRupture.fromTrace(
            xp0, yp0, xp1, yp1, zp[:-1], widths,
            dips, origin,
            reference='From J Smith, (personal communication)')
    print(str(e))

    # Error: invalid strike
    with pytest.raises(ShakeLibException) as e:
        rupture = QuadRupture.fromTrace(
            xp0, yp0, xp1, yp1, zp, widths,
            dips, origin, strike=[236.0, 250.0],
            reference='From J Smith, (personal communication)')
    print(str(e))

    # TODO: These write tests exercise code, but we don't check the results
    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths,
        dips, origin,
        reference='From J Smith, (personal communication)')
    fstr = io.StringIO()
    rupture.writeTextFile(fstr)

    tfile = tempfile.NamedTemporaryFile()
    tname = tfile.name
    tfile.close()
    rupture.writeTextFile(tname)
    os.remove(tname)

    tfile = tempfile.NamedTemporaryFile()
    tname = tfile.name
    tfile.close()
    rupture.writeGeoJson(tname)
    os.remove(tname)

    xp0 = [-121.81529, -121.82298]
    xp1 = [-121.82298, -121.83068]
    yp0 = [37.73707, 37.74233]
    yp1 = [37.74233, 37.74758]
    zp = [10, 15]
    widths = [15.0, 20.0]
    dips = [30.0, 45.0]
    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths,
        dips, origin,
        reference='From J Smith, (personal communication)')

    assert rupture.getReference() == 'From J Smith, (personal communication)'
    rorigin = rupture.getOrigin()
    assert rorigin.id == origin.id
    assert rorigin.mag == origin.mag
    assert rorigin.depth == origin.depth

    rx = rupture.getRuptureContext([])
    np.testing.assert_allclose([rx.strike, rx.dip, rx.ztor, rx.width],
                               [-49.183708644954905, 37.638322472702534,
                                9.999999999371358, 17.47024205615428])

    rhyp = rupture.computeRhyp(np.array([-121.5]), np.array([37.0]),
                               np.array([0]))
    repi = rupture.computeRepi(np.array([-121.5]), np.array([37.0]),
                               np.array([0]))
    np.testing.assert_allclose([rhyp[0], repi[0]], [86.709236, 86.564956])


def test_with_quakeml():
    np1 = NodalPlane(strike=259, dip=74, rake=10)
    np2 = NodalPlane(strike=166, dip=80, rake=164)
    nodal_planes = NodalPlanes(nodal_plane_1=np1, nodal_plane_2=np2)
    taxis = Axis(plunge=40, azimuth=70)
    naxis = Axis(plunge=50, azimuth=80)
    paxis = Axis(plunge=60, azimuth=90)
    paxes = PrincipalAxes(t_axis=taxis,
                          n_axis=naxis,
                          p_axis=paxis)
    focal = FocalMechanism(nodal_planes=nodal_planes,
                           principal_axes=paxes)
    event = Event(focal_mechanisms=[focal])
    catalog = Catalog(events=[event])
    event_text = '''<shakemap-data code_version="4.0" map_version="1">
<earthquake id="us2000cmy3" lat="56.046" lon="-149.073" mag="7.9"
time="2018-01-23T09:31:42Z"
depth="25.00" locstring="280km SE of Kodiak, Alaska" netid="us" network=""/>
</shakemap-data>'''
    try:
        tempdir = tempfile.mkdtemp()
        xmlfile = os.path.join(tempdir, 'quakeml.xml')
        catalog.write(xmlfile, format="QUAKEML")
        eventfile = os.path.join(tempdir, 'event.xml')
        f = open(eventfile, 'wt')
        f.write(event_text)
        f.close()
        params = read_moment_quakeml(xmlfile)
        assert params['moment']['NP1']['strike'] == 259.0
        assert params['moment']['NP1']['dip'] == 74.0
        assert params['moment']['NP1']['rake'] == 10.0
        assert params['moment']['NP2']['strike'] == 166.0
        assert params['moment']['NP2']['dip'] == 80.0
        assert params['moment']['NP2']['rake'] == 164.0
        origin = Origin.fromFile(eventfile, momentfile=xmlfile)
        assert origin.mag == 7.9
        assert origin.lat == 56.046
        assert origin.lon == -149.073
        assert origin.id == 'us2000cmy3'
    except Exception:
        assert False
    finally:
        shutil.rmtree(tempdir)


def test_fromOrientation():
    py = [0, 0.5]
    px = [0, 0.5]
    pz = [10, 20]
    dx = [5, 7]
    dy = [8, 5]
    width = [10, 40]
    length = [20, 50]
    strike = [0, 90]
    dip = [30, 20]

    # Rupture requires an origin even when not used:
    origin = Origin({'id': 'test',
                     'lon': 0, 'lat': 0,
                     'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                     'network': '', 'locstring': '',
                     'time': HistoricTime.utcfromtimestamp(time.time())})
    rupture = QuadRupture.fromOrientation(px, py, pz, dx, dy, length, width,
                                          strike, dip, origin)
    p1 = rupture._geojson['features'][0]['geometry']['coordinates'][0][0][0]
    p2 = rupture._geojson['features'][0]['geometry']['coordinates'][0][0][1]
    p3 = rupture._geojson['features'][0]['geometry']['coordinates'][0][0][2]
    p4 = rupture._geojson['features'][0]['geometry']['coordinates'][0][0][3]
    p5 = rupture._geojson['features'][0]['geometry']['coordinates'][0][1][0]
    p6 = rupture._geojson['features'][0]['geometry']['coordinates'][0][1][1]
    p7 = rupture._geojson['features'][0]['geometry']['coordinates'][0][1][2]
    p8 = rupture._geojson['features'][0]['geometry']['coordinates'][0][1][3]

    # Check depths
    np.testing.assert_allclose(p1[2], 6)
    np.testing.assert_allclose(p2[2], 6)
    np.testing.assert_allclose(p3[2], 11)
    np.testing.assert_allclose(p4[2], 11)
    np.testing.assert_allclose(p5[2], 18.2898992834)
    np.testing.assert_allclose(p6[2], 18.2898992834)
    np.testing.assert_allclose(p7[2], 31.9707050164)
    np.testing.assert_allclose(p8[2], 31.9707050164)

    # Exception raised if no origin
    with pytest.raises(Exception) as a:
        rupture = QuadRupture.fromOrientation(px, py, pz, dx, dy, length,
                                              width, strike, dip, None)
    print(str(a))

    # Exception raised if different lengths of arrays
    with pytest.raises(Exception) as a:
        py = [0, 2]
        px = [0]
        pz = [10]
        dx = [5]
        dy = [8]
        width = [10]
        length = [20]
        strike = [0]
        dip = [30]

        origin = Origin({'id': 'test',
                         'lon': 0, 'lat': 0,
                         'depth': 5.0, 'mag': 7.0, 'netid': 'us',
                         'network': '', 'locstring': '',
                         'time': HistoricTime.utcfromtimestamp(time.time())})
        rupture = QuadRupture.fromOrientation(px, py, pz, dx, dy, length,
                                              width, strike, dip, origin)
    print(str(a))


if __name__ == "__main__":
    test_no_origin()
    test_text_to_json()
    test_rupture_from_dict()
    test_EdgeRupture()
    test_QuadRupture()
    test_rupture_depth(interactive=True)
    test_slip()
    test_northridge()
    test_parse_complicated_rupture()
    test_incorrect()
    test_fromTrace()
    test_with_quakeml()
    test_fromOrientation()
