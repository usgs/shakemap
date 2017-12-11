#!/usr/bin/env python

# stdlib imports
import os
import os.path
import sys
import io
import copy

# third party
import numpy as np
import pytest
from openquake.hazardlib.geo.geodetic import azimuth
from mapio.geodict import GeoDict
import matplotlib.pyplot as plt

from shakelib.rupture.origin import Origin
from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.edge_rupture import EdgeRupture
from shakelib.rupture.factory import get_rupture
from shakelib.rupture.factory import rupture_from_dict

from shakelib.rupture.utils import get_local_unit_slip_vector
from shakelib.rupture.utils import get_quad_slip
from shakelib.rupture.factory import text_to_json

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_rupture_from_dict():
    # Grab an EdgeRupture
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})

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

    # Point rupture
    origin = Origin({'eventsourcecode': 'test',
                     'lon': -122.5, 'lat': 37.3,
                     'depth': 5.0, 'mag': 7.0})
    rup_original = get_rupture(origin)
    d = rup_original._geojson
    rup_from_dict = rupture_from_dict(d)
    assert rup_from_dict.lats == 37.3
    assert rup_from_dict.lons == -122.5


def test_EdgeRupture():
    # Rupture requires an origin even when not used:
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})

    file = os.path.join(homedir, 'rupture_data/cascadia.json')
    rup = get_rupture(origin, file)
    np.testing.assert_allclose(rup.getArea(), 105635.92827547337)

    # Force read Northridge as EdgeRupture
    file = os.path.join(homedir, 'rupture_data/northridge_fault.txt')
    d = text_to_json(file)
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
    d = text_to_json(file)
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
    rrup1 = qrup.computeRrup(lons, lats, deps)
    rrup2 = erup.computeRrup(lons, lats, deps)
    np.testing.assert_allclose(rrup1, rrup2, atol=2e-2)
    rjb1 = qrup.computeRjb(lons, lats, deps)
    rjb2 = erup.computeRjb(lons, lats, deps)
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
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})

    # First with json file
    file = os.path.join(homedir, 'rupture_data/izmit.json')
    rupj = get_rupture(origin, file)
    # Then with text file:
    file = os.path.join(homedir, 'rupture_data/Barkaetal02_fault.txt')
    rupt = get_rupture(origin, file)

    np.testing.assert_allclose(rupj.lats, rupt.lats, atol=1e-5)
    np.testing.assert_allclose(rupj.lons, rupt.lons, atol=1e-5)
    np.testing.assert_allclose(rupj._depth, rupt._depth, atol=1e-5)
    np.testing.assert_allclose(rupt.getArea(), 2391282.2653900268, atol=1e-5)

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
        origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                         'depth': 5.0, 'mag': 7.0})
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
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})
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

    slp = get_local_unit_slip_vector(22, 30, 86).getArray()
    slpd = np.array([0.82714003,  0.38830563,  0.49878203])
    np.testing.assert_allclose(slp, slpd)


def test_northridge():
    rupture_text = """# Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})
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


def test_parse_complicated_rupture():
    rupture_text = """# SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, and W. Page (2002). The Surface Rupture and Slip Distribution of the 17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. Seism. Soc. Am. 92, 43-60.
    40.70985 29.33760 0
    40.72733 29.51528 0
    40.72933 29.51528 20
    40.71185 29.33760 20
    40.70985 29.33760 0
    >
    40.70513 29.61152 0
    40.74903 29.87519 0
    40.75103 29.87519 20
    40.70713 29.61152 20
    40.70513 29.61152 0
    >
    40.72582 29.88662 0
    40.72336 30.11126 0
    40.73432 30.19265 0
    40.73632 30.19265 20
    40.72536 30.11126 20
    40.72782 29.88662 20
    40.72582 29.88662 0
    >
    40.71210 30.30494 0
    40.71081 30.46540 0
    40.70739 30.56511 0
    40.70939 30.56511 20
    40.71281 30.46540 20
    40.71410 30.30494 20
    40.71210 30.30494 0
    >
    40.71621 30.57658 0
    40.70068 30.63731 0
    40.70268 30.63731 20
    40.71821 30.57658 20
    40.71621 30.57658 0
    >
    40.69947 30.72900 0
    40.79654 30.93655 0
    40.79854 30.93655 20
    40.70147 30.72900 20
    40.69947 30.72900 0
    >
    40.80199 30.94688 0
    40.84501 31.01799 0
    40.84701 31.01799 20
    40.80399 30.94688 20
    40.80199 30.94688 0"""  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})
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
    rupture_text = """# Source: Ji, C., D. V. Helmberger, D. J. Wald, and K.-F. Ma (2003). Slip history and dynamic implications of the 1999 Chi-Chi, Taiwan, earthquake, J. Geophys. Res. 108, 2412, doi:10.1029/2002JB001764.
    24.27980 120.72300	0
    24.05000 121.00000	17
    24.07190 121.09300	17
    24.33120 121.04300	17
    24.33120 121.04300	17
    24.27980 120.72300	0
    >
    24.27980 120.72300	0
    23.70000 120.68000	0
    23.60400 120.97200	17
    24.05000 121.00000	17
    24.27980 120.72300	0
    >
    23.60400 120.97200	17
    23.70000 120.68000	0
    23.58850 120.58600	0
    23.40240 120.78900	17
    23.60400 120.97200	17"""  # noqa

    # Rupture requires an origin even when not used:
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})
    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        get_rupture(origin, cbuf)


def test_fromTrace():
    xp0 = [0.0]
    xp1 = [0.0]
    yp0 = [0.0]
    yp1 = [0.05]
    zp = [0.0]
    widths = [10.0]
    dips = [45.0]

    # Rupture requires an origin even when not used:
    origin = Origin({'eventsourcecode': 'test', 'lat': 0, 'lon': 0,
                     'depth': 5.0, 'mag': 7.0})
    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths,
        dips, origin,
        reference='From J Smith, (personal communication)')
    fstr = io.StringIO()
    rupture.writeTextFile(fstr)

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


if __name__ == "__main__":
    test_rupture_from_dict()
    test_rupture_depth(interactive=True)
    test_EdgeRupture()
    test_QuadRupture()
    test_slip()
    test_northridge()
    test_parse_complicated_rupture()
    test_incorrect()
    test_fromTrace()
