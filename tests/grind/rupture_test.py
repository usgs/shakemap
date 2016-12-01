#!/usr/bin/env python

# stdlib imports
import os
import os.path
import sys
import io

# third party
import numpy as np
import pytest

from shakemap.grind.rupture import QuadRupture
from shakemap.grind.rupture import EdgeRupture
from shakemap.grind.rupture import read_rupture_file
from shakemap.utils.exception import ShakeMapException
from impactutils.io.cmd import get_command_output
from shakemap.grind.rupture import get_local_unit_slip_vector
from shakemap.grind.rupture import get_quad_slip
from shakemap.grind.rupture import text_to_json

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

def test_EdgeRupture():
    file = 'tests/data/cascadia.json'
    rup = read_rupture_file(file)

    # Force read Northridge as EdgeRupture
    file = 'tests/data/eventdata/northridge/northridge_fault.txt'
    d = text_to_json(file)
    rupt = EdgeRupture.fromJson(d)
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
    file = 'tests/data/eventdata/northridge/northridge_fixed_fault.txt'
    d = text_to_json(file)
    rupt = EdgeRupture.fromJson(d)
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


def test_QuadRupture():
    # First with json file
    file = 'tests/data/izmit.json'
    rupj = read_rupture_file(file)
    # Then with text file:
    file = 'tests/data/Barkaetal02_fault.txt'
    rupt = read_rupture_file(file)


def test_misc():
    # Make a rupture
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])
    rup = QuadRupture.fromTrace(lon0, lat0, lon1, lat1, z, W, dip)
    fm = rup.getRuptureAsMesh()
    fa = rup.getRuptureAsArrays()



def test_slip():
    # Make a rupture
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])
    rup = QuadRupture.fromTrace(lon0, lat0, lon1, lat1, z, W, dip)

    slp = get_quad_slip(rup.getQuadrilaterals()[0], 30).getArray()
    slpd = np.array([0.80816457,  0.25350787,  0.53160491])
    np.testing.assert_allclose(slp, slpd)

    slp = get_local_unit_slip_vector(22, 30, 86).getArray()
    slpd = np.array([0.82714003,  0.38830563,  0.49878203])
    np.testing.assert_allclose(slp, slpd)


def test_northridge():
    rupture_text = """
    # Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """
    cbuf = io.StringIO(rupture_text)
    rupture = read_rupture_file(cbuf)
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
    lats = rupture.getLats()
    lats_d = np.array([34.315,  34.401,  34.261,  34.175,  np.nan])
    np.testing.assert_allclose(lats, lats_d, atol=0.001)
    lons = rupture.getLons()
    lons_d = np.array([-118.421, -118.587, -118.693, -118.527, np.nan])
    np.testing.assert_allclose(lons, lons_d, atol=0.001)


def parse_complicated_rupture():
    rupture_text = """#SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, and W. Page (2002). The Surface Rupture and Slip Distribution of the 17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. Seism. Soc. Am. 92, 43-60.
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
    40.80199 30.94688 0"""

    cbuf = io.StringIO(rupture_text)
    rupture = read_rupture_file(cbuf)
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
    lats = rupture.getLats()
    lats_d = np.array([
        40.70985,  40.72733,  40.72933,  40.71185,       np.nan,  40.70513,
        40.74903,  40.75103,  40.70713,       np.nan,  40.72582,  40.72336,
        40.73432,  40.73632,  40.72536,  40.72782,       np.nan,  40.7121 ,
        40.71081,  40.70739,  40.70939,  40.71281,  40.7141 ,       np.nan,
        40.71621,  40.70068,  40.70268,  40.71821,       np.nan,  40.69947,
        40.79654,  40.79854,  40.70147,       np.nan,  40.80199,  40.84501,
        40.84701,  40.80399,       np.nan])
    np.testing.assert_allclose(lats, lats_d, atol=0.001)
    lons = rupture.getLons()
    lons_d = np.array([
        29.3376 ,  29.51528,  29.51528,  29.3376 ,       np.nan,  29.61152,
        29.87519,  29.87519,  29.61152,       np.nan,  29.88662,  30.11126,
        30.19265,  30.19265,  30.11126,  29.88662,       np.nan,  30.30494,
        30.4654 ,  30.56511,  30.56511,  30.4654 ,  30.30494,       np.nan,
        30.57658,  30.63731,  30.63731,  30.57658,       np.nan,  30.729  ,
        30.93655,  30.93655,  30.729  ,       np.nan,  30.94688,  31.01799,
        31.01799,  30.94688,       np.nan])
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
    23.60400 120.97200	17"""

    cbuf = io.StringIO(rupture_text)
    with pytest.raises(Exception):
        rupture = read_rupture_file(cbuf)


def test_fromTrace():
    xp0 = [0.0]
    xp1 = [0.0]
    yp0 = [0.0]
    yp1 = [0.05]
    zp = [0.0]
    widths = [10.0]
    dips = [45.0]

    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths,
        dips, reference='From J Smith, (personal communication)')
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
        dips, reference='From J Smith, (personal communication)')


if __name__ == "__main__":
    test_misc()
    test_slip()
    test_northridge()
    parse_complicated_rupture()
    test_incorrect()
    test_fromTrace()
    
