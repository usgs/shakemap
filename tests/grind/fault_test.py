#!/usr/bin/env python

# stdlib imports
import os
import os.path
import sys
import io

# third party
import numpy as np
import pytest

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

from shakemap.grind.fault import Fault
from shakemap.utils.exception import ShakeMapException
from shakemap.utils.misc import get_command_output
from shakemap.grind.fault import get_local_unit_slip_vector
from shakemap.grind.fault import get_quad_slip

def test_pisgah_bullion_mtn(tmpdir):
    # a segment of this fault causes a division by zero error that
    # we trap for and are testing here.

    # make a temporary directory
    p = tmpdir.mkdir("sub")
    jsonfile = os.path.join(shakedir, 'tests/data/eventdata/UCERF3_EventSet_All.json')
    cmd = 'mkinputdir -f %s -i 46 -s %s' %(jsonfile, p)
    rc,so,se = get_command_output(cmd)
    if se != b'':
        print(so.decode())
        print(se.decode())
    assert se == b''

def test_misc():
    # Make a fault
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])
    flt = Fault.fromTrace(lon0, lat0, lon1, lat1, z, W, dip)
    fm = flt.getFaultAsMesh()
    fa = flt.getFaultAsArrays()
    ref = flt.getReference()


def test_slip():
    # Make a fault
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])
    flt = Fault.fromTrace(lon0, lat0, lon1, lat1, z, W, dip)

    slp = get_quad_slip(flt.getQuadrilaterals()[0], 30).getArray()
    slpd = np.array([0.80816457,  0.25350787,  0.53160491])
    np.testing.assert_allclose(slp, slpd)

    slp = get_local_unit_slip_vector(22, 30, 86).getArray()
    slpd = np.array([0.82714003,  0.38830563,  0.49878203])
    np.testing.assert_allclose(slp, slpd)


def test_northridge():
    fault_text = """
    # Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """
    cbuf = io.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
    strike = fault.getStrike()
    np.testing.assert_allclose(strike, 122.06408, atol=0.001)
    dip = fault.getDip()
    np.testing.assert_allclose(dip, 40.20979, atol=0.001)
    L = fault.getFaultLength()
    np.testing.assert_allclose(L, 17.99198, atol=0.001)
    W = fault.getWidth()
    np.testing.assert_allclose(W, 23.93699, atol=0.001)
    nq = fault.getNumQuads()
    np.testing.assert_allclose(nq, 1)
    ns = fault.getNumSegments()
    np.testing.assert_allclose(ns, 1)
    sind = fault._getSegmentIndex()
    np.testing.assert_allclose(sind, [0])
    ztor = fault.getTopOfRupture()
    np.testing.assert_allclose(ztor, 5, atol=0.001)
    itl = fault.getIndividualTopLengths()
    np.testing.assert_allclose(itl, 17.9919846, atol=0.001)
    iw = fault.getIndividualWidths()
    np.testing.assert_allclose(iw, 23.93699668, atol=0.001)
    lats = fault.getLats()
    lats_d = np.array([34.315,  34.401,  34.261,  34.175,  34.315])
    np.testing.assert_allclose(lats, lats_d, atol=0.001)
    lons = fault.getLons()
    lons_d = np.array([-118.421, -118.587, -118.693, -118.527, -118.421])
    np.testing.assert_allclose(lons, lons_d, atol=0.001)


def parse_complicated_fault():
    fault_text = """#SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, and W. Page (2002). The Surface Rupture and Slip Distribution of the 17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. Seism. Soc. Am. 92, 43-60.
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

    cbuf = io.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
    strike = fault.getStrike()
    np.testing.assert_allclose(strike, -100.464330, atol=0.001)
    dip = fault.getDip()
    np.testing.assert_allclose(dip, 89.3985, atol=0.001)
    L = fault.getFaultLength()
    np.testing.assert_allclose(L, 119.5578, atol=0.001)
    W = fault.getWidth()
    np.testing.assert_allclose(W, 20.001, atol=0.001)
    nq = fault.getNumQuads()
    np.testing.assert_allclose(nq, 9)
    ns = fault.getNumSegments()
    np.testing.assert_allclose(ns, 7)
    sind = fault._getSegmentIndex()
    np.testing.assert_allclose(sind, [0, 1, 2, 2, 3, 3, 4, 5, 6])
    ztor = fault.getTopOfRupture()
    np.testing.assert_allclose(ztor, 0, atol=0.001)
    itl = fault.getIndividualTopLengths()
    itl_d = np.array([15.13750778,  22.80237887,  18.98053425,   6.98263853,
                      13.55978731,   8.43444811,   5.41399812,  20.57788056,
                      7.66869463])
    np.testing.assert_allclose(itl, itl_d, atol=0.001)
    iw = fault.getIndividualWidths()
    iw_d = np.array([20.00122876,  20.00122608,  20.00120173,  20.00121028,
                     20.00121513,  20.00121568,  20.00107293,  20.00105498,
                     20.00083348])
    np.testing.assert_allclose(iw, iw_d, atol=0.001)
    lats = fault.getLats()
    lats_d = np.array([40.70985, 40.72733, 40.72933, 40.71185, 40.70985,
                       np.nan, 40.70513, 40.74903, 40.75103, 40.70713,
                       40.70513, np.nan, 40.72582, 40.72336, 40.73432,
                       40.73632, 40.72536, 40.72782, 40.72582, np.nan,
                       40.7121, 40.71081, 40.70739, 40.70939, 40.71281,
                       40.7141, 40.7121, np.nan, 40.71621, 40.70068,
                       40.70268, 40.71821, 40.71621, np.nan, 40.69947,
                       40.79654, 40.79854, 40.70147, 40.69947, np.nan,
                       40.80199, 40.84501, 40.84701, 40.80399, 40.80199])
    np.testing.assert_allclose(lats, lats_d, atol=0.001)
    lons = fault.getLons()
    lons_d = np.array([29.3376, 29.51528, 29.51528, 29.3376, 29.3376,
                       np.nan, 29.61152, 29.87519, 29.87519, 29.61152,
                       29.61152, np.nan, 29.88662, 30.11126, 30.19265,
                       30.19265, 30.11126, 29.88662, 29.88662, np.nan,
                       30.30494, 30.4654, 30.56511, 30.56511, 30.4654,
                       30.30494, 30.30494, np.nan, 30.57658, 30.63731,
                       30.63731, 30.57658, 30.57658, np.nan, 30.729,
                       30.93655, 30.93655, 30.729, 30.729, np.nan,
                       30.94688,  31.01799, 31.01799, 30.94688, 30.94688])
    np.testing.assert_allclose(lons, lons_d, atol=0.001)


def test_incorrect():
    fault_text = """# Source: Ji, C., D. V. Helmberger, D. J. Wald, and K.-F. Ma (2003). Slip history and dynamic implications of the 1999 Chi-Chi, Taiwan, earthquake, J. Geophys. Res. 108, 2412, doi:10.1029/2002JB001764.
    24.27980 120.72300	0 
    24.05000 121.00000	17
    24.07190 121.09300	17
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

    cbuf = io.StringIO(fault_text)
    with pytest.raises(ShakeMapException):
        fault = Fault.readFaultFile(cbuf)


def test_fromTrace():
    xp0 = [0.0]
    xp1 = [0.0]
    yp0 = [0.0]
    yp1 = [0.05]
    zp = [0.0]
    widths = [10.0]
    dips = [45.0]

    fault = Fault.fromTrace(xp0, yp0, xp1, yp1, zp, widths,
                            dips, reference='From J Smith, (personal communication)')
    fstr = io.StringIO()
    fault.writeFaultFile(fstr)

    xp0 = [-121.81529, -121.82298]
    xp1 = [-121.82298, -121.83068]
    yp0 = [37.73707, 37.74233]
    yp1 = [37.74233, 37.74758]
    zp = [10, 15]
    widths = [15.0, 20.0]
    dips = [30.0, 45.0]
    fault = Fault.fromTrace(xp0, yp0, xp1, yp1, zp, widths,
                            dips, reference='From J Smith, (personal communication)')
