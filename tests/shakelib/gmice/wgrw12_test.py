#!/usr/bin/env python

import numpy as np
import pytest

import os
import sys

from openquake.hazardlib.imt import PGA, PGV, SA, MMI

from shakelib.gmice.wgrw12 import WGRW12

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
sys.path.insert(0, shakedir)


# Inputs
amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5]))
mmi_in = np.array([2., 3., 4., 5., 6.2, 7.3, 8.1])
dists = np.array([22.2, 10., 1.1, 32., 120., 300.0, 450.5])

lfact = np.log10(np.e)
dmda_target = {PGA(): set([x * lfact for x in [2.08, 1.55, 3.70]]),
               PGV(): set([x * lfact for x in [2.17, 1.47, 3.16]]),
               SA(0.3): set([x * lfact for x in [1.92, 1.69, 4.14]]),
               SA(1.0): set([x * lfact for x in [2.17, 1.51, 2.90]]),
               SA(3.0): set([x * lfact for x in [3.45, 1.17, 3.01]])}
dadm_target = {PGA(): set([1.0 / (x * lfact) for x in [2.08, 1.55, 3.70]]),
               PGV(): set([1.0 / (x * lfact) for x in [2.17, 1.47, 3.16]]),
               SA(0.3): set([1.0 / (x * lfact) for x in [1.92, 1.69, 4.14]]),
               SA(1.0): set([1.0 / (x * lfact) for x in [2.17, 1.51, 2.90]]),
               SA(3.0): set([1.0 / (x * lfact) for x in [3.45, 1.17, 3.01]])}


def test_wgrw12():
    gmice = WGRW12()

    mi, dmda = gmice.getMIfromGM(amps_in, PGA(), dists=None, mag=None)
    mi_target = np.array(
        [3.31708696,  4.05662491,  5.76917533,  6.88298631,  7.53452397,
         7.9967973,  8.35536434])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[PGA()]) == set())

    mi, dmda = gmice.getMIfromGM(
        amps_in + np.log(100), PGV(), dists=dists, mag=None)
    mi_target = np.array(
        [3.78,  4.48136824,  6.05,  7.00125479,  7.55770316,
         7.95250957,  8.25874521])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[PGV()]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(0.3), dists=None, mag=3.0)
    mi_target = np.array(
        [2.93592062,  3.74225554,  4.62592062,  5.34177387,  6.07079169,
         6.58803805,  6.98924551])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(0.3)]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(1.0), dists=dists, mag=7.5)
    mi_target = np.array(
        [3.49070724,  4.3808733,  5.63884012,  6.26430362,  6.49369295,
         6.66102468,  6.94206372])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(1.0)]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(3.0), dists=None, mag=None)
    mi_target = np.array(
        [4.97492371,   6.41105869,   7.98492371,   8.891024,
         9.42105869,   9.79712429,  10.])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(3.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.log(np.array(
        [0.14134045,   0.6243495,   2.75796695,   6.19604804,
         13.07492182,  25.92605261,  42.65329774]) / 100.0)
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGA()]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGV(), dists=dists, mag=None)
    amps_target = np.log(np.array(
        [0.06153407,   0.29470517,   1.41143169,   4.65287643,
         11.15496866,  24.86392117,  44.53835548]))
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGV()]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(0.3), dists=None, mag=6.0)
    amps_target = np.log(np.array(
        [0.27938354,   1.09123124,   4.26218963,  16.53773087,
         32.23524628,  59.43352318,  92.74023471]) / 100.0)
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[SA(0.3)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(1.0), dists=dists, mag=8.0)
    amps_target = np.log(np.array(
        [1.02985637e-01,   3.65288187e-01,   1.67836836e+00,
         7.32931235e+00,   2.37600759e+01,   6.64349034e+01,
         1.25388693e+02]) / 100.0)
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[SA(1.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(3.0), dists=None, mag=5.0)
    amps_target = np.log(np.array(
        [2.89282689e-03,   2.07025242e-02,   1.48157675e-01,
         1.01936799e+00,   2.55271358e+00,   5.92175716e+00,
         1.09202184e+01]) / 100.0)
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[SA(3.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.log(np.array(
        [0.14134045,   0.6243495,   2.75796695,   6.19604804,
         13.07492182,  25.92605261,  42.65329774]) / 100.0)
    print(repr(100 * np.exp(amps)))
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGA()]) == set())

    sdd = gmice.getGM2MIsd()
    assert sdd[PGA()] == 0.66
    assert sdd[PGV()] == 0.63
    assert sdd[SA(0.3)] == 0.82
    assert sdd[SA(1.0)] == 0.75
    assert sdd[SA(3.0)] == 0.89

    sdd = gmice.getMI2GMsd()
    lnten = np.log(10.0)
    np.testing.assert_allclose(sdd[PGA()], np.log10(2.238721) * lnten)
    np.testing.assert_allclose(sdd[PGV()], np.log10(2.39883291) * lnten)
    np.testing.assert_allclose(sdd[SA(0.3)], np.log10(2.7542287) * lnten)
    np.testing.assert_allclose(sdd[SA(1.0)], np.log10(2.951209) * lnten)
    np.testing.assert_allclose(sdd[SA(3.0)], np.log10(4.365158) * lnten)

    nm = gmice.getName()
    assert nm == 'Worden et al. (2012)'

    sc = gmice.getScale()
    assert sc == 'scale_wgrw12.ps'

    mm = gmice.getMinMax()
    assert mm == (1.0, 10.0)

    dt = gmice.getDistanceType()
    assert dt == 'rrup'

    #
    # This should fail
    #
    with pytest.raises(ValueError) as e:  # noqa
        mi, dmda = gmice.getMIfromGM(amps_in, MMI(), dists=None, mag=None)


if __name__ == '__main__':
    test_wgrw12()
