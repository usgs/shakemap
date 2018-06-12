#!/usr/bin/env python

import numpy as np
import pytest

import os
import sys

from openquake.hazardlib.imt import PGA, PGV, SA, MMI

from shakelib.gmice.ak07 import AK07

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
sys.path.insert(0, shakedir)


# Inputs
amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5]))
mmi_in = np.array([2., 3., 4., 5., 6.2, 7.3, 8.1])
dists = np.array([22.2, 10., 1.1, 32., 120., 300.0, 450.5])

lfact = np.log10(np.e)
dmda_target = {PGA(): set([x * lfact for x in [1.55, 3.70]]),
               PGV(): set([x * lfact for x in [1.47, 3.16]]),
               SA(0.3): set([x * lfact for x in [1.69, 4.14]]),
               SA(1.0): set([x * lfact for x in [1.51, 2.90]]),
               SA(3.0): set([x * lfact for x in [1.17, 3.01]])}
dadm_target = {PGA(): set([1.0 / (x * lfact) for x in [1.55, 3.70]]),
               PGV(): set([1.0 / (x * lfact) for x in [1.47, 3.16]]),
               SA(0.3): set([1.0 / (x * lfact) for x in [1.69, 4.14]]),
               SA(1.0): set([1.0 / (x * lfact) for x in [1.51, 2.90]]),
               SA(3.0): set([1.0 / (x * lfact) for x in [1.17, 3.01]])}


def test_ak07():
    gmice = AK07()

    mi, dmda = gmice.getMIfromGM(amps_in, PGA(), dists=None, mag=None)
    mi_target = np.array(
        [4.02841992026,  4.69161846432,  6.23592624018,  7.46713892245,
         8.18735217199, 8.69835160472,  9.09471355792])
    np.testing.assert_allclose(mi, mi_target)

    mi, dmda = gmice.getMIfromGM(
        amps_in + np.log(100), PGV(), dists=dists, mag=None)
    mi_target = np.array(
        [4.37, 4.99980005623, 6.57, 7.48212088686, 8.0156774018,
         8.39424177372, 8.68787911314])
    np.testing.assert_allclose(mi, mi_target)

    mi, dmda = gmice.getMIfromGM(amps_in, SA(0.3), dists=None, mag=3.0)
    mi_target = np.array(
        [3.74866985004, 4.39755475646, 5.26034166627, 6.33200845084,
         6.95889333307, 7.4036752354, 7.74867488171])
    np.testing.assert_allclose(mi, mi_target)

    mi, dmda = gmice.getMIfromGM(amps_in, SA(1.0), dists=dists, mag=7.5)
    mi_target = np.array(
        [3.52702354769, 4.07617250928, 5.55842357177, 6.46666805811,
         7.00909852304, 7.39358539638, 7.67946993475])
    np.testing.assert_allclose(mi, mi_target)

    mi, dmda = gmice.getMIfromGM(amps_in, SA(3.0), dists=None, mag=None)
    mi_target = np.array(
        [4.99925301952, 6.3963707863, 7.96500702214, 8.86809700913,
         9.3963707863, 9.77118699612, 10.])
    np.testing.assert_allclose(mi, mi_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.log(np.array(
        [0.000347300247926,   0.00182024377197, 0.00954012388176,
         0.0498674941261, 0.0979977440572, 0.182039122534,
         0.285603651352]))
    np.testing.assert_allclose(amps, amps_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, PGV(), dists=dists, mag=None)
    amps_target = np.log(np.array(
        [0.0160156826445, 0.0916476244071, 0.524441401963,
         3.03283082017, 7.54897155245, 17.4150246057,
         31.9853049046]))
    np.testing.assert_allclose(amps, amps_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(0.3), dists=None, mag=6.0)
    amps_target = np.log(np.array(
        [0.000517861166862, 0.00281518839281, 0.0153038810286,
         0.0845026480232, 0.183632256171, 0.374056955014,
         0.627562284541]))
    np.testing.assert_allclose(amps, amps_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(1.0), dists=dists, mag=8.0)
    amps_target = np.log(np.array(
        [0.000508056775718, 0.00367375892063, 0.0258564132795,
         0.0636582041631, 0.159533013524, 0.371822895172,
         0.694260679175]))
    np.testing.assert_allclose(amps, amps_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(3.0), dists=None, mag=5.0)
    amps_target = np.log(np.array(
        [0.0000473148708829, 0.000281962568289, 0.00168029392098,
         0.0102722203276, 0.0258026508623, 0.0600248374471,
         0.110916883717]))
    np.testing.assert_allclose(amps, amps_target)

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.log(np.array(
        [0.000347300247926, 0.00182024377197, 0.00954012388176,
         0.0498674941261, 0.0979977440572, 0.182039122534,
         0.285603651352]))
    # print(repr(100 * np.exp(amps)))
    np.testing.assert_allclose(amps, amps_target)

    sdd = gmice.getGM2MIsd()
    assert sdd[PGA()] == 0.89
    assert sdd[PGV()] == 0.76
    assert sdd[SA(0.3)] == 0.79
    assert sdd[SA(1.0)] == 0.73
    assert sdd[SA(3.0)] == 0.72

    sdd = gmice.getMI2GMsd()
    assert abs(sdd[PGA()] - 1.312473503006606) < 0.0000001
    assert abs(sdd[PGV()] - 1.197344248356904) < 0.0000001
    assert abs(sdd[SA(0.3)] - 1.450628608586249) < 0.0000001
    assert abs(sdd[SA(1.0)] - 1.312473503006606) < 0.0000001
    assert abs(sdd[SA(3.0)] - 1.865093925325177) < 0.0000001

    nm = gmice.getName()
    assert nm == 'Atkinson and Kaka (2007)'

    sc = gmice.getScale()
    assert sc == 'scale_ak07.ps'

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
    test_ak07()
