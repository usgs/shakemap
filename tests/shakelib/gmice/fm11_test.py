#!/usr/bin/env python

import numpy as np
import pytest

import os
import sys

from openquake.hazardlib.imt import PGA, PGV, SA, MMI

from shakelib.gmice.fm11 import FM11

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
sys.path.insert(0, shakedir)


# Inputs
amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5]))
mmi_in = np.array([2., 3., 4., 5., 6.2, 7.3, 8.1])
dists = np.array([22.2, 10., 1.1, 32., 120., 300.0, 450.5])

lfact = np.log10(np.e)
dmda_target = {PGA(): set([x * lfact for x in [2.58]]),
               PGV(): set([x * lfact for x in [2.35]]),
               SA(0.3): set([x * lfact for x in [2.47]]),
               SA(1.0): set([x * lfact for x in [2.05]]),
               SA(3.0): set([x * lfact for x in [2.00]])}
dadm_target = {PGA(): set([1.0 / (x * lfact) for x in [2.58]]),
               PGV(): set([1.0 / (x * lfact) for x in [2.35]]),
               SA(0.3): set([1.0 / (x * lfact) for x in [2.47]]),
               SA(1.0): set([1.0 / (x * lfact) for x in [2.05]]),
               SA(3.0): set([1.0 / (x * lfact) for x in [2.00]])}


def test_fm11():
    gmice = FM11()

    mi, dmda = gmice.getMIfromGM(amps_in, PGA(), dists=None, mag=None)
    mi_target = np.array(
        [4.238506, 5.469479, 6.818506, 7.595163, 8.049479, 8.371821,
         8.621849])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[PGA()]) == set())

    mi, dmda = gmice.getMIfromGM(
        amps_in + np.log(100), PGV(), dists=dists, mag=None)
    mi_target = np.array(
        [5.11, 6.231235, 7.46, 8.16742, 8.581235, 8.874841,
         9.10258])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[PGV()]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(0.3), dists=None, mag=3.0)
    mi_target = np.array(
        [3.68942245, 4.86791195, 6.15942245, 6.90296654, 7.33791195,
         7.64651063, 7.88587836])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(0.3)]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(1.0), dists=dists, mag=7.5)
    mi_target = np.array(
        [5.152921, 6.13102, 7.202921, 7.820033, 8.18102, 8.437144,
         8.63581])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(1.0)]) == set())

    mi, dmda = gmice.getMIfromGM(amps_in, SA(3.0), dists=None, mag=None)
    mi_target = np.array(
        [6.293338, 7.247581, 8.293338, 8.895398, 9.247581, 9.497458,
         9.691278])
    np.testing.assert_allclose(mi, mi_target)
    assert((set(dmda) - dmda_target[SA(3.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.array([-6.60298051, -5.71050567, -4.81803083,
                            -3.92555598, -2.85458617, -1.87286385,
                            -1.15888397])
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGA()]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGV(), dists=dists, mag=None)
    amps_target = np.array(
        [-3.04725091, -2.06742747, -1.08760402, -0.10778058,  1.06800755,
         2.14581334, 2.9296721])
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGV()]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(0.3), dists=None, mag=6.0)
    amps_target = np.array(
        [-6.18008474, -5.24786405, -4.31564337, -3.38342268, -2.26475786,
         -1.23931511, -0.49353856])
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[SA(0.3)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(1.0), dists=dists, mag=8.0)
    amps_target = np.array([-8.14657017, -7.02335793, -5.90014569,
                            -4.77693345, -3.42907876, -2.19354529,
                            -1.2949755])
    assert((set(dadm) - dadm_target[SA(1.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, SA(3.0), dists=None, mag=5.0)
    amps_target = np.array([-9.54805824, -8.3967657,  -7.24547315, -6.0941806,
                            -4.71262955, -3.44620775, -2.52517371])
    np.testing.assert_allclose(amps, amps_target)
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[SA(3.0)]) == set())

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=7.0)
    amps_target = np.array([-6.60298051, -5.71050567, -4.81803083,
                            -3.92555598, -2.85458617, -1.87286385,
                            -1.15888397])
    np.testing.assert_allclose(amps, amps_target)
    assert((set(dadm) - dadm_target[PGA()]) == set())

    sdd1 = gmice.getGM2MIsd()
    assert sdd1[PGA()] == 0.18
    assert sdd1[PGV()] == 0.14
    assert sdd1[SA(0.3)] == 0.3
    assert sdd1[SA(1.0)] == 0.21
    assert sdd1[SA(3.0)] == 0.14

    sdd2 = gmice.getMI2GMsd()
    np.testing.assert_allclose(sdd2[PGA()], 0.7138013788281542)
    np.testing.assert_allclose(sdd2[PGV()], 0.5065687204586901)
    np.testing.assert_allclose(sdd2[SA(0.3)], 0.9670857390574993)
    np.testing.assert_allclose(sdd2[SA(1.0)], 0.7138013788281542)
    np.testing.assert_allclose(sdd2[SA(3.0)], 0.598672124178452)

    nm = gmice.getName()
    assert nm == 'Faenza and Michelini (2010, 2011)'

    sc = gmice.getScale()
    assert sc == 'scale_fm11.ps'

    mm = gmice.getMinMax()
    assert mm == (1.0, 10.0)

    dt = gmice.getDistanceType()
    assert dt == 'rrup'

    #
    # This should fail
    #
    with pytest.raises(ValueError):  # noqa
        mi, dmda = gmice.getMIfromGM(amps_in, MMI(), dists=None, mag=None)


if __name__ == '__main__':
    test_fm11()
