#!/usr/bin/env python

import numpy as np
import pytest

from openquake.hazardlib.imt import PGA, PGV, MMI

from shakelib.gmice.wald99 import Wald99


do_test = True

# Inputs
amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0]))
mmi_in = np.array([2., 3., 4., 5., 6.2, 7.3, 8.1, 9.0])

lfact = np.log10(np.e)
dmda_target = {PGA(): set([x * lfact for x in [3.66, 2.20]]),
               PGV(): set([x * lfact for x in [3.47, 2.10]])}
dadm_target = {PGA(): set([1.0 / (x * lfact) for x in [3.66, 2.20]]),
               PGV(): set([1.0 / (x * lfact) for x in [3.47, 2.10]])}


def test_wald99():
    gmice = Wald99()

    df = {'PGA': amps_in,
          'PGV': amps_in + np.log(100)}

    mi = gmice.getPreferredMI(df)
    mi_target = np.array(
        [3.18167182, 4.23133858, 5.68946656, 6.84666436, 7.47561075,
         7.90914817, 8.24542592, 9.29])

    if do_test is True:
        np.testing.assert_allclose(mi, mi_target)
    else:
        print(repr(mi))

    pga_amps_in_nan = np.log(np.array(
        [0.01, 0.03, 0.1, np.nan, 0.3, 0.4, 0.5, 1.0]))
    pgv_amps_in_nan = np.log(np.array(
        [0.01, 0.03, 0.1, 0.2, np.nan, 0.4, 0.5, 1.0])) + np.log(100)
    pga_amps_in_nan[5] = np.nan
    pgv_amps_in_nan[5] = np.nan
    df = {'PGA': pga_amps_in_nan,
          'PGV': pgv_amps_in_nan}

    mi = gmice.getPreferredMI(df)
    mi_target = np.array(
        [3.18167182, 4.23133858, 5.68946656, 6.86457408, 7.37577236,
         np.nan, 8.24542592, 9.29])

    if do_test is True:
        np.testing.assert_allclose(mi, mi_target)
    else:
        print(repr(mi))

    mi, dmda = gmice.getMIfromGM(amps_in, PGA(), dists=None, mag=None)
    mi_target = np.array(
        [3.18167182,  4.23133858,  5.62950857,  6.73127835,  7.37577236,
         7.83304814,  8.18773878,  9.28950857])

    if do_test is True:
        np.testing.assert_allclose(mi, mi_target)
        assert((set(dmda) - dmda_target[PGA()]) == set())
    else:
        print(repr(mi))

    mi, dmda = gmice.getMIfromGM(
        amps_in + np.log(100), PGV(), dists=None, mag=None)
    mi_target = np.array(
        [3.4,  4.40195463,  5.82,  6.86457408,  7.47561075,
         7.90914817,  8.24542592,  9.29])

    if do_test is True:
        np.testing.assert_allclose(mi, mi_target)
        assert((set(dmda) - dmda_target[PGV()]) == set())
    else:
        print(repr(mi))

    amps, dadm = gmice.getGMfromMI(mmi_in, PGA(), dists=None, mag=None)
    amps_target = np.array(
        [-5.84194287, -4.79531328, -3.7486837, -2.69862254, -1.9436766,
         -1.25164283, -0.74834554, -0.1821361])

    if do_test is True:
        np.testing.assert_allclose(amps, amps_target)
        assert((set(dadm) - dadm_target[PGA()]) == set())
    else:
        print(repr(amps))

    amps, dadm = gmice.getGMfromMI(mmi_in, PGV(), dists=None, mag=None)
    amps_target = np.array(
        [-1.53505673, -0.43858764,  0.65788146,  1.75845836,  2.55474139,
         3.2846675,  3.81552285,  4.41273512])

    if do_test is True:
        np.testing.assert_allclose(amps, amps_target)
        assert((set(dadm) - dadm_target[PGV()]) == set())
    else:
        print(repr(amps))

    if do_test is True:
        sdd = gmice.getGM2MIsd()
        assert sdd[PGA()] == 1.08
        assert sdd[PGV()] == 0.98

        sdd = gmice.getMI2GMsd()
        lnten = np.log(10.0)
        np.testing.assert_allclose(sdd[PGA()], 0.295 * lnten)
        np.testing.assert_allclose(sdd[PGV()], 0.282 * lnten)

        nm = gmice.getName()
        assert nm == 'Wald et al. (1999)'

        sc = gmice.getScale()
        assert sc == 'scale_wald99.ps'

        mm = gmice.getMinMax()
        assert mm == (1.0, 10.0)

        dt = gmice.getDistanceType()
        assert dt == 'rrup'

    #
    # This should fail because MMI() is not a valid argument to
    # getMIfromGM
    #
    with pytest.raises(ValueError) as e:  # noqa
        mi, dmda = gmice.getMIfromGM(amps_in, MMI(), dists=None, mag=None)


if __name__ == '__main__':
    test_wald99()
