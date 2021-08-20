#!/usr/bin/env python

# Standard library imports
import os.path
import sys

# Third party imports
import numpy as np
import pytest

# Local imports
from shakelib.conversions.imt.abrahamson_bhasin_2020 import (
    AbrahamsonBhasin2020, AbrahamsonBhasin2020PGA, AbrahamsonBhasin2020SA1)

from openquake.hazardlib import const


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..', '..'))
sys.path.insert(0, shakedir)


def test_abrahamsonbhasin2020():
    # Inputs
    mag = 8.0
    ab2020 = AbrahamsonBhasin2020(mag)
    Tref = ab2020.getTref()
    psa = np.array([2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1])
    sigma = np.array([0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72,
                      0.72, 0.72])
    tau = np.array([0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])
    phi = np.array([0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6])
    sdtypes = [const.StdDev.TOTAL, const.StdDev.INTER_EVENT,
               const.StdDev.INTRA_EVENT]
    rrup = np.array([2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                     90.0, 100.0])
    vs30 = np.array([200.0, 180.0, 420.0, 630.0, 740.0, 550.0, 360.0, 270.0,
                     480.0, 290.0, 600.0])

    pgv, stddevs = ab2020.getPGVandSTDDEVS(psa, [sigma, tau, phi], sdtypes,
                                           rrup, vs30)
    # print(repr(pgv))
    # print(repr(stddevs[0]))
    # print(repr(stddevs[1]))
    # print(repr(stddevs[2]))

    pgv_ref = \
        np.array([6.76516894, 6.47038352, 6.14539197, 5.88822518, 5.66883031,
                  5.4931866 , 5.33554182, 5.17621557, 4.96628769, 4.83086164,
                  4.68435124])
    np.testing.assert_almost_equal(pgv, pgv_ref)
    sig_ref = \
        np.array([0.57500259, 0.57500259, 0.57500259, 0.57500259, 0.57500259,
                  0.57500259, 0.57500259, 0.57500259, 0.57500259, 0.57500259,
                   0.57500259])
    np.testing.assert_almost_equal(sig_ref, stddevs[0])
    tau_ref = \
        np.array([0.30665055, 0.30665055, 0.30665055, 0.30665055, 0.30665055,
                  0.30665055, 0.30665055, 0.30665055, 0.30665055, 0.30665055,
                  0.30665055])
    np.testing.assert_almost_equal(tau_ref, stddevs[1])
    phi_ref = \
        np.array([0.48793213, 0.48793213, 0.48793213, 0.48793213, 0.48793213,
                  0.48793213, 0.48793213, 0.48793213, 0.48793213, 0.48793213,
                  0.48793213])
    np.testing.assert_almost_equal(phi_ref, stddevs[2])


def test_abrahamsonbhasin2020sa1():
    # Inputs
    mag = 7.4
    ab2020 = AbrahamsonBhasin2020SA1(mag)
    Tref = ab2020.getTref()
    psa = np.array([2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1])
    sigma = np.array([0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72,
                      0.72, 0.72])
    tau = np.array([0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])
    phi = np.array([0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6])
    sdtypes = [const.StdDev.TOTAL, const.StdDev.INTER_EVENT,
               const.StdDev.INTRA_EVENT]
    rrup = np.array([2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                     90.0, 100.0])
    vs30 = np.array([200.0, 180.0, 420.0, 630.0, 740.0, 550.0, 360.0, 270.0,
                     480.0, 290.0, 600.0])

    pgv, stddevs = ab2020.getPGVandSTDDEVS(psa, [sigma, tau, phi], sdtypes,
                                           rrup, vs30)
    # print(repr(pgv))
    # print(repr(stddevs[0]))
    # print(repr(stddevs[1]))
    # print(repr(stddevs[2]))

    pgv_ref = \
        np.array([5.61589861, 5.31341327, 4.86006093, 4.54874047, 4.31509445,
               4.1939768 , 4.11077608, 4.0071555 , 3.72850462, 3.68138081,
               3.43577395])
    np.testing.assert_almost_equal(pgv, pgv_ref)
    sig_ref = \
        np.array([0.58261055, 0.58261055, 0.58261055, 0.58261055, 0.58261055,
               0.58261055, 0.58261055, 0.58261055, 0.58261055, 0.58261055,
               0.58261055])
    np.testing.assert_almost_equal(sig_ref, stddevs[0])
    tau_ref = \
        np.array([0.28145952, 0.28145952, 0.28145952, 0.28145952, 0.28145952,
               0.28145952, 0.28145952, 0.28145952, 0.28145952, 0.28145952,
               0.28145952])
    np.testing.assert_almost_equal(tau_ref, stddevs[1])
    phi_ref = \
        np.array([0.50756161, 0.50756161, 0.50756161, 0.50756161, 0.50756161,
               0.50756161, 0.50756161, 0.50756161, 0.50756161, 0.50756161,
               0.50756161])
    np.testing.assert_almost_equal(phi_ref, stddevs[2])


def test_abrahamsonbhasin2020pga():
    # Inputs
    mag = 7.4
    ab2020 = AbrahamsonBhasin2020PGA(mag)
    Tref = ab2020.getTref()
    psa = np.array([2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.1])
    sigma = np.array([0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72,
                      0.72, 0.72])
    tau = np.array([0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4])
    phi = np.array([0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6])
    sdtypes = [const.StdDev.TOTAL, const.StdDev.INTER_EVENT,
               const.StdDev.INTRA_EVENT]
    rrup = np.array([2.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0,
                     90.0, 100.0])
    vs30 = np.array([200.0, 180.0, 420.0, 630.0, 740.0, 550.0, 360.0, 270.0,
                     480.0, 290.0, 600.0])

    pgv, stddevs = ab2020.getPGVandSTDDEVS(psa, [sigma, tau, phi], sdtypes,
                                           rrup, vs30)
    # print(repr(pgv))
    # print(repr(stddevs[0]))
    # print(repr(stddevs[1]))
    # print(repr(stddevs[2]))

    pgv_ref = \
        np.array([5.64284705, 5.40573036, 4.7921966 , 4.4158703 , 4.17001613,
                  4.13980353, 4.17526695, 4.15793101, 3.76625139, 3.85369198,
                  3.45235077])
    np.testing.assert_almost_equal(pgv, pgv_ref)
    sig_ref = \
        np.array([0.60554952, 0.60554952, 0.60554952, 0.60554952, 0.60554952,
                  0.60554952, 0.60554952, 0.60554952, 0.60554952, 0.60554952,
                  0.60554952])
    np.testing.assert_almost_equal(sig_ref, stddevs[0])
    tau_ref = \
        np.array([0.32660535, 0.32660535, 0.32660535, 0.32660535, 0.32660535,
                  0.32660535, 0.32660535, 0.32660535, 0.32660535, 0.32660535,
                  0.32660535])
    np.testing.assert_almost_equal(tau_ref, stddevs[1])
    phi_ref = \
        np.array([0.51411076, 0.51411076, 0.51411076, 0.51411076, 0.51411076,
                  0.51411076, 0.51411076, 0.51411076, 0.51411076, 0.51411076,
                  0.51411076])
    np.testing.assert_almost_equal(phi_ref, stddevs[2])


if __name__ == '__main__':
    test_abrahamsonbhasin2020()
    test_abrahamsonbhasin2020sa1()
    test_abrahamsonbhasin2020pga()
