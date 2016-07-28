
import numpy as np

import shakemap.grind.gmpe2shakemap as g2s
from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib.gsim.frankel_1996 import FrankelEtAl1996MblgAB1987NSHMP2008
from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008

# Inputs
sm_amps_in = np.array([0.01, 0.03, 0.1, 0.3]) * 100.
gmpe_amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.3]))
sigma_in_log = np.array([0.2, 0.3, 0.4, 0.5])
sigma_in_lin = np.array([1.23372505, 1.36350428, 1.50840426, 1.66936801])


def test_pga():
    gmpe = FrankelEtAl1996MblgAB1987NSHMP2008()
    imt = PGA()
    gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        gmpe_amps_out,
        np.array([-4.70048037, -3.60186808, -2.39789527, -1.29928298]))

    sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        sm_amps_out,
        np.array([1.1, 3.3, 11., 33.]))

    sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in_lin, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([0.2, 0.3, 0.4, 0.5]))

    sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in_log, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([1.23372505, 1.36350428, 1.50840426, 1.66936801]))


def test_pgv():
    gmpe = AbrahamsonEtAl2014()
    imt = PGV()
    gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        gmpe_amps_out,
        np.array([-0.09531018, 1.00330211, 2.20727491, 3.3058872]))

    sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        sm_amps_out,
        np.array([0.011, 0.033, 0.11, 0.33]))

    sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in_lin, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([0.2009975, 0.3006659, 0.4004997, 0.5003998]))

    sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in_log, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([1.23250054, 1.36258854, 1.50764041, 1.66869003]))


def test_sa03():
    gmpe = CampbellBozorgnia2008()
    imt = SA(0.3)
    gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        gmpe_amps_out,
        np.array([-4.7374321, -3.6388198, -2.434847, -1.3362347]))

    sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
    np.testing.assert_almost_equal(
        sm_amps_out,
        np.array([1.1414072, 3.4242217, 11.4140722, 34.2422167]))

    sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in_lin, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([0.2023046, 0.3015413, 0.4011573, 0.5009263]))

    sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in_log, gmpe, imt)
    np.testing.assert_almost_equal(
        sigma_out,
        np.array([1.2308798, 1.3613796, 1.5066329, 1.6677961]))
