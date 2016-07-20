
import numpy as np

import shakemap.grind.gmpe2shakemap as g2s
from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib.gsim.frankel_1996 import FrankelEtAl1996MblgAB1987NSHMP2008
from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2008 import CampbellBozorgnia2008

# Inputs
sm_amps_in = np.array([0.01, 0.03, 0.1, 0.3])*100.
gmpe_amps_in = np.log(np.array([0.01, 0.03, 0.1, 0.3]))
sigma_in = np.array([0.2, 0.3, 0.4, 0.5])

# PGA
gmpe = FrankelEtAl1996MblgAB1987NSHMP2008()
imt  = PGA()
gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    gmpe_amps_out,
    np.array([-4.70048037, -3.60186808, -2.39789527, -1.29928298]))

sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    sm_amps_out,
    np.array([1.1,   3.3,  11. ,  33.]))

sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.57711868,  1.17934718,  0.8969858 ,  0.67778574]))


sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.23372505,  1.36350428,  1.50840426,  1.66936801]))


# PGV
gmpe = AbrahamsonEtAl2014()
imt  = PGV()
gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    gmpe_amps_out,
    np.array([-0.09531018,  1.00330211,  2.20727491,  3.3058872]))

sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    sm_amps_out,
    np.array([ 0.011,  0.033,  0.11 ,  0.33]))

sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.57724549,  1.17951676,  0.89720874,  0.67808076]))


sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.23250054,  1.36258854,  1.50764041,  1.66869003]))

# SA
gmpe = CampbellBozorgnia2008()
imt  = SA()
gmpe_amps_out = g2s.ampShakeMapToGMPE(sm_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    gmpe_amps_out,
    np.array([-0.09531018,  1.00330211,  2.20727491,  3.3058872]))

sm_amps_out = g2s.ampGmpeToShakeMap(gmpe_amps_in, gmpe, imt)
np.testing.assert_almost_equal(
    sm_amps_out,
    np.array([ 0.011,  0.033,  0.11 ,  0.33]))

sigma_out = g2s.sigmaShakeMapToGMPE(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.57724549,  1.17951676,  0.89720874,  0.67808076]))


sigma_out = g2s.sigmaGmpeToShakeMap(sigma_in, gmpe, imt)
np.testing.assert_almost_equal(
    sigma_out,
    np.array([1.23250054,  1.36258854,  1.50764041,  1.66869003]))





