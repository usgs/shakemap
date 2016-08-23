import os
import sys
import time as time

import numpy as np
import pytest

from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import CampbellBozorgnia2014
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.gsim.campbell_2003 import Campbell2003
from openquake.hazardlib import imt, const

import shakemap.grind.multigmpe as mg
from shakemap.grind.sites import Sites
from shakemap.grind.source import Source
from shakemap.grind.fault import Fault
from shakemap.grind.distance import Distance
from shakemap.utils.timeutils import ShakeDateTime

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_multigmpe():
    # Define gmpes and their weights
    gmpes = [AbrahamsonEtAl2014(), BooreEtAl2014(),
             CampbellBozorgnia2014(), ChiouYoungs2014()]
    wts = [0.25, 0.25, 0.25, 0.25]

    # Make sites instance
    vs30file = os.path.join(shakedir, 'tests/data/Vs30_test.grd')
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    site = Sites.createFromCenter(cx, cy, xspan, yspan, dx, dy,
                                  vs30File=vs30file,
                                  padding=True, resample=False)
    sctx = site.getSitesContext()
    sctx.vs30 = np.reshape(sctx.vs30, (-1,))
    sctx.vs30measured = np.reshape(sctx.vs30measured, (-1,))
    sctx.z1pt0 = np.reshape(sctx.z1pt0, (-1,))

    # Need separate z1pt0 arrays
    sctx.z1pt0cy14 = mg._z1_from_vs30_cy14_cal(sctx.vs30)
    sctx.z1pt0ask14 = mg._z1_from_vs30_ask14_cal(sctx.vs30)
    sctx.z2pt5 = mg._z2p5_from_vs30_cb14_cal(sctx.vs30) / 1000.0

    # Make souce instance
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.])

    flt = Fault.fromTrace(lon0, lat0, lon1, lat1, z, W, dip)
    event = {'lat': 34.1, 'lon': -118.2, 'depth': 1, 'mag': 6,
             'id': '', 'locstring': '', 'rake': 30.3,
             'time': ShakeDateTime.utcfromtimestamp(int(time.time())),
             'timezone': 'UTC'}
    source = Source(event, flt)

    # Make a rupture context
    rupt = source.getRuptureContext(gmpes)

    # Make a distance context
    dctx = Distance.fromSites(gmpes, source, site).getDistanceContext()
    dctx.rhypo = np.reshape(dctx.rhypo, (-1,))
    dctx.rx = np.reshape(dctx.rx, (-1,))
    dctx.rjb = np.reshape(dctx.rjb, (-1,))
    dctx.ry0 = np.reshape(dctx.ry0, (-1,))
    dctx.rrup = np.reshape(dctx.rrup, (-1,))

    # Compute weighted GMPE
    iimt = imt.PGV()
    stddev_types = [const.StdDev.TOTAL]
    mgmpe = mg.MultiGMPE.from_list(gmpes, wts)
    lnmu, lnsd = mgmpe.get_mean_and_stddevs(
        sctx, rupt, dctx, iimt, stddev_types)

    lnmud = np.array(
        [3.44828531,  3.49829605,  3.61749432,  3.64343805,  3.7001028,
         3.7348924,  3.76927164,  3.78659955,  3.82600784,  3.46635007,
         3.53816879,  3.6486898,  3.67058155,  3.72223342,  3.75403094,
         3.79315031,  3.79871491,  3.82093027,  3.54889613,  3.57531437,
         3.64441687,  3.69915981,  3.74491289,  3.78931599,  3.80957828,
         3.80870754,  3.8731021,  3.5927326,  3.60764647,  3.66894024,
         3.72148551,  3.75742965,  3.82164661,  3.86341308,  3.87171115,
         3.79092594,  3.64153758,  3.61835381,  3.68166249,  3.7338161,
         3.82454214,  3.81543928,  3.81507658,  3.80006803,  3.77165695,
         3.65178742,  3.71324776,  3.70389969,  3.77034752,  3.78259432,
         3.78677497,  3.79838465,  3.79050287,  3.75066018,  3.52883328,
         3.67813977,  3.71754876,  3.65520574,  3.69463436,  3.72516445,
         3.7457098,  3.74672185,  3.72615784,  3.44535551,  3.61907294,
         3.58790363,  3.58068716,  3.61177983,  3.64349327,  3.66698468,
         3.67129902,  3.65483002]
    )

    lnsdd = np.array(
       [ 0.63560302,  0.63648101,  0.63610581,  0.6390135 ,  0.64203528,
         0.64624098,  0.64851812,  0.64640406,  0.64384305,  0.6361429 ,
         0.63677975,  0.63715381,  0.64040366,  0.64404005,  0.64782624,
         0.6476325 ,  0.64509458,  0.64297808,  0.63477576,  0.63727968,
         0.63899462,  0.64205578,  0.64604037,  0.64815296,  0.64609948,
         0.64402734,  0.63844724,  0.6343891 ,  0.63806041,  0.64043609,
         0.64406094,  0.64776777,  0.64717195,  0.64297191,  0.64011346,
         0.64110084,  0.63137566,  0.63864151,  0.64163093,  0.64588687,
         0.64714873,  0.64603694,  0.64397734,  0.64217431,  0.63958323,
         0.62883338,  0.63127469,  0.63961477,  0.64097303,  0.6442055 ,
         0.64376449,  0.64273526,  0.64112115,  0.63815862,  0.63575399,
         0.6291859 ,  0.63180644,  0.6394421 ,  0.63946545,  0.63947169,
         0.63935499,  0.63832598,  0.63664816,  0.63595663,  0.62755689,
         0.63523274,  0.63663489,  0.63631586,  0.63616589,  0.63597828,
         0.63542126,  0.63500847])

    np.testing.assert_allclose(lnmu, lnmud)
    np.testing.assert_allclose(lnsd[0], lnsdd)

    # Check for exception due to weights:
    with pytest.raises(Exception) as a:
        wts = [0.25, 0.25, 0.25, 0.25 + 1e-4]
        mgmpe = mg.MultiGMPE.from_list(gmpes, wts)

    # Check exception on GMPE check
    with pytest.raises(Exception) as a:
        wts = [1.0]
        mgmpe = mg.MultiGMPE.from_list(['BA08'], wts)

    # Check exception on tectonic region
    with pytest.raises(Exception) as a:
        gmpes = [BooreEtAl2014(), Campbell2003()]
        wts = [0.5, 0.5]
        mgmpe = mg.MultiGMPE.from_list(gmpes, wts)

    # Check exception on length of gmpe and weight lenghts
    with pytest.raises(Exception) as a:
        gmpes = [BooreEtAl2014(), Campbell2003()]
        wts = [1.0]
        mgmpe = mg.MultiGMPE.from_list(gmpes, wts)
    

    # Check PGV from a GMPE without PGV
    gmpes = [Campbell2003()]
    wts = [1.0]
    mgmpe = mg.MultiGMPE.from_list(gmpes, wts)
    lnmu, lnsd = mgmpe.get_mean_and_stddevs(
        sctx, rupt, dctx, iimt, stddev_types)

    lnmud = np.array(
      [ 3.09152212,  3.1524312 ,  3.20749883,  3.25431585,  3.29035521,
        3.31326677,  3.32116911,  3.31341321,  3.29819842,  3.12252648,
        3.18081138,  3.23208034,  3.27383205,  3.30358765,  3.319195  ,
        3.31916753,  3.30623521,  3.28938984,  3.15235911,  3.20745205,
        3.25429394,  3.29035582,  3.31328548,  3.32119931,  3.31344697,
        3.2982328 ,  3.27982759,  3.17945026,  3.23203088,  3.2738231 ,
        3.30360265,  3.31922869,  3.31921198,  3.30628471,  3.28944133,
        3.26955097,  3.18990634,  3.24351181,  3.28521502,  3.31195497,
        3.32124956,  3.3135073 ,  3.29830033,  3.27989827,  3.25860053,
        3.17942778,  3.23201703,  3.27282524,  3.29888607,  3.3078892 ,
        3.30156745,  3.2884687 ,  3.26964276,  3.24701758,  3.14910673,
        3.19888101,  3.23727522,  3.26163304,  3.2701699 ,  3.2690822 ,
        3.26201491,  3.24919602,  3.23101321,  3.10184816,  3.1475792 ,
        3.18259748,  3.20467529,  3.21444387,  3.21832088,  3.21671138,
        3.20966263,  3.19737325]
    )

    lnsdd = np.array(
       [ 0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518,  0.83458518,  0.83458518,  0.83458518,
         0.83458518,  0.83458518]
    )

    np.testing.assert_allclose(lnmu, lnmud)
    np.testing.assert_allclose(lnsd[0], lnsdd)
