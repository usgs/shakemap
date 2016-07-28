import os
import sys
import time as time

import numpy as np

from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import CampbellBozorgnia2014
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib import imt, const

import shakemap.gmpe.multigmpe as mg
from shakemap.grind.sites import Sites
from shakemap.grind.source import Source
from shakemap.grind.fault import Fault
from shakemap.grind.distance import Distance
from shakemap.utils.timeutils import ShakeDateTime

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..'))
# put this at the front of the system path, ignoring any installed mapio stuff
sys.path.insert(0, shakedir)


def test_multigmpe():
    # Define gmpes and their weights
    gmpes = [AbrahamsonEtAl2014(), BooreEtAl2014(),
             CampbellBozorgnia2014(), ChiouYoungs2014()]
    wts = [0.25, 0.25, 0.25, 0.25]

    # Make sites instance
    vs30file = os.path.join(shakedir, 'data/Vs30_test.grd')
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    site = Sites.createFromCenter(cx, cy, xspan, yspan, dx, dy,
                                  vs30File = vs30file,
                                  padding = True, resample = False)
    sctx = site.getSitesContext()
    sctx.vs30 = np.reshape(sctx.vs30, (-1,))
    sctx.vs30measured = np.reshape(sctx.vs30measured, (-1,))
    sctx.z1pt0 = np.reshape(sctx.z1pt0, (-1,))
    
    # Need separate z1pt0 arrays
    sctx.z1pt0cy14 = mg.z1_from_vs30_cy14_cal(sctx.vs30)
    sctx.z1pt0ask14 = mg.z1_from_vs30_ask14_cal(sctx.vs30)
    sctx.z2pt5 = mg.z2p5_from_vs30_cb14_cal(sctx.vs30) / 1000.0
    
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
      [ 3.44828531,  3.49829605,  3.61749432,  3.64343805,  3.7001028 ,
        3.7348924 ,  3.76927164,  3.78659955,  3.82600784,  3.46635007,
        3.53816879,  3.6486898 ,  3.67058155,  3.72223342,  3.75403094,
        3.79315031,  3.79871491,  3.82093027,  3.54889613,  3.57531437,
        3.64441687,  3.69915981,  3.74491289,  3.78931599,  3.80957828,
        3.80870754,  3.8731021 ,  3.5927326 ,  3.60764647,  3.66894024,
        3.72148551,  3.75742965,  3.82164661,  3.86341308,  3.87171115,
        3.79092594,  3.64153758,  3.61835381,  3.68166249,  3.7338161 ,
        3.82454214,  3.81543928,  3.81507658,  3.80006803,  3.77165695,
        3.65178742,  3.71324776,  3.70389969,  3.77034752,  3.78259432,
        3.78677497,  3.79838465,  3.79050287,  3.75066018,  3.52883328,
        3.67813977,  3.71754876,  3.65520574,  3.69463436,  3.72516445,
        3.7457098 ,  3.74672185,  3.72615784,  3.44535551,  3.61907294,
        3.58790363,  3.58068716,  3.61177983,  3.64349327,  3.66698468,
        3.67129902,  3.65483002]
                    )

    lnsdd = np.array(
      [[ 0.64164575,  0.64251245,  0.64214907,  0.64502258,  0.64801656,
         0.65218297,  0.65443907,  0.65234344,  0.64981076,  0.64217808,
         0.64280892,  0.64318596,  0.64639957,  0.65000285,  0.65375482,
         0.65356255,  0.65104703,  0.64895333,  0.64083309,  0.6433045 ,
         0.64500401,  0.64803674,  0.65198555,  0.65407996,  0.65204472,
         0.64999085,  0.64449041,  0.64045403,  0.64407781,  0.64643149,
         0.65002344,  0.65369733,  0.65311033,  0.64896129,  0.64613641,
         0.64709103,  0.63748584,  0.6446528 ,  0.64761461,  0.65183228,
         0.65309205,  0.65198359,  0.64994219,  0.64815471,  0.64558726,
         0.63498016,  0.63739414,  0.64562305,  0.64697806,  0.65017003,
         0.64973003,  0.6487104 ,  0.64711108,  0.64417615,  0.64179696,
         0.63533018,  0.63791922,  0.64544734,  0.64546954,  0.64547626,
         0.64536096,  0.64434165,  0.64268083,  0.64199403,  0.63371941,
         0.64128579,  0.64266503,  0.64234883,  0.64220083,  0.64201544,
         0.64146392,  0.64105492]]
                    )

    np.testing.assert_allclose(lnmu, lnmud)
    np.testing.assert_allclose(lnsd, lnsdd)
