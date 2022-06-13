#!/usr/bin/env python

# stdlib imports
import os.path
import sys

# third party imports
from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.berge_thierry_2003 import BergeThierryEtAl2003SIGMA
import numpy as np
import pytest
import time
from impactutils.rupture.distance import Distance, get_distance
from impactutils.rupture.origin import Origin
from impactutils.rupture.point_rupture import PointRupture
from impactutils.rupture.quad_rupture import QuadRupture
from impactutils.time.ancient_time import HistoricTime

# local imports
from shakelib.sites import Sites


do_tests = True

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, "..", ".."))
sys.path.insert(0, shakedir)


def test_exceptions():
    vs30file = os.path.join(homedir, "distance_data/Vs30_test.grd")
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    site = Sites.fromCenter(
        cx, cy, xspan, yspan, dx, dy, vs30File=vs30file, padding=True, resample=False
    )
    # Make souce instance
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.0])

    # Rupture requires an origin even when not used:
    origin = Origin(
        {
            "id": "test",
            "lat": 0,
            "lon": 0,
            "depth": 5.0,
            "mag": 7.0,
            "netid": "",
            "network": "",
            "locstring": "",
            "time": HistoricTime.utcfromtimestamp(int(time.time())),
        }
    )
    rup = QuadRupture.fromTrace(lon0, lat0, lon1, lat1, z, W, dip, origin)

    event = {
        "lat": 34.1,
        "lon": -118.2,
        "depth": 1,
        "mag": 6,
        "id": "",
        "locstring": "",
        "mech": "RS",
        "rake": 90,
        "netid": "",
        "network": "",
        "time": HistoricTime.utcfromtimestamp(int(time.time())),
    }
    origin = Origin(event)

    gmpelist = ["Primate"]
    with pytest.raises(Exception) as e:  # noqa
        Distance.fromSites(gmpelist, origin, site, rup)

    gmpelist = [AbrahamsonEtAl2014()]
    sctx = site.getSitesContext()
    dist_types = ["repi", "rhypo", "rjb", "rrup", "rx", "ry", "ry0", "U", "V"]
    with pytest.raises(Exception) as e:  # noqa
        get_distance(dist_types, sctx.lats, sctx.lons, np.zeros_like(sctx.lons), rup)

    dist_types = ["repi", "rhypo", "rjb", "rrup", "rx", "ry", "ry0", "U", "T"]
    with pytest.raises(Exception) as e:  # noqa
        get_distance(
            dist_types,
            sctx.lats,
            sctx.lons[
                0:4,
            ],
            np.zeros_like(sctx.lons),
            rup,
        )
    # Exception when not a GMPE subclass
    with pytest.raises(Exception) as e:  # noqa
        Distance([None], [-118.2], [34.1], [1], rupture=None)


def test_distance_no_rupture():
    event = {
        "lat": 34.1,
        "lon": -118.2,
        "depth": 1,
        "mag": 6,
        "id": "",
        "locstring": "",
        "mech": "RS",
        "rake": 90,
        "netid": "",
        "network": "",
        "time": HistoricTime.utcfromtimestamp(int(time.time())),
    }
    origin = Origin(event)
    origin.setMechanism("ALL")
    # Make sites instance
    vs30file = os.path.join(homedir, "distance_data/Vs30_test.grd")
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    site = Sites.fromCenter(
        cx, cy, xspan, yspan, dx, dy, vs30File=vs30file, padding=True, resample=False
    )
    # Make souce instance
    #  - Unknown/no tectonic region
    #  - Mech is ALL

    gmpe = AbrahamsonEtAl2014()
    rupture = PointRupture(origin)
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.19350211e00,
                1.01453734e00,
                8.94306248e-01,
                8.51431703e-01,
                8.94306248e-01,
                1.01453734e00,
                1.19350211e00,
            ],
            [
                9.23698454e-01,
                6.97204114e-01,
                5.32067867e-01,
                4.69137288e-01,
                5.32067867e-01,
                6.97204114e-01,
                9.23698454e-01,
            ],
            [
                7.28251778e-01,
                4.44114326e-01,
                2.60572550e-01,
                1.94977658e-01,
                2.60572550e-01,
                4.44114326e-01,
                7.28251778e-01,
            ],
            [
                6.54236979e-01,
                3.39249542e-01,
                1.57170497e-01,
                1.98278110e-05,
                1.57170497e-01,
                3.39249542e-01,
                6.54236979e-01,
            ],
            [
                7.28338531e-01,
                4.44167697e-01,
                2.60583985e-01,
                1.94977658e-01,
                2.60583985e-01,
                4.44167697e-01,
                7.28338531e-01,
            ],
            [
                9.23844143e-01,
                6.97283640e-01,
                5.32091716e-01,
                4.69137288e-01,
                5.32091716e-01,
                6.97283640e-01,
                9.23844143e-01,
            ],
            [
                1.19368104e00,
                1.01462773e00,
                8.94331130e-01,
                8.51431703e-01,
                8.94331130e-01,
                1.01462773e00,
                1.19368104e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                4.0129619,
                3.93137849,
                3.87656959,
                3.85702467,
                3.87656959,
                3.93137849,
                4.0129619,
            ],
            [
                3.88996841,
                3.78671803,
                3.71143853,
                3.68275081,
                3.71143853,
                3.78671803,
                3.88996841,
            ],
            [
                3.80087151,
                3.67134376,
                3.60166506,
                3.58311968,
                3.60166506,
                3.67134376,
                3.80087151,
            ],
            [
                3.7671309,
                3.62390909,
                3.57243062,
                3.53580973,
                3.57243062,
                3.62390909,
                3.7671309,
            ],
            [
                3.80091105,
                3.67136809,
                3.60166829,
                3.58311968,
                3.60166829,
                3.67136809,
                3.80091105,
            ],
            [
                3.89003482,
                3.78675428,
                3.7114494,
                3.68275081,
                3.7114494,
                3.78675428,
                3.89003482,
            ],
            [
                4.01304347,
                3.9314197,
                3.87658093,
                3.85702467,
                3.87658093,
                3.9314197,
                4.01304347,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region unsupported
    #  - Mech is ALL
    origin._tectonic_region = "Volcano"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbt = np.array(
        [
            [
                1.19350211e00,
                1.01453734e00,
                8.94306248e-01,
                8.51431703e-01,
                8.94306248e-01,
                1.01453734e00,
                1.19350211e00,
            ],
            [
                9.23698454e-01,
                6.97204114e-01,
                5.32067867e-01,
                4.69137288e-01,
                5.32067867e-01,
                6.97204114e-01,
                9.23698454e-01,
            ],
            [
                7.28251778e-01,
                4.44114326e-01,
                2.60572550e-01,
                1.94977658e-01,
                2.60572550e-01,
                4.44114326e-01,
                7.28251778e-01,
            ],
            [
                6.54236979e-01,
                3.39249542e-01,
                1.57170497e-01,
                1.98278110e-05,
                1.57170497e-01,
                3.39249542e-01,
                6.54236979e-01,
            ],
            [
                7.28338531e-01,
                4.44167697e-01,
                2.60583985e-01,
                1.94977658e-01,
                2.60583985e-01,
                4.44167697e-01,
                7.28338531e-01,
            ],
            [
                9.23844143e-01,
                6.97283640e-01,
                5.32091716e-01,
                4.69137288e-01,
                5.32091716e-01,
                6.97283640e-01,
                9.23844143e-01,
            ],
            [
                1.19368104e00,
                1.01462773e00,
                8.94331130e-01,
                8.51431703e-01,
                8.94331130e-01,
                1.01462773e00,
                1.19368104e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjbt, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    # Souce instance
    #  - Tectonic region: active
    #  - Mech is ALL

    origin.setMechanism("ALL")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.19350211e00,
                1.01453734e00,
                8.94306248e-01,
                8.51431703e-01,
                8.94306248e-01,
                1.01453734e00,
                1.19350211e00,
            ],
            [
                9.23698454e-01,
                6.97204114e-01,
                5.32067867e-01,
                4.69137288e-01,
                5.32067867e-01,
                6.97204114e-01,
                9.23698454e-01,
            ],
            [
                7.28251778e-01,
                4.44114326e-01,
                2.60572550e-01,
                1.94977658e-01,
                2.60572550e-01,
                4.44114326e-01,
                7.28251778e-01,
            ],
            [
                6.54236979e-01,
                3.39249542e-01,
                1.57170497e-01,
                1.98278110e-05,
                1.57170497e-01,
                3.39249542e-01,
                6.54236979e-01,
            ],
            [
                7.28338531e-01,
                4.44167697e-01,
                2.60583985e-01,
                1.94977658e-01,
                2.60583985e-01,
                4.44167697e-01,
                7.28338531e-01,
            ],
            [
                9.23844143e-01,
                6.97283640e-01,
                5.32091716e-01,
                4.69137288e-01,
                5.32091716e-01,
                6.97283640e-01,
                9.23844143e-01,
            ],
            [
                1.19368104e00,
                1.01462773e00,
                8.94331130e-01,
                8.51431703e-01,
                8.94331130e-01,
                1.01462773e00,
                1.19368104e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                4.0129619,
                3.93137849,
                3.87656959,
                3.85702467,
                3.87656959,
                3.93137849,
                4.0129619,
            ],
            [
                3.88996841,
                3.78671803,
                3.71143853,
                3.68275081,
                3.71143853,
                3.78671803,
                3.88996841,
            ],
            [
                3.80087151,
                3.67134376,
                3.60166506,
                3.58311968,
                3.60166506,
                3.67134376,
                3.80087151,
            ],
            [
                3.7671309,
                3.62390909,
                3.57243062,
                3.53580973,
                3.57243062,
                3.62390909,
                3.7671309,
            ],
            [
                3.80091105,
                3.67136809,
                3.60166829,
                3.58311968,
                3.60166829,
                3.67136809,
                3.80091105,
            ],
            [
                3.89003482,
                3.78675428,
                3.7114494,
                3.68275081,
                3.7114494,
                3.78675428,
                3.89003482,
            ],
            [
                4.01304347,
                3.9314197,
                3.87658093,
                3.85702467,
                3.87658093,
                3.9314197,
                4.01304347,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: active
    #  - Mech is RS

    origin.setMechanism("RS")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                7.76090807e-01,
                6.49225734e-01,
                5.63995966e-01,
                5.33602932e-01,
                5.63995966e-01,
                6.49225734e-01,
                7.76090807e-01,
            ],
            [
                5.84831599e-01,
                4.24273624e-01,
                3.07211355e-01,
                2.62600941e-01,
                3.07211355e-01,
                4.24273624e-01,
                5.84831599e-01,
            ],
            [
                4.46282784e-01,
                2.44862590e-01,
                1.32264468e-01,
                9.99797788e-02,
                1.32264468e-01,
                2.44862590e-01,
                4.46282784e-01,
            ],
            [
                3.93814955e-01,
                1.70987945e-01,
                8.13717378e-02,
                1.03958777e-05,
                8.13717378e-02,
                1.70987945e-01,
                3.93814955e-01,
            ],
            [
                4.46344282e-01,
                2.44900424e-01,
                1.32270097e-01,
                9.99797788e-02,
                1.32270097e-01,
                2.44900424e-01,
                4.46344282e-01,
            ],
            [
                5.84934876e-01,
                4.24329999e-01,
                3.07228262e-01,
                2.62600941e-01,
                3.07228262e-01,
                4.24329999e-01,
                5.84934876e-01,
            ],
            [
                7.76217650e-01,
                6.49289812e-01,
                5.64013604e-01,
                5.33602932e-01,
                5.64013604e-01,
                6.49289812e-01,
                7.76217650e-01,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                3.42235562,
                3.338452,
                3.28208435,
                3.26198358,
                3.28208435,
                3.338452,
                3.42235562,
            ],
            [
                3.29586422,
                3.18967743,
                3.112257,
                3.08275341,
                3.112257,
                3.18967743,
                3.29586422,
            ],
            [
                3.20423343,
                3.07102195,
                2.99912626,
                2.97986242,
                2.99912626,
                3.07102195,
                3.20423343,
            ],
            [
                3.16953325,
                3.02223204,
                2.96875925,
                2.92616469,
                2.96875925,
                3.02223204,
                3.16953325,
            ],
            [
                3.2042741,
                3.07104698,
                2.99912962,
                2.97986242,
                2.99912962,
                3.07104698,
                3.2042741,
            ],
            [
                3.29593253,
                3.18971471,
                3.11226818,
                3.08275341,
                3.11226818,
                3.18971471,
                3.29593253,
            ],
            [
                3.42243951,
                3.33849438,
                3.28209601,
                3.26198358,
                3.28209601,
                3.33849438,
                3.42243951,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: active
    #  - Mech is NM

    origin.setMechanism("NM")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                8.32771820e-01,
                6.96170087e-01,
                6.04399092e-01,
                5.71673449e-01,
                6.04399092e-01,
                6.96170087e-01,
                8.32771820e-01,
            ],
            [
                6.26833822e-01,
                4.53953319e-01,
                3.27906737e-01,
                2.79872556e-01,
                3.27906737e-01,
                4.53953319e-01,
                6.26833822e-01,
            ],
            [
                4.77651641e-01,
                2.60772819e-01,
                1.38685718e-01,
                1.03235484e-01,
                1.38685718e-01,
                2.60772819e-01,
                4.77651641e-01,
            ],
            [
                4.21157003e-01,
                1.81206068e-01,
                8.28029065e-02,
                1.03958777e-05,
                8.28029065e-02,
                1.81206068e-01,
                4.21157003e-01,
            ],
            [
                4.77717859e-01,
                2.60813557e-01,
                1.38691898e-01,
                1.03235484e-01,
                1.38691898e-01,
                2.60813557e-01,
                4.77717859e-01,
            ],
            [
                6.26945025e-01,
                4.54014020e-01,
                3.27924941e-01,
                2.79872556e-01,
                3.27924941e-01,
                4.54014020e-01,
                6.26945025e-01,
            ],
            [
                8.32908398e-01,
                6.96239083e-01,
                6.04418084e-01,
                5.71673449e-01,
                6.04418084e-01,
                6.96239083e-01,
                8.32908398e-01,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                3.3192606,
                3.22072248,
                3.15452316,
                3.13091641,
                3.15452316,
                3.22072248,
                3.3192606,
            ],
            [
                3.17070653,
                3.0459986,
                2.95507447,
                2.92042485,
                2.95507447,
                3.0459986,
                3.17070653,
            ],
            [
                3.06309346,
                2.90664719,
                2.82107391,
                2.79752673,
                2.82107391,
                2.90664719,
                3.06309346,
            ],
            [
                3.02234086,
                2.84931729,
                2.78395476,
                2.73772697,
                2.78395476,
                2.84931729,
                3.02234086,
            ],
            [
                3.06314123,
                2.90667658,
                2.82107802,
                2.79752673,
                2.82107802,
                2.90667658,
                3.06314123,
            ],
            [
                3.17078675,
                3.04604238,
                2.9550876,
                2.92042485,
                2.9550876,
                3.04604238,
                3.17078675,
            ],
            [
                3.31935913,
                3.22077225,
                3.15453686,
                3.13091641,
                3.15453686,
                3.22077225,
                3.31935913,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: active
    #  - Mech is SS

    origin.setMechanism("SS")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.95958776e00,
                1.66988434e00,
                1.47525745e00,
                1.40585328e00,
                1.47525745e00,
                1.66988434e00,
                1.95958776e00,
            ],
            [
                1.52283677e00,
                1.15619376e00,
                8.88875589e-01,
                7.87005240e-01,
                8.88875589e-01,
                1.15619376e00,
                1.52283677e00,
            ],
            [
                1.20645289e00,
                7.46498734e-01,
                4.23057706e-01,
                2.95503135e-01,
                4.23057706e-01,
                7.46498734e-01,
                1.20645289e00,
            ],
            [
                1.08663970e00,
                5.76051478e-01,
                2.21984054e-01,
                1.98278110e-05,
                2.21984054e-01,
                5.76051478e-01,
                1.08663970e00,
            ],
            [
                1.20659332e00,
                7.46585130e-01,
                4.23079943e-01,
                2.95503135e-01,
                4.23079943e-01,
                7.46585130e-01,
                1.20659332e00,
            ],
            [
                1.52307261e00,
                1.15632249e00,
                8.88914196e-01,
                7.87005240e-01,
                8.88914196e-01,
                1.15632249e00,
                1.52307261e00,
            ],
            [
                1.95987741e00,
                1.67003067e00,
                1.47529773e00,
                1.40585328e00,
                1.47529773e00,
                1.67003067e00,
                1.95987741e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                2.54969772,
                2.27038241,
                2.08273439,
                2.01581889,
                2.08273439,
                2.27038241,
                2.54969772,
            ],
            [
                2.12860763,
                1.77511159,
                1.51737884,
                1.41916133,
                1.51737884,
                1.77511159,
                2.12860763,
            ],
            [
                1.82356854,
                1.38010729,
                1.08693739,
                0.97911408,
                1.08693739,
                1.38010729,
                1.82356854,
            ],
            [
                1.70805158,
                1.21626476,
                0.91696757,
                0.78911491,
                0.91696757,
                1.21626476,
                1.70805158,
            ],
            [
                1.82370394,
                1.38019059,
                1.08695619,
                0.97911408,
                1.08695619,
                1.38019059,
                1.82370394,
            ],
            [
                2.12883501,
                1.77523571,
                1.51741606,
                1.41916133,
                1.51741606,
                1.77523571,
                2.12883501,
            ],
            [
                2.54997699,
                2.27052349,
                2.08277323,
                2.01581889,
                2.08277323,
                2.27052349,
                2.54997699,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is all

    origin.setMechanism("ALL")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.49285078e00,
                1.26359361e00,
                1.10957536e00,
                1.05465228e00,
                1.10957536e00,
                1.26359361e00,
                1.49285078e00,
            ],
            [
                1.14722732e00,
                8.57083889e-01,
                6.45541307e-01,
                5.64926073e-01,
                6.45541307e-01,
                8.57083889e-01,
                1.14722732e00,
            ],
            [
                8.96856520e-01,
                5.32871196e-01,
                2.99662245e-01,
                2.17185537e-01,
                2.99662245e-01,
                5.32871196e-01,
                8.96856520e-01,
            ],
            [
                8.02042196e-01,
                3.98587924e-01,
                1.69648145e-01,
                1.98278110e-05,
                1.69648145e-01,
                3.98587924e-01,
                8.02042196e-01,
            ],
            [
                8.96967653e-01,
                5.32939565e-01,
                2.99676623e-01,
                2.17185537e-01,
                2.99676623e-01,
                5.32939565e-01,
                8.96967653e-01,
            ],
            [
                1.14741395e00,
                8.57185764e-01,
                6.45571858e-01,
                5.64926073e-01,
                6.45571858e-01,
                8.57185764e-01,
                1.14741395e00,
            ],
            [
                1.49308000e00,
                1.26370940e00,
                1.10960724e00,
                1.05465228e00,
                1.10960724e00,
                1.26370940e00,
                1.49308000e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                4.17967552,
                4.07332411,
                4.00187571,
                3.97639713,
                4.00187571,
                4.07332411,
                4.17967552,
            ],
            [
                4.01934229,
                3.88474601,
                3.78661232,
                3.74921526,
                3.78661232,
                3.88474601,
                4.01934229,
            ],
            [
                3.90319636,
                3.73434515,
                3.64558217,
                3.62308648,
                3.64558217,
                3.73434515,
                3.90319636,
            ],
            [
                3.85921241,
                3.67256434,
                3.61012056,
                3.57133422,
                3.61012056,
                3.67256434,
                3.85921241,
            ],
            [
                3.90324792,
                3.73437686,
                3.64558609,
                3.62308648,
                3.64558609,
                3.73437686,
                3.90324792,
            ],
            [
                4.01942887,
                3.88479327,
                3.7866265,
                3.74921526,
                3.7866265,
                3.88479327,
                4.01942887,
            ],
            [
                4.17978186,
                4.07337783,
                4.0018905,
                3.97639713,
                4.0018905,
                4.07337783,
                4.17978186,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is RS

    origin.setMechanism("RS")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.11052523e00,
                9.25877479e-01,
                8.01828481e-01,
                7.57592465e-01,
                8.01828481e-01,
                9.25877479e-01,
                1.11052523e00,
            ],
            [
                8.32154030e-01,
                5.98467416e-01,
                4.28087307e-01,
                3.63158382e-01,
                4.28087307e-01,
                5.98467416e-01,
                8.32154030e-01,
            ],
            [
                6.30500991e-01,
                3.37340822e-01,
                1.69925286e-01,
                1.20068361e-01,
                1.69925286e-01,
                3.37340822e-01,
                6.30500991e-01,
            ],
            [
                5.54135870e-01,
                2.29725567e-01,
                9.13321474e-02,
                1.03958777e-05,
                9.13321474e-02,
                2.29725567e-01,
                5.54135870e-01,
            ],
            [
                6.30590499e-01,
                3.37395888e-01,
                1.69933978e-01,
                1.20068361e-01,
                1.69933978e-01,
                3.37395888e-01,
                6.30590499e-01,
            ],
            [
                8.32304345e-01,
                5.98549467e-01,
                4.28111914e-01,
                3.63158382e-01,
                4.28111914e-01,
                5.98549467e-01,
                8.32304345e-01,
            ],
            [
                1.11070985e00,
                9.25970743e-01,
                8.01854154e-01,
                7.57592465e-01,
                8.01854154e-01,
                9.25970743e-01,
                1.11070985e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                3.4885951,
                3.37216961,
                3.29395331,
                3.26606128,
                3.29395331,
                3.37216961,
                3.4885951,
            ],
            [
                3.3130744,
                3.16572856,
                3.05829921,
                3.01735974,
                3.05829921,
                3.16572856,
                3.3130744,
            ],
            [
                3.18592661,
                3.00108105,
                2.90341742,
                2.87839095,
                2.90341742,
                3.00108105,
                3.18592661,
            ],
            [
                3.1377763,
                2.9334351,
                2.86396637,
                2.81798622,
                2.86396637,
                2.9334351,
                3.1377763,
            ],
            [
                3.18598305,
                3.00111577,
                2.90342178,
                2.87839095,
                2.90342178,
                3.00111577,
                3.18598305,
            ],
            [
                3.31316918,
                3.16578029,
                3.05831472,
                3.01735974,
                3.05831472,
                3.16578029,
                3.31316918,
            ],
            [
                3.48871151,
                3.37222842,
                3.29396949,
                3.26606128,
                3.29396949,
                3.37222842,
                3.48871151,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is NM

    origin.setMechanism("NM")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.12678662e00,
                9.39133949e-01,
                8.13066202e-01,
                7.68110298e-01,
                8.13066202e-01,
                9.39133949e-01,
                1.12678662e00,
            ],
            [
                8.43885262e-01,
                6.06395679e-01,
                4.33242838e-01,
                3.67257274e-01,
                4.33242838e-01,
                6.06395679e-01,
                8.43885262e-01,
            ],
            [
                6.38950562e-01,
                3.41019564e-01,
                1.70913434e-01,
                1.20272659e-01,
                1.70913434e-01,
                3.41019564e-01,
                6.38950562e-01,
            ],
            [
                5.61342691e-01,
                2.31653894e-01,
                9.10846554e-02,
                1.03958777e-05,
                9.10846554e-02,
                2.31653894e-01,
                5.61342691e-01,
            ],
            [
                6.39041527e-01,
                3.41075526e-01,
                1.70922263e-01,
                1.20272659e-01,
                1.70922263e-01,
                3.41075526e-01,
                6.39041527e-01,
            ],
            [
                8.44038024e-01,
                6.06479066e-01,
                4.33267846e-01,
                3.67257274e-01,
                4.33267846e-01,
                6.06479066e-01,
                8.44038024e-01,
            ],
            [
                1.12697424e00,
                9.39228730e-01,
                8.13092292e-01,
                7.68110298e-01,
                8.13092292e-01,
                9.39228730e-01,
                1.12697424e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                3.42781739,
                3.30181908,
                3.21717161,
                3.18698623,
                3.21717161,
                3.30181908,
                3.42781739,
            ],
            [
                3.23786489,
                3.07840387,
                2.96214139,
                2.91783576,
                2.96214139,
                3.07840387,
                3.23786489,
            ],
            [
                3.10026266,
                2.9002186,
                2.79362772,
                2.76581535,
                2.79362772,
                2.9002186,
                3.10026266,
            ],
            [
                3.0481533,
                2.82698693,
                2.74978504,
                2.70136713,
                2.74978504,
                2.82698693,
                3.0481533,
            ],
            [
                3.10032374,
                2.90025617,
                2.79363257,
                2.76581535,
                2.79363257,
                2.90025617,
                3.10032374,
            ],
            [
                3.23796746,
                3.07845986,
                2.96215818,
                2.91783576,
                2.96215818,
                3.07845986,
                3.23796746,
            ],
            [
                3.42794337,
                3.30188272,
                3.21718913,
                3.18698623,
                3.21718913,
                3.30188272,
                3.42794337,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))

    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is SS

    origin.setMechanism("SS")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjb = np.array(
        [
            [
                1.80104893e00,
                1.52092305e00,
                1.33273049e00,
                1.26562081e00,
                1.33273049e00,
                1.52092305e00,
                1.80104893e00,
            ],
            [
                1.37873685e00,
                1.02421498e00,
                7.65734302e-01,
                6.67231768e-01,
                7.65734302e-01,
                1.02421498e00,
                1.37873685e00,
            ],
            [
                1.07281256e00,
                6.28064399e-01,
                3.42919369e-01,
                2.41987662e-01,
                3.42919369e-01,
                6.28064399e-01,
                1.07281256e00,
            ],
            [
                9.56960370e-01,
                4.63980672e-01,
                1.83813296e-01,
                1.98278110e-05,
                1.83813296e-01,
                4.63980672e-01,
                9.56960370e-01,
            ],
            [
                1.07294835e00,
                6.28147939e-01,
                3.42936965e-01,
                2.41987662e-01,
                3.42936965e-01,
                6.28147939e-01,
                1.07294835e00,
            ],
            [
                1.37896489e00,
                1.02433946e00,
                7.65771633e-01,
                6.67231768e-01,
                7.65771633e-01,
                1.02433946e00,
                1.37896489e00,
            ],
            [
                1.80132901e00,
                1.52106454e00,
                1.33276944e00,
                1.26562081e00,
                1.33276944e00,
                1.52106454e00,
                1.80132901e00,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rjb))

    rrup = np.array(
        [
            [
                2.85894272,
                2.62140075,
                2.46181667,
                2.4049088,
                2.46181667,
                2.62140075,
                2.85894272,
            ],
            [
                2.50082927,
                2.20020077,
                1.98101356,
                1.89748509,
                1.98101356,
                2.20020077,
                2.50082927,
            ],
            [
                2.24141069,
                1.86427183,
                1.65402932,
                1.59405522,
                1.65402932,
                1.86427183,
                2.24141069,
            ],
            [
                2.14317001,
                1.72596453,
                1.55948774,
                1.48557451,
                1.55948774,
                1.72596453,
                2.14317001,
            ],
            [
                2.24152584,
                1.86434267,
                1.65403978,
                1.59405522,
                1.65403978,
                1.86434267,
                2.24152584,
            ],
            [
                2.50102265,
                2.20030633,
                1.98104522,
                1.89748509,
                1.98104522,
                2.20030633,
                2.50102265,
            ],
            [
                2.85918022,
                2.62152073,
                2.46184969,
                2.4049088,
                2.46184969,
                2.62152073,
                2.85918022,
            ],
        ]
    )

    if do_tests is True:
        np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)
    else:
        print(repr(dctx.rrup))


def test_distance_from_sites_origin():
    # Make sites instance
    vs30file = os.path.join(homedir, "distance_data/Vs30_test.grd")
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    site = Sites.fromCenter(
        cx, cy, xspan, yspan, dx, dy, vs30File=vs30file, padding=True, resample=False
    )
    # Make souce instance
    lat0 = np.array([34.1])
    lon0 = np.array([-118.2])
    lat1 = np.array([34.2])
    lon1 = np.array([-118.15])
    z = np.array([1.0])
    W = np.array([3.0])
    dip = np.array([30.0])

    event = {
        "lat": 34.1,
        "lon": -118.2,
        "depth": 1,
        "mag": 6,
        "id": "",
        "locstring": "",
        "mech": "ALL",
        "netid": "",
        "network": "",
        "time": HistoricTime.utcfromtimestamp(int(time.time())),
    }
    origin = Origin(event)

    rup = QuadRupture.fromTrace(lon0, lat0, lon1, lat1, z, W, dip, origin)
    gmpelist = [AbrahamsonEtAl2014(), BergeThierryEtAl2003SIGMA()]
    dists = Distance.fromSites(gmpelist, site, rup)
    dctx = dists.getDistanceContext()

    rhypo = np.array(
        [
            [
                3.74498133,
                3.32896405,
                3.05225679,
                2.95426722,
                3.05225679,
                3.32896405,
                3.74498133,
            ],
            [
                3.11965436,
                2.60558436,
                2.24124201,
                2.10583262,
                2.24124201,
                2.60558436,
                3.11965436,
            ],
            [
                2.67523213,
                2.05265767,
                1.564393,
                1.36331682,
                1.564393,
                2.05265767,
                2.67523213,
            ],
            [
                2.50973226,
                1.83166664,
                1.26045653,
                1.0,
                1.26045653,
                1.83166664,
                2.50973226,
            ],
            [
                2.67542717,
                2.05277065,
                1.56443006,
                1.36331682,
                1.56443006,
                2.05277065,
                2.67542717,
            ],
            [
                3.11998886,
                2.60576236,
                2.24129374,
                2.10583262,
                2.24129374,
                2.60576236,
                3.11998886,
            ],
            [
                3.74539929,
                3.32917303,
                3.05231378,
                2.95426722,
                3.05231378,
                3.32917303,
                3.74539929,
            ],
        ]
    )
    np.testing.assert_allclose(rhypo, dctx.rhypo, rtol=0, atol=0.01)

    rx = np.array(
        [
            [
                -3.18894050e00,
                -2.48001769e00,
                -1.77111874e00,
                -1.06224366e00,
                -3.53392480e-01,
                3.55434794e-01,
                1.06423815e00,
            ],
            [
                -2.83506890e00,
                -2.12607622e00,
                -1.41710740e00,
                -7.08162466e-01,
                7.58576362e-04,
                7.09655709e-01,
                1.41852892e00,
            ],
            [
                -2.48119723e00,
                -1.77213470e00,
                -1.06309603e00,
                -3.54081243e-01,
                3.54909645e-01,
                1.06387662e00,
                1.77281967e00,
            ],
            [
                -2.12732550e00,
                -1.41819312e00,
                -7.09084619e-01,
                2.56774082e-12,
                7.09060719e-01,
                1.41809752e00,
                2.12711040e00,
            ],
            [
                -1.77345370e00,
                -1.06425151e00,
                -3.55073182e-01,
                3.54081255e-01,
                1.06321179e00,
                1.77231841e00,
                2.48140110e00,
            ],
            [
                -1.41958186e00,
                -7.10309855e-01,
                -1.06172493e-03,
                7.08162516e-01,
                1.41736285e00,
                2.12653927e00,
                2.83569175e00,
            ],
            [
                -1.06570997e00,
                -3.56368176e-01,
                3.52949744e-01,
                1.06224377e00,
                1.77151390e00,
                2.48076010e00,
                3.18998236e00,
            ],
        ]
    )

    np.testing.assert_allclose(rx, dctx.rx, rtol=0, atol=0.01)

    rjb = np.array(
        [
            [
                3.19372137e00,
                2.48373511e00,
                1.77377308e00,
                1.06383562e00,
                3.53925643e-01,
                2.25816823e-03,
                2.45009861e-03,
            ],
            [
                2.83931844e00,
                2.12926243e00,
                1.41923064e00,
                7.09223517e-01,
                1.57594916e-03,
                1.86044244e-03,
                2.05239165e-03,
            ],
            [
                2.48510934e00,
                1.77479025e00,
                1.06468863e00,
                3.54611655e-01,
                1.04375185e-03,
                1.32827303e-03,
                1.52024106e-03,
            ],
            [
                2.30690967e00,
                1.53793979e00,
                7.68969896e-01,
                5.88918451e-12,
                3.77111295e-04,
                6.61660373e-04,
                8.53647223e-04,
            ],
            [
                2.48531877e00,
                1.79442084e00,
                1.20242597e00,
                8.54793253e-01,
                5.62052963e-01,
                2.69254693e-01,
                5.26105100e-05,
            ],
            [
                2.95646628e00,
                2.40489915e00,
                2.00231070e00,
                1.70958533e00,
                1.41681634e00,
                1.12398937e00,
                8.63761551e-01,
            ],
            [
                3.60741953e00,
                3.17112489e00,
                2.85711592e00,
                2.56437623e00,
                2.27157856e00,
                1.97872291e00,
                1.78518260e00,
            ],
        ]
    )

    np.testing.assert_allclose(rjb, dctx.rjb, rtol=0, atol=0.01)

    ry0 = np.array(
        [
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                2.29490054e-02,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                8.79341006e-01,
                5.86285236e-01,
                2.93171565e-01,
                6.21003581e-12,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                1.73573289e00,
                1.44264826e00,
                1.14950573e00,
                8.56305300e-01,
                5.63046975e-01,
                2.69730762e-01,
                0.00000000e00,
            ],
            [
                2.59212463e00,
                2.29901116e00,
                2.00583977e00,
                1.71261048e00,
                1.41932329e00,
                1.12597821e00,
                8.32575235e-01,
            ],
            [
                3.44851622e00,
                3.15537391e00,
                2.86217367e00,
                2.56891553e00,
                2.27559947e00,
                1.98222553e00,
                1.68879368e00,
            ],
        ]
    )

    np.testing.assert_allclose(ry0, dctx.ry0, rtol=0, atol=0.01)

    rrup = np.array(
        [
            [
                3.34678672,
                2.67788811,
                2.03697073,
                1.46129187,
                1.06271102,
                1.06352692,
                1.40073832,
            ],
            [
                3.01030105,
                2.3526499,
                1.73673635,
                1.22706347,
                1.00157564,
                1.22283363,
                1.57764099,
            ],
            [
                2.67858182,
                2.03712377,
                1.46095502,
                1.06170931,
                1.06220616,
                1.39958479,
                1.75442695,
            ],
            [2.51415965, 1.8343632, 1.26143652, 1.0, 1.2212501, 1.57621925, 1.9310962],
            [
                2.67877609,
                2.05412785,
                1.56384179,
                1.3617346,
                1.50608502,
                1.77308319,
                2.10764873,
            ],
            [
                3.12078859,
                2.6043486,
                2.23799413,
                2.09885629,
                2.11696797,
                2.23191013,
                2.4299612,
            ],
            [
                3.74318473,
                3.32482368,
                3.04635272,
                2.9183523,
                2.86659485,
                2.88815116,
                2.98141559,
            ],
        ]
    )

    np.testing.assert_allclose(rrup, dctx.rrup, rtol=0, atol=0.01)


if __name__ == "__main__":
    test_exceptions()
    test_distance_no_rupture()
    test_distance_from_sites_origin()
