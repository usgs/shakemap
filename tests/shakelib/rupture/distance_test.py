#!/usr/bin/env python

# stdlib imports
import os
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

    datadir = os.path.join(homedir, "data")

    if do_tests is False:
        os.makedirs(datadir, exist_ok=True)

    # TEST #1
    # Make souce instance
    #  - Unknown/no tectonic region
    #  - Mech is ALL

    gmpe = AbrahamsonEtAl2014()
    rupture = PointRupture(origin)
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test1_rjb.npy")
    rrupfile = os.path.join(datadir, "test1_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #2
    # Souce instance
    #  - Tectonic region unsupported
    #  - Mech is ALL
    origin._tectonic_region = "Volcano"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test2_rjb.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)

    # TEST #3
    # Souce instance
    #  - Tectonic region: active
    #  - Mech is ALL

    origin.setMechanism("ALL")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test3_rjb.npy")
    rrupfile = os.path.join(datadir, "test3_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #4
    # Souce instance
    #  - Tectonic region: active
    #  - Mech is RS

    origin.setMechanism("RS")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test4_rjb.npy")
    rrupfile = os.path.join(datadir, "test4_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #5
    # Souce instance
    #  - Tectonic region: active
    #  - Mech is NM

    origin.setMechanism("NM")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test5_rjb.npy")
    rrupfile = os.path.join(datadir, "test5_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #6
    # Souce instance
    #  - Tectonic region: active
    #  - Mech is SS

    origin.setMechanism("SS")
    origin._tectonic_region = "Active Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test6_rjb.npy")
    rrupfile = os.path.join(datadir, "test6_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #7
    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is all

    origin.setMechanism("ALL")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test7_rjb.npy")
    rrupfile = os.path.join(datadir, "test7_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #8
    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is RS

    origin.setMechanism("RS")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test8_rjb.npy")
    rrupfile = os.path.join(datadir, "test8_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #9
    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is NM

    origin.setMechanism("NM")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test9_rjb.npy")
    rrupfile = os.path.join(datadir, "test9_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)

    # TEST #10
    # Souce instance
    #  - Tectonic region: stable
    #  - Mech is SS

    origin.setMechanism("SS")
    origin._tectonic_region = "Stable Shallow Crust"
    dists = Distance.fromSites(gmpe, site, rupture)
    dctx = dists.getDistanceContext()

    rjbfile = os.path.join(datadir, "test10_rjb.npy")
    rrupfile = os.path.join(datadir, "test10_rrup.npy")
    if do_tests is True:
        np.testing.assert_allclose(np.load(rjbfile), dctx.rjb, rtol=0, atol=0.01)
        np.testing.assert_allclose(np.load(rrupfile), dctx.rrup, rtol=0, atol=0.01)
    else:
        np.save(rjbfile, dctx.rjb, allow_pickle=False)
        np.save(rrupfile, dctx.rrup, allow_pickle=False)


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
