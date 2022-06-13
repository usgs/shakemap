#!/usr/bin/env python
import os
import sys
import numpy as np
import time


from impactutils.rupture.edge_rupture import EdgeRupture, QuadRupture
from impactutils.rupture.origin import Origin
from impactutils.time.ancient_time import HistoricTime

from shakelib.sites import Sites


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, "..", ".."))
sys.path.insert(0, shakedir)


def test_EdgeRupture_vs_QuadRupture():
    # Sites stuff
    sites = Sites.fromCenter(-122.15, 37.15, 1.5, 1.5, 0.01, 0.01)
    sm_dict = sites._GeoDict
    west = sm_dict.xmin
    east = sm_dict.xmax
    south = sm_dict.ymin
    north = sm_dict.ymax
    nx = sm_dict.nx
    ny = sm_dict.ny
    lats = np.linspace(north, south, ny)
    lons = np.linspace(west, east, nx)
    lon, lat = np.meshgrid(lons, lats)
    dep = np.zeros_like(lon)

    # Construct QuadRupture
    xp0 = np.array([-122.0, -122.5])
    yp0 = np.array([37.1, 37.4])
    xp1 = np.array([-121.7, -122.3])
    yp1 = np.array([37.2, 37.2])
    zp = np.array([0, 6])
    widths = np.array([30, 20])
    dips = np.array([30, 40])

    origin = Origin(
        {
            "lat": 33.15,
            "lon": -122.15,
            "depth": 0,
            "mag": 7.2,
            "id": "",
            "netid": "",
            "network": "",
            "locstring": "",
            "time": HistoricTime.utcfromtimestamp(time.time()),
        }
    )
    qrup = QuadRupture.fromTrace(xp0, yp0, xp1, yp1, zp, widths, dips, origin)
    rrup_q, _ = qrup.computeRrup(lon, lat, dep)
    rjb_q, _ = qrup.computeRjb(lon, lat, dep)

    # Construct equivalent EdgeRupture
    toplons = np.array([-122.0, -121.7, -122.5, -122.3])
    toplats = np.array([37.1, 37.2, 37.4, 37.2])
    topdeps = np.array([0, 0, 6, 6])
    botlons = np.array([-121.886864, -121.587568, -122.635467, -122.435338])
    botlats = np.array([36.884527, 36.984246, 37.314035, 37.114261])
    botdeps = np.array([15.0000, 14.9998, 18.8558, 18.8559])
    group_index = [0, 0, 1, 1]

    erup = EdgeRupture.fromArrays(
        toplons, toplats, topdeps, botlons, botlats, botdeps, origin, group_index
    )
    rrup_e, _ = erup.computeRrup(lon, lat, dep)
    rjb_e, _ = erup.computeRjb(lon, lat, dep)

    # Check that QuadRupture and EdgeRupture give the same result
    # (we check the absolute values of QuadRupture elsewhere)
    np.testing.assert_allclose(rrup_e, rrup_q, atol=0.35)
    np.testing.assert_allclose(rjb_e, rjb_q, atol=0.35)


if __name__ == "__main__":
    test_EdgeRupture_vs_QuadRupture()
