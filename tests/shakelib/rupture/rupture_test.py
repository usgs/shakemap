#!/usr/bin/env python

# stdlib imports
import os
import os.path
import time

# third party
import numpy as np
from openquake.hazardlib.geo.geodetic import azimuth
from mapio.geodict import GeoDict
import matplotlib.pyplot as plt

from impactutils.rupture.origin import Origin
from impactutils.rupture.quad_rupture import QuadRupture
from impactutils.time.ancient_time import HistoricTime

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?

do_tests = True


def test_rupture_depth(interactive=False):
    DIP = 17.0
    WIDTH = 20.0
    GRIDRES = 0.1

    names = ["single", "double", "triple", "concave", "concave_simple", "ANrvSA"]
    means = [
        3.1554422780092461,
        2.9224454569459781,
        3.0381968625073563,
        2.0522694624400271,
        2.4805390352818755,
        2.8740121776209673,
    ]
    stds = [
        2.1895293825074575,
        2.0506459673526174,
        2.0244588429154402,
        2.0112565876976416,
        2.1599789955270019,
        1.6156220309120068,
    ]
    xp0list = [
        np.array([118.3]),
        np.array([10.1, 10.1]),
        np.array([10.1, 10.1, 10.3]),
        np.array([10.9, 10.5, 10.9]),
        np.array([10.9, 10.6]),
        np.array(
            [
                -76.483,
                -76.626,
                -76.757,
                -76.99,
                -77.024,
                -76.925,
                -76.65,
                -76.321,
                -75.997,
                -75.958,
            ]
        ),
    ]
    xp1list = [
        np.array([118.3]),
        np.array([10.1, 10.3]),
        np.array([10.1, 10.3, 10.1]),
        np.array([10.5, 10.9, 11.3]),
        np.array([10.6, 10.9]),
        np.array(
            [
                -76.626,
                -76.757,
                -76.99,
                -77.024,
                -76.925,
                -76.65,
                -76.321,
                -75.997,
                -75.958,
                -76.006,
            ]
        ),
    ]
    yp0list = [
        np.array([34.2]),
        np.array([34.2, 34.5]),
        np.array([34.2, 34.5, 34.8]),
        np.array([34.2, 34.5, 34.8]),
        np.array([35.1, 35.2]),
        np.array(
            [
                -52.068,
                -51.377,
                -50.729,
                -49.845,
                -49.192,
                -48.507,
                -47.875,
                -47.478,
                -47.08,
                -46.422,
            ]
        ),
    ]
    yp1list = [
        np.array([34.5]),
        np.array([34.5, 34.8]),
        np.array([34.5, 34.8, 35.1]),
        np.array([34.5, 34.8, 34.6]),
        np.array([35.2, 35.4]),
        np.array(
            [
                -51.377,
                -50.729,
                -49.845,
                -49.192,
                -48.507,
                -47.875,
                -47.478,
                -47.08,
                -46.422,
                -45.659,
            ]
        ),
    ]

    for i in range(0, len(xp0list)):
        xp0 = xp0list[i]
        xp1 = xp1list[i]
        yp0 = yp0list[i]
        yp1 = yp1list[i]
        name = names[i]
        mean_value = means[i]
        std_value = stds[i]

        zp = np.zeros(xp0.shape)
        strike = azimuth(xp0[0], yp0[0], xp1[-1], yp1[-1])
        widths = np.ones(xp0.shape) * WIDTH
        dips = np.ones(xp0.shape) * DIP
        strike = [strike]

    origin = Origin(
        {
            "id": "test",
            "lon": 0,
            "lat": 0,
            "depth": 5.0,
            "mag": 7.0,
            "netid": "us",
            "network": "",
            "locstring": "",
            "time": HistoricTime.utcfromtimestamp(time.time()),
        }
    )

    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths, dips, origin, strike=strike
    )

    # make a grid of points over both quads, ask for depths
    ymin = np.nanmin(rupture.lats)
    ymax = np.nanmax(rupture.lats)
    xmin = np.nanmin(rupture.lons)
    xmax = np.nanmax(rupture.lons)

    xmin = np.floor(xmin * (1 / GRIDRES)) / (1 / GRIDRES)
    xmax = np.ceil(xmax * (1 / GRIDRES)) / (1 / GRIDRES)
    ymin = np.floor(ymin * (1 / GRIDRES)) / (1 / GRIDRES)
    ymax = np.ceil(ymax * (1 / GRIDRES)) / (1 / GRIDRES)
    geodict = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax, GRIDRES, GRIDRES)
    nx = geodict.nx
    ny = geodict.ny
    depths = np.zeros((ny, nx))
    for row in range(0, ny):
        for col in range(0, nx):
            lat, lon = geodict.getLatLon(row, col)
            depth = rupture.getDepthAtPoint(lat, lon)
            depths[row, col] = depth

    np.testing.assert_almost_equal(np.nanmean(depths), mean_value)
    np.testing.assert_almost_equal(np.nanstd(depths), std_value)

    if interactive:
        fig, axes = plt.subplots(nrows=2, ncols=1)
        ax1, ax2 = axes
        xdata = np.append(xp0, xp1[-1])
        ydata = np.append(yp0, yp1[-1])
        plt.sca(ax1)
        plt.plot(xdata, ydata, "b")
        plt.sca(ax2)
        im = plt.imshow(depths, cmap="viridis_r")  # noqa
        ch = plt.colorbar()  # noqa
        fname = os.path.join(os.path.expanduser("~"), f"quad_{name}_test.png")
        print(f"Saving image for {name} quad test... {fname}")
        plt.savefig(fname)
        plt.close()


if __name__ == "__main__":
    test_rupture_depth(interactive=True)
