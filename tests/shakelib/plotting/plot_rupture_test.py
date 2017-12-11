#!/usr/bin/env python

import os.path
import sys

from openquake.hazardlib.geo.geodetic import azimuth
import numpy as np
import matplotlib.pyplot as plt

from shakelib.rupture.quad_rupture import QuadRupture
from shakelib.rupture.origin import Origin
from shakelib.plotting.plotrupture import plot_rupture_wire3d, map_rupture

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

MAX_DEPTH = 70
DIP = 17


def test_plot_rupture(interactive=False):
    xp0 = np.array([-90.898000])
    xp1 = np.array([-91.308000])
    yp0 = np.array([12.584000])
    yp1 = np.array([12.832000])
    zp = [0.0]
    strike = azimuth(yp0[0], xp0[0], yp1[0], xp1[0])
    origin = Origin({'lat': 0.0,
                     'lon': 0.0,
                     'depth': 0.0,
                     'mag': 5.5,
                     'eventsourcecode': 'abcd'})
    interface_width = MAX_DEPTH / np.sin(np.radians(DIP))
    widths = np.ones(xp0.shape) * interface_width
    dips = np.ones(xp0.shape) * DIP
    strike = [strike]
    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths, dips, origin, strike=strike)
    plot_rupture_wire3d(rupture)
    if interactive:
        fname = os.path.join(os.path.expanduser('~'), 'rupture_wire_plot.png')
        plt.savefig(fname)
        print('Wire 3D plot saved to %s.  Delete this file if you wish.'
              % fname)


def test_map_rupture(interactive=False):
    xp0 = np.array([-90.898000])
    xp1 = np.array([-91.308000])
    yp0 = np.array([12.584000])
    yp1 = np.array([12.832000])
    zp = [0.0]
    strike = azimuth(yp0[0], xp0[0], yp1[0], xp1[0])
    origin = Origin({'lat': 0.0,
                     'lon': 0.0,
                     'depth': 0.0,
                     'mag': 5.5,
                     'eventsourcecode': 'abcd'})
    interface_width = MAX_DEPTH / np.sin(np.radians(DIP))
    widths = np.ones(xp0.shape) * interface_width
    dips = np.ones(xp0.shape) * DIP
    strike = [strike]
    rupture = QuadRupture.fromTrace(
        xp0, yp0, xp1, yp1, zp, widths, dips, origin, strike=strike)
    map_rupture(rupture)
    if interactive:
        fname = os.path.join(os.path.expanduser('~'), 'rupture_map.png')
        plt.savefig(fname)
        print('Rupture map plot saved to %s.  Delete this file if you wish.'
              % fname)


if __name__ == '__main__':
    test_plot_rupture(interactive=True)
    test_map_rupture(interactive=True)
