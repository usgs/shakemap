#!/usr/bin/env python

# stdlib imports
import sys
import os.path

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..'))
# put this at the front of the system path, ignoring any installed mapio stuff
sys.path.insert(0, shakedir)

import numpy as np

# local imports
from shakemap.grind.sites import Sites


def test(vs30file=None):
    vs30file = os.path.join(shakedir, 'data/Vs30_test.grd')
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 0.0083 * 5
    yspan = 0.0083 * 5
    mysite = Sites.createFromCenter(cx, cy, xspan, yspan, dx, dy,
                                    vs30File=vs30file, padding=True,
                                    resample=False)
    sc = mysite.getSitesContext()
    scsamp = mysite.sampleFromSites(np.array([34.1, 34.111]),
                                    np.array([-118.2, -118.222]))

    xmin = -118.234
    xmax = -118.110
    ymin = 34.01
    ymax = 34.21
    dx = 0.0083
    dy = 0.0083
    mysite = Sites.createFromBounds(xmin, xmax, ymin, ymax, dx, dy,
                                    vs30File=vs30file, padding=False,
                                    resample=False)

if __name__ == '__main__':
    vs30file = None
    if len(sys.argv) > 1:
        vs30file = sys.argv[1]
    test(vs30file=vs30file)
