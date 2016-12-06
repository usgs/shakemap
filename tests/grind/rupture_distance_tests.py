import os
import sys
import numpy as np

import matplotlib.pyplot as plt

import openquake.hazardlib.geo as geo

from shakemap.grind.rupture import EdgeRupture
from shakemap.grind.rupture import QuadRupture
from shakemap.grind.origin import Origin
from shakemap.grind.sites import Sites
from shakemap.grind.distance import get_distance

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_multisegment_discordant():
    # The one thing that isn't check above is discordancy for segments
    # with multiple quads. For this, we need a synthetic example. 
    x0 = np.array([0,   1, -1, 10,   9,  7])
    y0 = np.array([0,  10, 20, 40,  35, 30])
    z0 = np.array([0,   0,  0,  0,   0,  0])
    x1 = np.array([1,  -1,  0,  9,   7,  6])
    y1 = np.array([10, 20, 30, 35,  30, 25])
    z1 = np.array([0,   0,  0,  0,   0,  0])
    x2 = np.array([3,   1,  2,  7,   5,  4])
    y2 = np.array([10, 20, 30, 35,  30, 25])
    z2 = np.array([10, 10, 10, 10,  10, 10])
    x3 = np.array([2,   3,  1,  8,   7,  5])
    y3 = np.array([0,  10, 20, 40,  35, 30])
    z3 = np.array([10, 10, 10, 10,  10, 10])

    epilat = 32.15270
    epilon = -115.30500
    proj = geo.utils.get_orthographic_projection(epilon-1, epilon+1, epilat+1, epilat-1)
    lon0,lat0 = proj(x0, y0, reverse = True)
    lon1,lat1 = proj(x1, y1, reverse = True)
    lon2,lat2 = proj(x2, y2, reverse = True)
    lon3,lat3 = proj(x3, y3, reverse = True)

    rup = QuadRupture.fromVertices(
        lon0, lat0, z0, lon1, lat1, z1, lon2, lat2, z2, lon3, lat3, z3,
        group_index = [0, 0, 0, 1, 1, 1])
    # Make an Origin object; most of the 'event' values don't matter for this example
    event = {'lat': 0,  'lon': 0, 'depth':0, 'mag': 7.2, 
             'id':'', 'locstring':'', 'type':'ALL'}
    origin = Origin(event)

    # Sites
    buf = 0.25
    lat = np.linspace(np.nanmin(rup.lats)-buf, np.nanmax(rup.lats)+buf, 350)
    lon = np.linspace(np.nanmin(rup.lons)-buf, np.nanmax(rup.lons)+buf, 350)
    lons, lats = np.meshgrid(lon, lat)
    dep = np.zeros_like(lons)
    x,y = proj(lon, lat)
    rupx,rupy = proj(rup.lons, rup.lats)
    # Calculate U and T
    dtypes = ['U', 'T']
    dists = get_distance(dtypes, lats, lons, dep, origin, rup)

    if False:
        # Plot T
        fig = plt.figure(figsize=(6,6))
        CS = plt.contourf(x, y, dists['T'], 
                          cmap = plt.cm.coolwarm,
                          levels = np.arange(-120, 120, 5),
                          extend = "both")
        CS2 = plt.contour(CS, hold = 'on', colors = 'k')
        plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=14)
        plt.plot(rupx, rupy, color = 'g', linewidth = 3)
        plt.xlim([np.min(x), np.max(x)])
        plt.ylim([np.min(y), np.max(y)])
        cbar = plt.colorbar(CS)

        # Plot U
        fig = plt.figure(figsize=(6,6))
        CS = plt.contourf(x, y, dists['U'], 
                          cmap = plt.cm.coolwarm,
                          levels = np.arange(-30, 70, 5),
                          extend = "both")
        CS2 = plt.contour(CS, hold = 'on', colors = 'k')
        plt.clabel(CS2, fmt='%2.2f', colors='k', fontsize=14)
        plt.plot(rupx, rupy, color = 'g', linewidth = 3)
        plt.xlim([np.min(x), np.max(x)])
        plt.ylim([np.min(y), np.max(y)])
        cbar = plt.colorbar(CS)



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

    qrup = QuadRupture.fromTrace(xp0, yp0, xp1, yp1, zp, widths, dips)
    rrup_q = qrup.computeRrup(lon, lat, dep)
    rjb_q = qrup.computeRjb(lon, lat, dep)

    # Construct equivalent EdgeRupture
    toplons = np.array([-122.0, -121.7, -122.5, -122.3])
    toplats = np.array([37.1, 37.2, 37.4, 37.2])
    topdeps = np.array([0, 0, 6, 6])
    botlons = np.array([-121.886864, -121.587568, -122.635467, -122.435338])
    botlats = np.array([36.884527, 36.984246, 37.314035,  37.114261])
    botdeps = np.array([15.0000, 14.9998, 18.8558, 18.8559])
    group_index = [0, 0, 1, 1]

    erup = EdgeRupture(toplons, toplats, topdeps, botlons, botlats, botdeps, group_index)
    rrup_e = erup.computeRrup(lon, lat, dep, mesh_dx=0.5)
    rjb_e = erup.computeRjb(lon, lat, dep, mesh_dx=0.5)

    # Check that QuadRupture and EdgeRupture give the same result
    # (we check the absolute values of QuadRupture elsewhere)
    np.testing.assert_allclose(rrup_e, rrup_q, atol=0.35)
    np.testing.assert_allclose(rjb_e, rjb_q, atol=0.35)


    # For ploting
    #plt.imshow(rjb_q, interpolation="none")
    #plt.imshow(rjb_e, interpolation="none")

    #fig = plt.contourf(lon, lat, rrup_e,
    #                   levels = range(0, 100, 1),
    #                   cmap=plt.cm.spectral)
    #cbar = plt.colorbar(fig)
    #for q in qrup.getQuadrilaterals():
    #    x = [a.longitude for a in q]+[q[0].longitude]
    #    y = [a.latitude for a in q]+[q[0].latitude]
    #    plt.plot(x, y, 'r')


if __name__ == "__main__":
    test_multisegment_discordant()
    test_EdgeRupture_vs_QuadRupture()
