#!/usr/bin/env python

#stdlib imports
import os.path
import sys
from collections import OrderedDict

#third party
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as time

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.distance import *
from shakemap.grind.vector import Vector
from shakemap.grind.fault import Fault
from shakemap.utils.timeutils import ShakeDateTime
from shakemap.grind.source import Source

from openquake.hazardlib.geo.geodetic import npoints_towards
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection
from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014

def test_repi():
    print('Testing high level get_distance() method with Repi...')
    #this is the coordinate of a station used in Northridge
    slat = np.array([34.0700])
    slon = np.array([-118.1500])
    sdep = np.array([0.0])

    hypolat = 34.2057
    hypolon = -118.5539
    hypodep = 17.5

    # NGA database Repi value
    output = 40.15
    
    # quadlist doesn't get used for Repi, fill with zeros
    quadlist = []
    P0 = Point(0, 0, 0.0)
    P1 = Point(0, 0, 0.0)
    P2 = Point(0, 0, 0) 
    P3 = Point(0, 0, 0)
    quad = (P0, P1, P2, P3)
    quadlist.append(quad)
    repi = get_distance('repi', slat, slon, sdep, quadlist = quadlist,
                        hypo = Point(hypolon, hypolat, hypodep))['repi']
    np.testing.assert_almost_equal(repi, output, decimal = 1)
    print('Passed high level get_distance() method with Repi...')

def test_rhypo():
    print('Testing high level get_distance() method with Rhypo...')
    # this is the coordinate of a station used in Northridge
    slat = np.array([34.0700])
    slon = np.array([-118.1500])
    sdep = np.array([0.0])

    hypolat = 34.2057
    hypolon = -118.5539
    hypodep = 17.5

    # NGA database Rhypo value
    output = 43.80
    
    # quadlist doesn't get used for Rhypo, fill with zeros
    quadlist = []
    P0 = Point(0, 0, 0.0)
    P1 = Point(0, 0, 0.0)
    P2 = Point(0, 0, 0) 
    P3 = Point(0, 0, 0)
    quad = (P0,P1,P2,P3)
    quadlist.append(quad)
    rhypo = get_distance('rhypo', slat, slon, sdep, quadlist = quadlist,
                         hypo = Point(hypolon, hypolat, hypodep))['rhypo']
    np.testing.assert_almost_equal(rhypo,output,decimal=1)
    print('Passed high level get_distance() method with Rhypo...')

def test_rrup():
    print('Testing low level calc_rrup() method...')
    quads = []
    quad = (Point(-118.598000, 34.387000, 5.0000),
            Point(-118.431819, 34.301105, 5.0000),
            Point(-118.537983, 34.160984, 20.4270),
            Point(-118.703995, 34.246737, 20.4269))
    quads.append(quad)
    P0, P1, P2, P3 = quad
    # this is the coordinate of a station used in Northridge
    slat = np.array([34.0700])
    slon = np.array([-118.1500])
    sdep = np.array([0.0])
    p0 = Vector.fromPoint(Point(slon, slat, sdep))
    points = np.array([p0.x, p0.y, p0.z]).reshape(1, 3)
    rrup = calc_rupture_distance(P0, P1, P2, P3, points)
    output = 36.77
    np.testing.assert_almost_equal(rrup, output, decimal = 0)
    print('Passed low level calc_rrup() method...')

    print('Testing high level get_distance() method with Rrup...')
    rrupd = get_distance('rrup', slat, slon, sdep, quadlist = quads)['rrup']
    np.testing.assert_almost_equal(rrupd,output,decimal = 0)
    print('Passed high level get_distance() method with Rrup...')

def test_rjb():
    print('Testing low level calc_rjb() method...')
    quads = []
    quad = (Point(-118.598000, 34.387000, 0.0000),
            Point(-118.431819, 34.301105, 0.0000),
            Point(-118.537983, 34.160984, 0.0),
            Point(-118.703995, 34.246737, 0.0))
    quads.append(quad)
    P0, P1, P2, P3 = quad
    # this is the coordinate of a station used in Northridge
    slat = np.array([34.0700])
    slon = np.array([-118.1500])
    sdep = np.array([0.0])
    p0 = Vector.fromPoint(Point(slon, slat, sdep))
    points = np.array([p0.x, p0.y, p0.z]).reshape(1, 3)
    rjb = calc_rupture_distance(P0, P1, P2, P3, points)
    output = 35.66
    np.testing.assert_almost_equal(rjb, output, decimal = 0)
    print('Passed low level calc_rjb() method...')

    print('Testing high level get_distance() method with Rjb...')
    rjbd = get_distance('rjb', slat, slon, sdep, quadlist = quads)['rjb']
    np.testing.assert_almost_equal(rjbd, output, decimal = 0)
    print('Passed high level get_distance() method with Rjb...')

    
def test_rx():
    print('Testing low level calc_rx() method...')
    # this is the coordinate of a station used in Northridge
    slat = np.array([34.0700])
    slon = np.array([-118.1500])
    sdep = np.array([0.0])

    # this is the starting point of the top edge of rupture for the Northridge fault
    p0lon = -118.598000
    p0lat = 34.387000

    # this is the ending point of the top edge of rupture for the Northridge fault
    p1lon = -118.431819
    p1lat = 34.301105
    
    P0 = Point(p0lon, p0lat, 0.0)
    P1 = Point(p1lon, p1lat, 0.0)

    rx = calc_t_i(P0, P1, slat, slon)

    output = 7.98
    np.testing.assert_almost_equal(rx, output, decimal = 1)
    print('Passed low level calc_rx() method...')

    print('Testing high level get_distance() method with Rx...')
#    mesh = Mesh(np.array([slon]), np.array([slat]), np.array([0.0]))
    quadlist = []
    P0 = Point(p0lon, p0lat, 0.0)
    P1 = Point(p1lon, p1lat, 0.0)
    P2 = Point(0, 0, 0) #it doesn't matter for Rx, as it only uses the top edge
    P3 = Point(0, 0, 0) #it doesn't matter for Rx, as it only uses the top edge
    quad = (P0, P1, P2, P3)
    quadlist.append(quad)
    rxd = get_distance('rx', slat, slon, sdep, quadlist = quadlist)['rx']
    np.testing.assert_almost_equal(rxd, output, decimal = 1)
    print('Passed high level get_distance() method with Rx...')

def test_chichi():
    print('Testing Chi-Chi...')
    # read in fault file
    f = '../data/0137A.POL'
    i0 = np.arange(0, 9*11*3, 11)
    i1 = i0 + 10
    cs = zip(i0, i1)
    df = pd.read_fwf(f, cs, skiprows = 2, nrows = 5, header = None)
    mat = df.as_matrix()
    ix = np.arange(0, 9*3, 3)
    iy = ix + 1
    iz = ix + 2
    x0 = mat[0, ix]
    x1 = mat[1, ix]
    x2 = mat[2, ix]
    x3 = mat[3, ix]
    y0 = mat[0, iy]
    y1 = mat[1, iy]
    y2 = mat[2, iy]
    y3 = mat[3, iy]
    # Depth, positive down
    z0 = np.abs(mat[0, iz])
    z1 = np.abs(mat[1, iz])
    z2 = np.abs(mat[2, iz])
    z3 = np.abs(mat[3, iz])
    epilat = 23.85
    epilon = 120.82
    proj = get_orthographic_projection(
        epilon-1, epilon+1, epilat+1, epilat-1)
    lon0,lat0 = proj(x0, y0, reverse = True)
    lon1,lat1 = proj(x1, y1, reverse = True)
    lon2,lat2 = proj(x2, y2, reverse = True)
    lon3,lat3 = proj(x3, y3, reverse = True)
    flt = Fault.fromVertices(
        lon0, lat0, z0, lon1, lat1, z1, lon2, lat2, z2, lon3, lat3, z3)
    ask14 = AbrahamsonEtAl2014()
    # event information doesn't matter...
    event = {'lat': 0,  'lon': 0, 'depth':0, 'mag': 7, 
             'id':'', 'locstring':'', 'type':'U', 
             'time':ShakeDateTime.utcfromtimestamp(int(time.time())), 
             'timezone':'UTC'}
    source = Source(event, flt)
    
    # Get NGA distances
    distfile = '../data/NGAW2_distances.csv'
    df = pd.read_csv(distfile)
    df2 = df.loc[df['EQID'] == 137]
    slat = df2['Station Latitude'].as_matrix()
    slon = df2['Station Longitude'].as_matrix()
    sdep = np.zeros(slat.shape)
    nga_repi = df2['EpiD (km)'].as_matrix()
    nga_rhypo = df2['HypD (km)'].as_matrix()
    nga_rrup = df2['ClstD (km)'].as_matrix()
    nga_rjb = df2['Joyner-Boore Dist. (km)'].as_matrix()
    nga_rx = df2['T'].as_matrix()
    
    dist = Distance(ask14, source, slat, slon, sdep)
    dctx = dist.getDistanceContext()
    fig = plt.figure(figsize=(8,8))
    plt.scatter(nga_rjb, dctx.rjb, alpha = 0.5, facecolors='none')
    plt.plot([0, nga_rjb.max()], [0, dctx.rjb.max()], 'b');
    plt.savefig('Chi-Chi_Rjb.png')
    fig = plt.figure(figsize=(8,8))
    plt.scatter(nga_rrup, dctx.rrup, alpha = 0.5, facecolors='none')
    plt.plot([0, nga_rrup.max()], [0, dctx.rrup.max()], 'b');
    plt.savefig('Chi-Chi_Rrup.png')
    fig = plt.figure(figsize=(8,8))
    plt.scatter(nga_rx, dctx.rx, alpha = 0.5, facecolors='none')
    plt.plot([nga_rx.min(), nga_rx.max()],
             [dctx.rx.min(), dctx.rx.max()], 'b');
    plt.savefig('Chi-Chi_Rx.png')

if __name__ == '__main__':
    test_repi()
    test_rhypo()
    test_rrup()
    test_rjb()
    test_rx()
    test_chichi()
