#!/usr/bin/env python

#stdlib imports
import os.path
import sys
from collections import OrderedDict

#third party
import numpy as np

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.distance import *
from shakemap.grind.vector import Vector
from shakemap.grind.fault import Fault

from openquake.hazardlib.geo.geodetic import npoints_towards
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection

def test_repi():
    print('Testing high level get_distance() method with Repi...')
    #this is the coordinate of a station used in Northridge
    slat = 34.0700
    slon = -118.1500
    sdep = 0.0

    hypolat = 34.2057
    hypolon = -118.5539
    hypodep = 17.5

    #NGA database Repi value
    output = 40.15
    
    mesh = Mesh(np.array([slon]),np.array([slat]),np.array([0.0]))

    #quadlist doesn't get used for Repi, fill with zeros
    quadlist = []
    P0 = Point(0,0,0.0)
    P1 = Point(0,0,0.0)
    P2 = Point(0,0,0) 
    P3 = Point(0,0,0)
    quad = (P0,P1,P2,P3)
    quadlist.append(quad)
    repi = get_distance('repi',mesh,quadlist=quadlist,mypoint=Point(hypolon,hypolat,hypodep))[0]
    np.testing.assert_almost_equal(repi,output,decimal=1)
    print('Passed high level get_distance() method with Repi...')

def test_rhypo():
    print('Testing high level get_distance() method with Rhypo...')
    #this is the coordinate of a station used in Northridge
    slat = 34.0700
    slon = -118.1500
    sdep = 0.0

    hypolat = 34.2057
    hypolon = -118.5539
    hypodep = 17.5

    #NGA database Rhypo value
    output = 43.80
    
    mesh = Mesh(np.array([slon]),np.array([slat]),np.array([0.0]))

    #quadlist doesn't get used for Rhypo, fill with zeros
    quadlist = []
    P0 = Point(0,0,0.0)
    P1 = Point(0,0,0.0)
    P2 = Point(0,0,0) 
    P3 = Point(0,0,0)
    quad = (P0,P1,P2,P3)
    quadlist.append(quad)
    rhypo = get_distance('rhypo',mesh,quadlist=quadlist,mypoint=Point(hypolon,hypolat,hypodep))[0]
    np.testing.assert_almost_equal(rhypo,output,decimal=1)
    print('Passed high level get_distance() method with Rhypo...')

def test_rrup():
    print('Testing low level calc_rrup() method...')
    quads = []
    quad = (Point(-118.598000,34.387000,5.0000),
            Point(-118.431819,34.301105,5.0000),
            Point(-118.537983,34.160984,20.4270),
            Point(-118.703995,34.246737,20.4269))
    quads.append(quad)
    P0 = Vector.fromPoint(quad[0])
    P1 = Vector.fromPoint(quad[1])
    P2 = Vector.fromPoint(quad[2])
    P3 = Vector.fromPoint(quad[3])
    #this is the coordinate of a station used in Northridge
    slat = 34.0700
    slon = -118.1500
    sdep = 0.0
    p0 = Vector.fromPoint(Point(slon,slat,sdep))
    points = np.array([p0.x,p0.y,p0.z]).reshape(1,3)

    rrup = calc_rupture_distance(P0,P1,P2,P3,points)[0]
    output = 36.77

    np.testing.assert_almost_equal(rrup,output,decimal=0)
    print('Passed low level calc_rrup() method...')

    print('Testing high level get_distance() method with Rrup...')
    mesh = Mesh(np.array([slon]),np.array([slat]),np.array([sdep]))
    
    rrupd = get_distance('rrup',mesh,quadlist=quads)[0]
    np.testing.assert_almost_equal(rrupd,output,decimal=0)
    print('Passed high level get_distance() method with Rrup...')

def test_rjb():
    print('Testing low level calc_rjb() method...')
    quads = []
    quad = (Point(-118.598000,34.387000,0.0000),
            Point(-118.431819,34.301105,0.0000),
            Point(-118.537983,34.160984,0.0),
            Point(-118.703995,34.246737,0.0))
    quads.append(quad)
    P0 = Vector.fromPoint(quad[0])
    P1 = Vector.fromPoint(quad[1])
    P2 = Vector.fromPoint(quad[2])
    P3 = Vector.fromPoint(quad[3])
    #this is the coordinate of a station used in Northridge
    slat = 34.0700
    slon = -118.1500
    sdep = 0.0
    p0 = Vector.fromPoint(Point(slon,slat,sdep))
    points = np.array([p0.x,p0.y,p0.z]).reshape(1,3)

    rjb = calc_rupture_distance(P0,P1,P2,P3,points)[0]
    output = 35.66

    np.testing.assert_almost_equal(rjb,output,decimal=0)
    print('Passed low level calc_rjb() method...')

    print('Testing high level get_distance() method with Rjb...')
    mesh = Mesh(np.array([slon]),np.array([slat]),np.array([sdep]))
    
    rjbd = get_distance('rjb',mesh,quadlist=quads)[0]
    np.testing.assert_almost_equal(rjbd,output,decimal=0)
    print('Passed high level get_distance() method with Rjb...')

    
def test_rx():
    print('Testing low level calc_rx() method...')
    #this is the coordinate of a station used in Northridge
    slat = 34.0700
    slon = -118.1500
    sdep = 0.0

    #this is the starting point of the top edge of rupture for the Northridge fault
    p0lon = -118.598000
    p0lat = 34.387000

    #this is the ending point of the top edge of rupture for the Northridge fault
    p1lon = -118.431819
    p1lat = 34.301105

    west = min(p0lon,p1lon)
    east = max(p0lon,p1lon)
    south = min(p0lat,p1lat)
    north = max(p0lat,p1lat)
    proj = get_orthographic_projection(west,east,north,south) #projected coordinates are in km
    p0x,p0y = proj(p0lon,p0lat)
    p1x,p1y = proj(p1lon,p1lat)
    P0 = Vector(p0x*1000,p0y*1000,0.0)
    P1 = Vector(p1x*1000,p1y*1000,0.0)

    sx,sy = proj(slon,slat)
    station = np.array([sx*1000,sy*1000,0.0]).reshape(1,3)
    rx = calc_rx_distance(P0,P1,station)[0][0]

    output = 7.98
    np.testing.assert_almost_equal(rx,output,decimal=1)
    print('Passed low level calc_rx() method...')

    print('Testing high level get_distance() method with Rx...')
    mesh = Mesh(np.array([slon]),np.array([slat]),np.array([0.0]))
    quadlist = []
    P0 = Point(p0lon,p0lat,0.0)
    P1 = Point(p1lon,p1lat,0.0)
    P2 = Point(0,0,0) #it doesn't matter for Rx, as it only uses the top edge
    P3 = Point(0,0,0) #it doesn't matter for Rx, as it only uses the top edge
    quad = (P0,P1,P2,P3)
    quadlist.append(quad)
    rxd = get_distance('rx',mesh,quadlist=quadlist)[0]
    np.testing.assert_almost_equal(rxd,output,decimal=1)
    print('Passed high level get_distance() method with Rx...')
    

if __name__ == '__main__':
    test_repi()
    test_rhypo()
    test_rrup()
    test_rjb()
    test_rx()
    
