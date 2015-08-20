#!/usr/bin/env python

import numpy as np

DEGREES_TO_RADIANS = np.pi/180.0

a = 6378137.0 # radius
e = 8.1819190842622e-2  #eccentricity

asq = np.power(a,2)
esq = np.power(e,2)

DEGREES_TO_RADIANS = np.PI/180.0

def ecef2latlon(x,y,z):
    b = np.sqrt( asq * (1-esq) )
    bsq = np.power(b,2)
    ep = np.sqrt( (asq - bsq)/bsq)
    p = np.sqrt( np.power(x,2) + np.power(y,2) )
    th = np.atan2(a*z, b*p)
    
    lon = np.atan2(y,x)
    lat = np.atan2( (z + np.power(ep,2)*b*np.power(np.sin(th),3) ), (p - esq*a*np.power(np.cos(th),3)) )
    N = a/( np.sqrt(1-esq*np.power(np.sin(lat),2)) )
    alt = p / np.cos(lat) - N
    
    #mod lat to 0-2pi
    lon = lon % (2*np.PI)
    
    #correction for altitude near poles left out.
    
    return (lat, lon, alt)

def latlon2ecef(lat,lon,alt):
    clat = np.cos(lat * DEGREES_TO_RADIANS)
    slat = np.sin(lat * DEGREES_TO_RADIANS)
    clon = np.cos(lon * DEGREES_TO_RADIANS)
    slon = np.sin(lon * DEGREES_TO_RADIANS)
    
    N = WGS84_A / np.sqrt(1.0 - WGS84_E * WGS84_E * slat * slat)

    x = (N + alt) * clat * clon
    y = (N + alt) * clat * slon
    z = (N * (1.0 - WGS84_E * WGS84_E) + alt) * slat

    return (x,y,z)
