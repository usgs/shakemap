#!/usr/bin/env python

import numpy as np

DEGREES_TO_RADIANS = np.pi/180.0

a = 6378137.0 # radius
e = 8.1819190842622e-2  #eccentricity

asq = np.power(a,2)
esq = np.power(e,2)

# void latLonToEcef(double lat, double lon, double alt, double *x, double *y, double *z)
# {   
#     double clat = cos(lat * DEGREES_TO_RADIANS);
#     Double slat = sin(lat * DEGREES_TO_RADIANS);
#     double clon = cos(lon * DEGREES_TO_RADIANS);
#     double slon = sin(lon * DEGREES_TO_RADIANS);

#     double N = WGS84_A / sqrt(1.0 - WGS84_E * WGS84_E * slat * slat);

#     *X = (N + alt) * clat * clon;
#     *y = (N + alt) * clat * slon;
#     *z = (N * (1.0 - WGS84_E * WGS84_E) + alt) * slat;
# }
# /*
# *
# *  ECEF - Earth Centered Earth Fixed
# *   
# *  LLA - Lat Lon Alt
# *
# *  ported from matlab code at
# *  https://gist.github.com/1536054
# *     and
# *  https://gist.github.com/1536056
# */

# // WGS84 ellipsoid constants
# private final double a = 6378137; // radius
# private final double e = 8.1819190842622e-2;  // eccentricity

# private final double asq = Math.pow(a,2);
# private final double esq = Math.pow(e,2);

# private double[] ecef2lla(double[] ecef){
#   double x = ecef[0];
#   double y = ecef[1];
#   double z = ecef[2];

#   double b = Math.sqrt( asq * (1-esq) );
#   double bsq = Math.pow(b,2);
#   double ep = Math.sqrt( (asq - bsq)/bsq);
#   double p = Math.sqrt( Math.pow(x,2) + Math.pow(y,2) );
#   double th = Math.atan2(a*z, b*p);

#   double lon = Math.atan2(y,x);
#   double lat = Math.atan2( (z + Math.pow(ep,2)*b*Math.pow(Math.sin(th),3) ), (p - esq*a*Math.pow(Math.cos(th),3)) );
#   double N = a/( Math.sqrt(1-esq*Math.pow(Math.sin(lat),2)) );
#   double alt = p / Math.cos(lat) - N;

#   // mod lat to 0-2pi
#   lon = lon % (2*Math.PI);

#   // correction for altitude near poles left out.

#   double[] ret = {lat, lon, alt};

#   return ret;
# }

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
