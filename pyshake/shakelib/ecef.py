#!/usr/bin/env python

import numpy as np

DEGREES_TO_RADIANS = np.PI/180.0

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
# // Coverts ECEF to ENU coordinates centered at given lat, lon
# void ecefToEnu(double lat, double lon, double x, double y, double z, double xr, double yr, double zr, double *e, double *n, double *u)
# {
#     double clat = cos(lat * DEGREES_TO_RADIANS);
#     double slat = sin(lat * DEGREES_TO_RADIANS);
#     double clon = cos(lon * DEGREES_TO_RADIANS);
#     double slon = sin(lon * DEGREES_TO_RADIANS);
#     double dx = x - xr;
#     double dy = y - yr;
#     double dz = z - zr;

#     *e = -slon*dx  + clon*dy;
#     *n = -slat*clon*dx - slat*slon*dy + clat*dz;
#     *u = clat*clon*dx + clat*slon*dy + slat*dz;
# }

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
