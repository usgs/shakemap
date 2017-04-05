#!/usr/bin/env python

import numpy as np

DEGREES_TO_RADIANS = np.pi / 180.0
RADIANS_TO_DEGREES = 180.0 / np.pi

WGS84_A = 6378137.0  # radius
WGS84_E = 8.1819190842622e-2  # eccentricity

asq = np.power(WGS84_A, 2)
esq = np.power(WGS84_E, 2)

# adapted from Matlab code: https://gist.github.com/klucar/1536056


def ecef2latlon(x, y, z):
    """
    Convert Earth-Centered-Earth-Fixed (ECEF) cartesian coordinates
    to lat,lon,depth.

    :parameter x:
        A numpy array (or scalar value) of x coordinates (meters).
    :parameter y:
        A numpy array (or scalar value) of y coordinates (meters).
    :parameter z:
        A numpy array (or scalar value) of y coordinates (meters), positive UP.
    :return:
        Tuple of lat,lon,depth numpy arrays, where lat,lon are in dd and depth 
        is in km and positive DOWN.
    """
    inputIsScalar = False
    if not isinstance(x, np.ndarray):
        x = np.array([x])
        y = np.array([y])
        z = np.array([z])
        inputIsScalar = True
    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622e-2

    # calculations:
    b = np.sqrt(a**2 * (1 - e**2))
    ep = np.sqrt((a**2 - b**2) / b**2)
    p = np.sqrt(x**2 + y**2)
    th = np.arctan2(a * z, b * p)
    lon = np.arctan2(y, x)
    lat = np.arctan2(
        (z + ep**2 * b * np.sin(th)**3),
        (p - e**2 * a * np.cos(th)**3))
    N = a / np.sqrt(1 - e**2 * np.sin(lat)**2)
    alt = p / np.cos(lat) - N

    # return lon in range [0,2*pi)
    lon = np.mod(lon, 2 * np.pi)

    # correct for numerical instability in altitude near exact poles:
    # (after this correction, error is about 2 millimeters, which is about
    # the same as the numerical precision of the overall function)
    k = (np.abs(x) < 1) & (np.abs(y) < 1)
    alt[k] = np.abs(z[k]) - b

    # convert lat,lon to dd, and alt to depth positive DOWN in km
    lat = lat * RADIANS_TO_DEGREES
    lon = lon * RADIANS_TO_DEGREES
    lon[lon > 180] = lon[lon > 180] - 360.0
    dep = -alt / 1000.0
    # if input values were scalar, give that back to them
    if inputIsScalar:
        lat = lat[0]
        lon = lon[0]
        dep = dep[0]
    return (lat, lon, dep)


def latlon2ecef(lat, lon, dep):
    """
    Convert lat,lon,depth to Earth-Centered-Earth-Fixed (ECEF) 
    cartesian coordinates.

    :parameter lat:
        A numpy array (or scalar value) of latitude values (decimal degrees).
    :parameter lon:
        A numpy array (or scalar value) of longitude values (decimal degrees).
    :parameter dep:
        A numpy array (or scalar value) of depths (km), positive DOWN.
    :return:
        Tuple of x,y,z numpy arrays of ECEF coordinates.
    """
    alt = -dep * 1000.0  # convert to altitude (positive UP) in meters
    clat = np.cos(lat * DEGREES_TO_RADIANS)
    slat = np.sin(lat * DEGREES_TO_RADIANS)
    clon = np.cos(lon * DEGREES_TO_RADIANS)
    slon = np.sin(lon * DEGREES_TO_RADIANS)

    N = WGS84_A / np.sqrt(1.0 - WGS84_E * WGS84_E * slat * slat)

    x = (N + alt) * clat * clon
    y = (N + alt) * clat * slon
    z = (N * (1.0 - WGS84_E * WGS84_E) + alt) * slat

    return (x, y, z)
