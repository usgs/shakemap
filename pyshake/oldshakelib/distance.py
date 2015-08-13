#!/usr/bin/env python

#third party imports
import numpy as np

def geoToECEF(phi,lam,h):
    phi = np.radians(phi)
    lam = np.radians(lam)
    f = 1/298.257223563
    a = 6378137.0
    esq = 2*f - np.power(f,2)
    N = a/(np.sqrt(1-esq*np.power(np.sin(phi),2)))
    X = (N + h) * (np.cos(phi) * np.cos(lam))
    Y = (N + h) * (np.cos(phi) * np.sin(lam))
    Z = (N * (1 - esq) + h) * np.sin(phi)
    return (X,Y,Z)

def ECEFToGeo(X,Y,Z):
    f = 1/298.257223563
    tol = 1e-3
    a = 6378137.0
    esq = 2*f - np.power(f,2)
    kappa0 = 1/(1 - esq)
    p = np.sqrt(np.power(X,2) + np.power(Y,2))
    kappai = kappa0
    kappa = kappai*1000
    while (kappa/kappai - 1) > tol:
        ci = (np.power(np.power(p,2) + (1 - esq) * np.power(Z,2) * np.power(kappai,3),3.0/2.0))/(a*esq)
        kappa = 1 + ((np.power(p,2) + (1 - esq) * np.power(Z,2) * np.power(kappai,3))/(ci - np.power(p,2)))
    h = (1/esq) * ((1/kappa) - (1/kappa0)) * np.sqrt(np.power(p,2) + (np.power(Z,2) * np.power(kappa,2)))
    phi = np.arctan((kappa * Z)/p)
    lam = np.arcsin(Y/(N + h) * np.cos(phi))
    phi = np.degrees(phi)
    lam = np.degrees(lam)
    return (phi,lam,h)

def getCartesianDistance(lat1,lon1,h1,lat2,lon2,h2):
    x1,y1,z1 = geoToECEF(lat1,lon1,h1)
    x2,y2,z2 = geoToECEF(lat2,lon2,h2)
    dist = np.sqrt(np.power(x1-x2,2) + np.power(y1-y2,2) + np.power(z1-z2,2))
    return dist
