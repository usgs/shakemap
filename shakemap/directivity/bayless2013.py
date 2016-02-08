#!/usr/bin/env python

import numpy as np
import openquake.hazardlib.geo as geo

import shakemap.shakelib.fault as fault
import shakemap.shakelib.ecef as ecef
from shakemap.shakelib.vector import Vector
from shakemap.shakelib.distance import getDistance


"""
Implements the Bayless and Somerville (2013) directivity model. This is an update to 
Somerville et al. (1997), which defines some of the input parameters. 

Fd is the directivity effect parameter, and it is evaluated as

Fd = (c0 + c1 * Fgeom) * Tcd * Tmw * Taz

Model is seprated into strike-slip and dip-slip categories:
SS: (abs(rake) >  0 & abs(rake) < 30) | (abs(rake) > 150 & abs(rake) < 180)
DS:  abs(rake) > 60 & abs(rake) < 120

Notes from Somerville et al. (1997): 
 * d is the same as Wrup in Rowshandel (2013)
 * Y = d/W
 * s is the same as Ls in Rowshandel (2013), except it should be computed for
     each quadrilateral separately. 
 * X = s/L

Multi-segment ruptures need to be treated somewhat differently than in Rowshandel (2013). 
Bayless and Somerville explicitly state that each quadrilateral should have a pseudo-
hypocenter: 
    'Define the pseudo-hypocenter for rupture of successive segments as the point on the 
     side edge of the fault segment that is closest to the side edge of the previous segment, 
     and that lies half way between the top and bottom of the fault. We assume that the fault 
     is segmented along strike, not updip. All geometric parameters are computed relative to 
     the pseudo-hypocenter.'

"""

class bayless2013(object):
    """
    Class for Bayless and Somerville (2013) directvity model. 
    """
    #-----------------------------------------------------------------------------
    # C0 adn C1 are for RotD50. One set for each mechanism (SS vs DS).
    # FN and FP are also available
    #-----------------------------------------------------------------------------    
    __T = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10])
    __C0ss = np.array([0.0, 0.0, -0.12, -0.175, -0.21, -0.235, -0.255, -0.275, -0.29, -0.3])
    __C1ss = np.array([0.0, 0.0, 0.075, 0.090, 0.095, 0.099, 0.103, 0.108, 0.112, 0.115])
    __C0ds = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -0.033, -0.089, -0.133, -0.16, -0.176])
    __C1ds = np.array([0.0, 0.0, 0.0, 0.0, 0.034, 0.093, 0.128, 0.15, 0.165, 0.179])
    
    def __init__(self, flt, hyp, sites, rake, M, T):
        """
        Constructor for bayless2013.
        :param flt:
            pyhsake fault object.
        :param hyp:
            hypocenter; list of three elements: lon, lat, depth (km). 
        :param sites:
            two-element list; first element is a 2d numpy array of longitudes, 
            second element is a 2d array of latitudes. 
             ** TODO: change so that this argument is a SitesContext
        :param rake:
            Float (scalar) for direction of slip on fault; motion of hangingwall
            relative to footwall; clockwise from strike vector. 
             ** Potentially allow for variable rake with subfault. 
        :param M:
            Float (scalar) for moment magnitude. 
        :param T:
            Period; Currently, only acceptable values are 0.5 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10
             ** TODO: interpolate for other periods. 
        """
        self.flt = flt
        self.hyp = hyp
        self.sites = sites
        self.rake = rake
        self.M = M
        self.T = T
        self.mesh = geo.mesh.Mesh(self.sites[0], self.sites[1], np.zeros_like(self.sites[1]))
        self.Rrup = np.reshape(getDistance('rrup', self.mesh, self.flt.Quadrilaterals), self.sites[0].shape)
        self.Rx = np.reshape(getDistance('rx', self.mesh, self.flt.Quadrilaterals), self.sites[0].shape)
        self.Ry = np.reshape(getDistance('ry0', self.mesh, self.flt.Quadrilaterals), self.sites[0].shape)
        self.W = self.flt.getIndividualWidths()
        self.L = self.flt.getIndividualTopLengths()
        self.getSlipCategory()
        # d is the length of dipping fault rupturing toward site; max[(Y*W),exp(0)]
        self.computeD()
        # s is the length of striking fault rupturing toward site; max[(X*L),exp(1)]
        # theta (see Figure 5 in SSGA97)
        self.computeThetaAndS()
        self.computeAz()
        
        # Magnitude taper (does not depend on mechanism)
        if self.M <= 5: 
            T_Mw = 0.0
        elif self.M > 5 and self.M < 6.5:
            T_Mw = 1 - (6.5 - M)/1.5
        else:
            T_Mw = 1.0
        
        if self.SlipCategory == 'SS':
            # need to add loop over quads, for now assuming one quad
            
            # Geometric directivity predictor:
            f_geom = np.log(self.s[0]) * (0.5 * np.cos(2*self.theta[0]) + 0.5)
            
            # Distance taper
            T_CD = np.ones_like(self.sites[0])
            ix = [(self.Rrup/self.L > 0.5) & (self.Rrup/self.L < 1)]
            T_CD[ix] = 1 - (self.Rrup[ix]/self.L - 0.5)/0.5
            T_CD[self.Rrup/self.L >= 1.0 ] = 0

            # Azimuth taper
            T_Az = 1.0
            
            # Select Coefficients
            ix = [self.T == self.__T]
            C0 = self.__C0ss[ix]
            C1 = self.__C1ss[ix]
            
        elif self.SlipCategory == 'DS':
            # need to add loop over quads, for now assuming one quad
            
            # Geometric directivity predictor:            
            f_geom = np.log(self.d[0]) * np.cos(self.Rx/self.W[0])
            
            # Distance taper
            T_CD = np.ones_like(self.sites[0])
            ix = [(self.Rrup/self.W > 1.5) & (self.Rrup/self.W < 2)]
            T_CD[ix] = 1 - (self.Rrup[ix]/self.W - 1.5)/0.5
            T_CD[self.Rrup/self.W >= 2.0 ] = 0
            
            # Azimuth taper
            T_Az = np.sin(np.abs(self.Az))**2
            
            # Select Coefficients
            ix = [self.T == self.__T]
            C0 = self.__C0ds[ix]
            C1 = self.__C1ds[ix]
            
#        else:

        self.fd = (C0 + C1*f_geom) * T_CD * T_Mw * T_Az
        

    
    def computeAz(self):
        Az = np.ones_like(self.Rx) * np.pi/2.0
        Az = Az * np.sign(self.Rx)
        ix = [self.Ry > 0]
        Az[ix] = np.arctan(self.Rx[ix]/self.Ry[ix])
        self.Az = Az
        

    def computeD(self):
        """
        Y = d/W, where d is the portion (in km) of the width of the fault which 
        ruptures up-dip from the hypocenter to the top of the fault.

        Unlike Rowshandel's Wrup, this is computed for each quad. 
        """
        hyp_ecef = Vector.fromPoint(geo.point.Point(self.hyp[0], self.hyp[1], self.hyp[2]))
        nquad = len(self.flt.Quadrilaterals)
        d = np.zeros(nquad)
        for i in range(0, nquad):
            P0,P1,P2,P3 = self.flt.Quadrilaterals[i]
            p0 = Vector.fromPoint(P0) # convert to ECEF
            p3 = Vector.fromPoint(P3)
            hp0 = p0 - hyp_ecef
            p3p0n = (p0 - p3).norm()
            d[i] = Vector.dot(p3p0n, hp0)/1000
        self.d = d.clip(min = 0)
    
    def computeThetaAndS(self):
        # Use an orthographic projection
        latmin = self.sites[1].min()
        latmax = self.sites[1].max()
        lonmin = self.sites[0].min()
        lonmax = self.sites[0].max()
        proj = geo.utils.get_orthographic_projection(lonmin, lonmax, latmax, latmin)
        
        # Get epi projection
        epi_x,epi_y = proj(self.hyp[0], self.hyp[1])
        
        # Get the lines for the top edge of the fault
        qds = self.flt.Quadrilaterals
        nq = len(qds)
        top_lat = np.array([[]])
        top_lon = np.array([[]])
        top_dep = np.array([[]])
        for j in range(0, nq):
            top_lat = np.append(top_lat, qds[j][0].latitude)
            top_lat = np.append(top_lat, qds[j][1].latitude)
            top_lon = np.append(top_lon, qds[j][0].longitude)
            top_lon = np.append(top_lon, qds[j][1].longitude)
            top_dep = np.append(top_dep, qds[j][0].depth)
            top_dep = np.append(top_dep, qds[j][1].depth)
        top_x,top_y = proj(top_lon, top_lat)
        
        slat = self.sites[1]
        slon = self.sites[0]
        nquad = len(self.flt.Quadrilaterals)
        self.s = [None]*nquad
        self.theta = [None]*nquad
        ni, nj = slon.shape
        for h in range(nquad):
            s = np.zeros_like(slat)
            theta = np.zeros_like(slat)
            for i in range(ni): 
                for j in range(nj):
                    # Convert to local orthographic
                    site_x,site_y = proj(slon[i, j], slat[i, j])
                    
                    # Shift so center is at epicenter
                    site_x2 = site_x - epi_x
                    site_y2 = site_y - epi_y
                    top_x2 = top_x - epi_x
                    top_y2 = top_y - epi_y
                    
                    # Angle to rotate to put site on x-axis
                    alpha = np.arctan2(site_y2, site_x2) # "theta" in Bayless and Somerville notation
                    axis = [0, 0, 1]
                    rmat = _rotation_matrix(axis, -alpha)
                    llr = [None] * len(top_lat)
                    
                    # Apply rotation to each point on trace
                    for k in range(len(top_lat)):
                        llr[k] = np.dot(rmat, [top_x2[k], top_y2[k], 0])
                    site3 = np.dot(rmat, [site_x2, site_y2, 0])
                    s[i, j] = np.min([np.max([a[0] for a in llr]), site3[0]])
                    theta[i, j] = alpha
            self.s[h] = s.clip(min = np.exp(1))
            # theta is defined to be between 0 and 90 deg
            sintheta = np.abs(np.sin(theta))
            costheta = np.abs(np.cos(theta))
            self.theta[h] = np.arctan2(costheta, sintheta)
    
    def getSlipCategory(self):
        """
        Sets self.SlipCategory based on rake angle. Can be SS for 
        strike slip, DS for dip slip, or Unspecified. 
        """
        arake = np.abs(self.rake)
        self.SlipCategory = 'Unspecified'
        if (arake >=  0 & arake <= 30) | (arake >= 150 & arake <= 180):
            self.SlipCategory = 'SS'
        if (arake >= 60) & (arake <= 120):
            self.SlipCategory = 'DS'

def _rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    Source: Response by 'unutbu' in this thread: 
    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2)
    b, c, d = -axis*np.sin(theta/2)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
