#!/usr/bin/env python

import numpy as np
import openquake.hazardlib.geo as geo

import shakemap.shakelib.fault as fault

from shakemap.shakelib.vector import Vector


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
        self.setPredictors()
    
        
    def setPredictors(self):
        """
        Set the input/predictor variables Fgeom, Tcd, Tmw, Taz. 
        """
        self.W = self.flt.getIndividualWidths()
        self.L = self.flt.getIndividualTopLengths()
        self.computeY()
        
#        theta
#        d
#        Rx
#        W
#        Rrup
#        Az

#        if(self.SlipCategory == 'SS'):
#            
#        elif(self.SlipCategory == 'DS'):
#            
#        else:
            
    def computeX(self):
        """
        X = s/L, where s is Rowshandel's Ls except for each quad individually. 
        Since this is a function of site, it returns a list of arrays, each array 
        is for each quad. 
        """
        # Get lat/lon arrays
        slat = self.sites[1]
        slon = self.sites[0]
        # Convert site to ECEF
        site_ecef_x = np.ones_like(slat)
        site_ecef_y = np.ones_like(slat)
        site_ecef_z = np.ones_like(slat)
        
        # Make a 2x(#number of sites) matrix of site locations (rows are x, y, z) in ECEF
        site_ecef_x, site_ecef_y, site_ecef_z = ecef.latlon2ecef(slat, slon, np.zeros(slon.shape) )
        site_mat = np.array([np.reshape(site_ecef_x, (-1,)),
                             np.reshape(site_ecef_y, (-1,))])
        
        hyp_ecef = Vector.fromPoint(geo.point.Point(self.hyp[0], self.hyp[1]), self.hyp[2])
        epi_ecef_col = np.array([[hyp_ecef[0]],[hyp_ecef[1]]])
        
        nquad = len(self.flt.Quadrilaterals)
        X = [None] * nquad
        for i in range(0, nquad):
            P0,P1,P2,P3 = self.flt.Quadrilaterals[i]
            p0 = Vector.fromPoint(P0) # convert to ECEF
            p1 = Vector.fromPoint(P1) 
            p3 = Vector.fromPoint(P3)
            str_vec = (p1 - p0).norm()
            str_vec_col = np.array([[str_vec.x], [str_vec.y], [str_vec.z]])
            str_dot_site = np.sum(str_vec_col * (site_mat-epi_ecef_col), axis = 0)
            
            X[i] = s/self.L[i]
        self.Y = Y



    def computeY(self):
        """
        Y = d/W, where d is the portion (in km) of the width of the fault which 
        ruptures up-dip from the hypocenter to the top of the fault.

        Unlike Rowshandel's Wrup, this is computed for each quad. 
        """
        hyp_ecef = Vector.fromPoint(geo.point.Point(self.hyp[0], self.hyp[1], self.hyp[2]))
        nquad = len(self.flt.Quadrilaterals)
        Y = np.zeros(nquad)
        for i in range(0, nquad):
            P0,P1,P2,P3 = self.flt.Quadrilaterals[i]
            p0 = Vector.fromPoint(P0) # convert to ECEF
            p3 = Vector.fromPoint(P3)
            hp0 = p0 - hyp_ecef
            p3p0n = (p0 - p3).norm()
            d = Vector.dot(p3p0n, hp0)/1000
            Y[i] = d/self.W[i]
        self.Y = Y
        
    def getSlipCategory(self):
        """
        Sets self.SlipCategory based on rake angle. Can be SS for 
        strike slip, DS for dip slip, or Unspecified. 
        """
        arake = np.abs(rake)
        self.SlipCategory = 'Unspecified'
        if (arake >  0 & arake < 30) | (arake > 150 & arake < 180):
            self.SlipCategory = 'SS'
        if arake > 60 & arake < 120:
            self.SlipCategory = 'DS'
        

