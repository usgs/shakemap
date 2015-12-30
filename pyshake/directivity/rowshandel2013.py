#!/usr/bin/env python


import numpy as np
import openquake.hazardlib.geo as geo
import pyshake.shakelib.ecef as ecef
import pyshake.shakelib.fault as fault
from pyshake.shakelib.vector import Vector
from pyshake.shakelib.distance import calcRuptureDistance
from pyshake.shakelib.distance import getDistance
#from numba import jit

"""
Implements the Rowshandel (2013) directivity model. 

To do: 
    * Add checks on function arguments (e.g., mtype) for valid values. 
    * Add a validation function. 
    * Use np.clip for applying the minimum of zero to s-dot-q and p-dot-q
    * Optimize with numba?
"""

class rowshandel2013(object):
    """
    Class for Rowshandel (2013) directvity model. 
    """
    #-----------------------------------------------------------------------------
    # C1 adn C2 are GMPE-specific coefficients. Also, these are for Model-I option. 
    #-----------------------------------------------------------------------------
    __c1 = {'as': np.array([0.375,0.561,0.833,0.873,1.019,1.387,1.551]), 
            'ba': np.array([0.257,0.699,0.824,0.834,0.986,1.378,1.588]),
            'cb': np.array([0.149,0.377,0.570,0.596,0.772,1.152,np.nan]), 
            'cy': np.array([0.079,0.321,0.629,0.693,0.898,1.258,1.369]),
            'id': np.array([0.325,0.803,1.544,1.576,1.823,2.089,2.250]),
            'T': np.array([1, 2, 3, 4, 5, 7.5, 10])}
    __c1['mean'] = (__c1['as'] + __c1['ba'] + __c1['cb'] + __c1['cy'] + __c1['id'])/5
    __c2 = {'as': np.array([0.028,-0.071,-0.130,-0.146,-0.178,-0.250,-0.268]), 
            'ba': np.array([0.024,-0.097,-0.133,-0.145,-0.182,-0.273,-0.319]),
            'cb': np.array([-0.007,-0.057,-0.095,-0.106,-0.144,-0.230,np.nan]), 
            'cy': np.array([-0.003,-0.404,-0.098,-0.115,-0.156,-0.224,-0.235]),
            'id': np.array([-0.019,-0.082,-0.228,-0.255,-0.312,-0.358,-0.369]),
            'T': np.array([1, 2, 3, 4, 5, 7.5, 10])}
    __c2['mean'] = (__c2['as'] + __c2['ba'] + __c2['cb'] + __c2['cy'] + __c2['id'])/5

    def __init__(self, flt, hyp, sites, rake, dx, M, T, simpleDT = False):
        """
        Constructor for rowshandel2013.
        :param fault:
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
        :param dx:
            Float (scalar) for target mesh spacing for subfaults in km. The mesh
            snaps to the edges of the quadrilaterals so the actual mesh spacing 
            will not equal this exactly, and will also be variable. 
        :param M:
            Float (scalar) for moment magnitude. 
        :param T:
            Period; Currently, only acceptable values are 1, 2, 3, 4, 5, 7.5
             ** TODO: interpolate for other periods. 
        :param simpleDT:
            Boolean; should the simpler DT equation be used? Usually False. 
        """
        self.flt = flt
        self.hyp = hyp
        self.sites = sites
        self.rake = rake
        self.dx = dx
        self.M = M
        self.T = T
        self.simpleDT = simpleDT
        
        self.computeWrup()
        self.computeLD()
        
        self.computeDT()
        self.computeWP()
        self.computeXiPrime()
        self.computeFd()

    
    def computeFd(self):
        """
        Fd is the term that modifies the GMPE output as an additive term (in log space). 
        ln(Y) = f(M, R, ...) + Fd
        """
        c1sel = self.__c1['mean'][self.__c1['T'] == self.T]
        c2sel = self.__c2['mean'][self.__c2['T'] == self.T]
        # Eqn 3.13
        self.Fd = (c1sel*self.xi_prime + c2sel) * self.LD * self.DT * self.WP


    def computeWrup(self):
        """
        Wrup is the portion (in km) of the width of the fault which 
        ruptures up-dip from the hypocenter to the top of the fault.
          * This is ambiguous for faults with varible top of rupture (not 
            allowed in NGA). For now, lets just compute this for the 
            quad where the hypocenter is located.
          * Alternative is to compute max Wrup for the different quads. 
        """
        nquad = len(self.flt.Quadrilaterals)
        
        #-------------------------------------------
        # First find which quad the hypocenter is on
        #-------------------------------------------
        
        x,y,z = ecef.latlon2ecef(self.hyp[1],self.hyp[0],self.hyp[2])
        hyp_ecef = np.array([[x, y, z]])
        qdist = np.zeros_like([nquad])
        for i in range(0, nquad):
            P0,P1,P2,P3 = self.flt.Quadrilaterals[i]
            p0 = Vector.fromPoint(P0) # convert to ECEF
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            qdist[i] = calcRuptureDistance(p0, p1, p2, p3, hyp_ecef)
        ind = int(np.where(qdist == np.min(qdist))[0])  ## check that this doesn't break with more than one quad
        q = self.flt.Quadrilaterals[ind]
        
        #-------------------------------------------
        # Compute Wrup on that quad
        #-------------------------------------------
        
        pp0 = Vector.fromPoint(geo.point.Point(q[0].longitude, q[0].latitude, q[0].depth) )
        pp3 = Vector.fromPoint(geo.point.Point(q[3].longitude, q[3].latitude, q[3].depth) )
        hyp_ecef = Vector.fromPoint(geo.point.Point(self.hyp[0], self.hyp[1], self.hyp[2]))
        hp0 = pp0 - hyp_ecef
        p3p0n = (pp0 - pp3).norm()
        self.Wrup = Vector.dot(p3p0n, hp0)/1000

    
    def computeXiPrime(self, mtype = 1, a_weight = 0.5):
        """
        Computes the xi' value. 
        :param mtype:
            Integer, either 1 or 2; 1 for adding only positive components of dot products, 
            2 for adding all components of dot products. 
        :param a_weight:
            Float (scalar) weighting factor: xi' = a_weight*xi_s' + (1-a_weight)*xi_p'
            a_weight = 1 means only use s-dot-q; a_weight = 0 means only use p-dot-q. 
        """
        hypo_ecef = Vector.fromPoint(geo.point.Point(self.hyp[0], self.hyp[1], self.hyp[2]))
        epi_ll = Vector(self.hyp[0], self.hyp[1], 0)
        epi_ecef = Vector.fromPoint(geo.point.Point(epi_ll.x, epi_ll.y, 0))
        
        slat = self.sites[1]
        slon = self.sites[0]
        
        # Convert site to ECEF:
        site_ecef_x = np.ones_like(slat)
        site_ecef_y = np.ones_like(slat)
        site_ecef_z = np.ones_like(slat)
        
        # Make a 3x(#number of sites) matrix of site locations (rows are x, y, z) in ECEF
        site_ecef_x, site_ecef_y, site_ecef_z = ecef.latlon2ecef(slat, slon, np.zeros(slon.shape) )
        site_mat = np.array([np.reshape(site_ecef_x, (-1,)),
                             np.reshape(site_ecef_y, (-1,)),
                             np.reshape(site_ecef_z, (-1,))])
        
        xi_prime_unnormalized = np.zeros_like(slat)        
        
        # Normalize by total number of subfaults, can't know till loop is finished
        n_sub_faults = 0
        
        for k in range(len(self.flt.Quadrilaterals)):
            # Selecte a quad
            q = self.flt.Quadrilaterals[k]
            
            # Quad mesh (ECEF coords)
            mesh = fault.getQuadMesh(q, self.dx)
            
            # Rupture plane normal vector (ECEF coords)
            rpnv = fault.getQuadNormal(q)
            
            # Unit slip vector (ECEF coords)
            slpv = fault.getQuadSlip(q, self.rake)
            scol = np.array([[slpv.x], [slpv.y], [slpv.z]]) # column vector
            
            # Make 3x(i*j) matrix of cp
            ni, nj = mesh['llx'].shape
            
            cp_mat = np.array([np.reshape(mesh['cpx'], (-1,)),
                               np.reshape(mesh['cpy'],(-1,)),
                               np.reshape(mesh['cpz'], (-1,))])
            
            # Compute matrix of p vectors
            hypcol = np.array([[hypo_ecef.x], [hypo_ecef.y], [hypo_ecef.z]])
            pmat = cp_mat - hypcol
            mag = np.sqrt(np.sum(pmat*pmat, axis = 0))
            pmatnorm = pmat/mag
            
            xi_prime = np.zeros([site_mat.shape[1]])
            n_sub_faults = n_sub_faults + cp_mat.shape[1]
            
            # Loop over sites
            for i in range(site_mat.shape[1]):
                sitecol = np.array([[site_mat[0, i]], [site_mat[1, i]], [site_mat[2, i]]])
                qmat = sitecol - cp_mat  # 3x(ni*nj)
                mag = np.sqrt(np.sum(qmat*qmat, axis = 0))
                qmatnorm = qmat/mag
                # Dot products
                if mtype == 1:
                    pdotq = np.sum(pmatnorm * qmatnorm, axis = 0).clip(min = 0)
                    sdotq = np.abs(np.sum(scol * qmatnorm, axis = 0))
                elif mtype == 2:
                    pdotq = np.sum(pmatnorm * qmatnorm, axis = 0)
                    sdotq = np.abs(np.sum(scol * qmatnorm, axis = 0))
                xi_prime[i] =a_weight*np.sum(sdotq) + (1-a_weight)*np.sum(pdotq)
            
            # Reshape xi_prime array and add it to the initialized version (one for every quad)
            xi_prime_unnormalized = xi_prime_unnormalized + np.reshape(xi_prime, site_ecef_x.shape)

        # Normalize
        xi_prime_unscaled = xi_prime_unnormalized / n_sub_faults
        
        # Scale so that xi_prime has range (0, 1)
        if mtype == 1: 
            xi_prime = xi_prime_unscaled
        elif mtype == 2:
            xi_prime = 0.5*(xi_prime_unscaled + 1)
        
        self.xi_prime = xi_prime
    

    def computeLD(self):
        """
        Computes LD -- the rupture length de-normalization term. 
        """
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
    
        Lrup_max = 400
        slat = self.sites[1]
        slon = self.sites[0]
        LD = np.zeros_like(slat)
        Ls = np.zeros_like(slat)
        
        ni, nj = LD.shape
        for i in range(ni): 
            for j in range(nj):
                #-----------------------
                # Compute Ls
                #-----------------------
                # Convert to local orthographic
                site_x,site_y = proj(slon[i, j], slat[i, j])
                
                # Shift so center is at epicenter
                site_x2 = site_x - epi_x
                site_y2 = site_y - epi_y
                top_x2 = top_x - epi_x
                top_y2 = top_y - epi_y
                
                # Angle to rotate to put site on x-axis
                alpha = np.arctan2(site_y2, site_x2)
                axis = [0, 0, 1]
                rmat = _rotation_matrix(axis, -alpha)
                llr = [None] * len(top_lat)
                
                # Apply rotation to each point on trace
                for k in range(len(top_lat)):
                    llr[k] = np.dot(rmat, [top_x2[k], top_y2[k], 0])
                    site3 = np.dot(rmat, [site_x2, site_y2, 0])            
                Li = np.min([np.max([a[0] for a in llr]), site3[0]])
                
                #-----------------------------------
                # Compute LD and save results into matrices
                #-----------------------------------
            
                Lrup = np.sqrt(Li*Li + self.Wrup*self.Wrup)
                LD[i, j] = np.log(Lrup)/np.log(Lrup_max)
                Ls[i, j] = Li
        self.LD = LD
        self.Ls = Ls


    def computeDT(self):
        """
        Computes DT -- the distance taper term. 
        """
        slat = self.sites[1]
        slon = self.sites[0]
        site_z = np.zeros_like(slat)
        mesh = geo.mesh.Mesh(slon, slat, site_z)
        Rrup = np.reshape(getDistance('rrup', mesh, self.flt.Quadrilaterals), (-1,))
        nsite = len(Rrup)
        
        if self.simpleDT == True:   # eqn 3.10
            R1 = 35
            R2 = 70
            DT = np.ones(nsite)
            DT[(Rrup > R1) & (Rrup < R2)] = 2 - Rrup[(Rrup > R1) & (Rrup < R2)]/(20 + 10*np.log(self.T))
            DT[Rrup >= R2] = 0
        else:                  # eqn 3.9
            if self.T >= 1:
                R1 = 20 + 10 * np.log(self.T)
                R2 = 2*(20 + 10*np.log(self.T))
            else: 
                R1 = 20
                R2 = 40
            DT = np.ones(nsite)
            DT[(Rrup > R1) & (Rrup < R2)] = 2 - Rrup[(Rrup > R1) & (Rrup < R2)]/35
            DT[Rrup >= R2] = 0
            DT = np.reshape(DT, slat.shape)
        self.DT = DT
        

    def computeWP(self):
        """
        Computes WP -- the narrow-band multiplier. 
        """
        sig = 0.6
        Tp = np.exp(1.27*self.M - 7.28)
        self.WP = np.exp(-(np.log10(self.T/Tp)*np.log10(self.T/Tp))/(2*sig*sig))


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


