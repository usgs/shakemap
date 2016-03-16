#!/usr/bin/env python

#stdlib imports
import struct
from datetime import datetime
import copy
import warnings

#third party imports
from .ecef import latlon2ecef
from .vector import Vector
from .fault import get_quad_length
from openquake.hazardlib.geo import point, geodetic
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.utils import get_orthographic_projection
from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.gsim import base
import numpy as np
import matplotlib.pyplot as plt

#local imports
from shakemap.utils.exception import ShakeMapException


class Distance(object):
    """
    And object to hold distance information. 
    """
    def __init__(self, gmpe, source, lat, lon, dep, use_median_distance = True):
        """
        Construct a Distance object. 
        :param gmpe:
            Concrete subclass of GMPE 
            https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py
        :param source:
            Source object.
        :param lat: 
            A numpy array of site latitudes. 
        :param lon: 
            A numpy array of site longitudes. 
        :param dep: 
            A numpy array of site depths (km). 
        :param use_median_distance: 
            Boolean; only used if GMPE requests fault distances and not fault is 
            availalbe. Default is True, meaning that point-source distances are 
            adjusted based on magnitude to get the median fault distance. 
        :returns:
            Distance object.
        """
        self.source = source
        
        self._distance_context = self._calcDistanceContext(
            gmpe, lat, lon, dep, use_median_distance)
        
        # Place holder for additional sigma due to point-to-fault conversion
        self.delta_sigma = 0.0
    
    @classmethod
    def fromSites(gmpe, source, sites, use_median_distance = True):
        """
        Convenience class method to construct a Distance object from a sites object. 
        :param gmpe:
            Concrete subclass of GMPE 
            https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py
        :param source:
            Source object.
        :param sites: 
            Site object.
        :param use_median_distance: 
            Boolean; only used if GMPE requests fault distances and not fault is 
            availalbe. Default is True, meaning that point-source distances are 
            adjusted based on magnitude to get the median fault distance. 
        :returns:
            Distance object.
        """
        self.source = source
        sm_dict = sites.GeoDict
        west = sm_dict.xmin
        east = sm_dict.xmax
        south = sm_dict.ymin
        north = sm_dict.ymax
        nx = sm_dict.nx
        ny = sm_dict.ny
        lats = np.linspace(north, south, ny)
        lons = np.linspace(west, east, nx)
        lon, lat = np.meshgrid(lons, lats)
        dep = np.zeros_like(lon)
#        mesh = Mesh(lon, lat, dep)
    
    @classmethod
    def fromPoints(gmpe, source, lat, lon, dep, use_median_distance = True):
        self.source = source
        pass
    
    def getDistanceContext(self):
        return self._distance_context
    
    def _calcDistanceContext(self, gmpe, lat, lon, dep, use_median_context = True):
        """
        Create a DistancesContext object
        :param gmpe: 
            Concrete subclass of GMPE 
            (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
        :param lat:
            Numpy array of latitudes. 
        :param lon:
            Numpy array of longitudes. 
        :param dep:
            Numpy array of depths (km). 
        :param use_median_distance: 
            Boolean; only used if GMPE requests fault distances and not fault is 
            availalbe. Default is True, meaning that point-source distances are 
            adjusted based on magnitude to get the median fault distance. 
        :returns:
            DistancesContext object with distance grids required by input gmpe.
        :raises TypeError:
            if gmpe is not a subclass of GMPE
        """
        if not isinstance(gmpe, GMPE):
            raise TypeError('getDistanceContext() cannot work with objects of type "%s"' % type(gmpe))
        
        if self.source.Fault is not None:
            quadlist = self.source.Fault.getQuadrilaterals()
        else:
            quadlist = None
        
        hyplat = self.source.getEventParam('lat')
        hyplon = self.source.getEventParam('lon')
        hypdepth = self.source.getEventParam('depth')
        hyppoint = point.Point(hyplon, hyplat, hypdepth)
        
        context = base.DistancesContext()
        
        ddict = get_distance(list(gmpe.REQUIRES_DISTANCES), lat, lon, dep,
                             quadlist = quadlist, hypo = hyppoint)
        
        for method in gmpe.REQUIRES_DISTANCES:
            (context.__dict__)[method] = ddict[method]

        return context
    
    ## add class methods for fromPoint, fromSites, fromArray


def get_distance(methods, lat, lon, dep, quadlist = None, hypo = None):
    """
    Calculate distance using any one of a number of distance measures. 
    One of quadlist OR hypo must be specified.
    :param methods:
       List of strings (or just a string) of distances to compute; 
       can include: 'rjb', 'rx', 'rrup', 'ry0', 'repi', 'rhypo'
    :param lat:
       A numpy array of latitudes.
    :param lon:
       A numpy array of longidues.
    :param dep:
       A numpy array of depths (km).
    :param quadlist:
       optional list of quadrilaterals (see Fault.py)
    :param hypo:
       Optional Point object of hypocenter. 
    https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
    :returns:
       dictionary of numpy array of distances, size of lon.shape
    :raises ShakeMapException:
       if a fault distance method is called without quadlist
    :raises NotImplementedError:
       for unknown distance measures or ones not yet implemented.
    """
    
    # Dictionary for holding the distances
    distdict = dict()
    
    if not isinstance(methods, list):
        methods = [methods]
    
    methods_available = set(['rjb', 'rx', 'rrup', 'ry0', 'repi', 'rhypo'])
    if not set(methods).issubset(methods_available):
        raise NotImplementedError('One or more requested distance method is not '\
                                  'valid or is not implemented yet')

    if (lat.shape == lon.shape) and (lat.shape == dep.shape):
        pass
    else:
        raise ShakeMapException('lat, lon, and dep must have the same shape.')
    
    oldshape = lon.shape
    
    if len(oldshape) == 2:
        newshape = (oldshape[0]*oldshape[1],1)
    else:
        newshape = (oldshape[0],1)
    
    # Need to integrate ps2ff around here...
    
    if ('rrup' in methods) or ('rjb' in methods):
        x, y, z = latlon2ecef(lat, lon, dep)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))
    
    # ---------------------------------------------
    # Distances that do not require loop over quads
    # ---------------------------------------------
    
    if ('repi' in methods) or \
       (('rjb' in methods) and (quadlist is None)) or \
       (('ry0' in methods) and (quadlist is None)) or \
       (('rx' in methods) and (quadlist is None)):
        if hypo is None:
            raise ShakeMapException('Cannot calculate epicentral distance '\
                                    'without a point object')
        epi = point.Point(hypo.longitude, hypo.latitude, 0.0)
        mesh = Mesh(lon, lat, dep)
        repidist = epi.distance_to_mesh(mesh)
        repidist = repidist.reshape(oldshape)
        distdict['repi'] = repidist
    
    if ('rhypo' in methods) or \
       (('rrup' in methods) and (quadlist is None)):
        if hypo is None:
            raise ShakeMapException('Cannot calculate epicentral distance '\
                                    'without a point object')
        hyppoint = point.Point(hypo.longitude, hypo.latitude, hypo.depth)
        mesh = Mesh(lon, lat, dep)
        rhypodist = hyppoint.distance_to_mesh(mesh)
        rhypodist = rhypodist.reshape(oldshape)
        distdict['rhypo'] = rhypodist

    # ry0 does not require a loop, but it does require quads exist
    if 'ry0' in methods:
        if quadlist is not None:
            # This is just a bandaid on ry0 to get reasonably sensible results
            # for multiple segments ruptures. The 'correct' way to do this is to
            # compare the source-to-site azimuth and use Ry0 = Rx/abs(azimuth)
            # as given in Abrahamson et al. (2014).
            # Note that Ry0 is not the same as:
            #    Ry (in the database summary paper), or
            #    Ry2 that is given in the NGA W2 flatfile, or
            #    U in the GC2 coordinate system.
            # Note: this assumes that the quadrilaters are sorted in order along
            #       strike.
            FP0, FP1, FP2, FP3 = quadlist[0]
            LP0, LP1, LP2, LP3 = quadlist[-1]
            ry0dist = calc_ry0_distance(FP0, LP1, lat, lon, dep)
            distdict['ry0'] = ry0dist
        else:
            warnings.warn('No fault; Replacing ry0 with repi')
            distdict['ry0'] = distdict['repi']
    
    # --------------------------------------------------------
    # Loop over quadlist for those distances that require loop    
    # --------------------------------------------------------
    if 'rrup' in methods:
        minrrup = np.ones(newshape, dtype = lon.dtype)*1e16
    if 'rjb' in methods:
        minrjb = np.ones(newshape, dtype = lon.dtype)*1e16
    if 'rx' in methods:
        totweight = np.zeros(newshape, dtype = lon.dtype)
        meanrxdist = np.zeros(newshape, dtype = lon.dtype)
    
    if quadlist is not None: 
        for quad in quadlist:
            P0, P1, P2, P3 = quad
            
            if 'rrup' in methods:
                rrupdist = calc_rupture_distance(P0, P1, P2, P3, sites_ecef)
                minrrup = np.minimum(minrrup, rrupdist)
            
            if 'rjb' in methods:
                S0 = copy.deepcopy(P0)
                S1 = copy.deepcopy(P1)
                S2 = copy.deepcopy(P2)
                S3 = copy.deepcopy(P3)
                S0.depth = 0.0
                S1.depth = 0.0
                S2.depth = 0.0
                S3.depth = 0.0
                rjbdist = calc_rupture_distance(S0, S1, S2, S3, sites_ecef)
                minrjb = np.minimum(minrjb, rjbdist)
            
            if 'rx' in methods:
                # This is currently just a bandaid to get resonable and 'graceful'
                # results for multiple segment ruptures. The 'correct' wayt to do
                # this is to use the unpublished equations for GC2 (2nd version of
                # Paul Spudich's generalized coordinate system). Rx appears to be
                # defined the same as T in GC2, but gives slightly different values
                # in the NGA W2 flatfile. 
                
                # Compute u_i and t_i for this segment
                t_i = calc_rx_distance(P0, P1, lat, lon)
#                u_i = calc_u_i()
                point0 = point.Point(P0.longitude, P0.latitude, 0.0)
                point1 = point.Point(P1.longitude, P1.latitude, 0.0)
                mesh = Mesh(lon, lat, dep)
                r0 = point0.distance_to_mesh(mesh)
                r1 = point1.distance_to_mesh(mesh)
                li = get_quad_length(quad)
                dweight = (0.5*(1.0/(r0**2) + 1.0/(r1**2)) * li).reshape(newshape)
                totweight = totweight + dweight
                rxdist = t_i
                meanrxdist = meanrxdist + rxdist*dweight
                
        # Collect distances from loop into the distance dict
        if 'rjb' in methods:
            minrjb = minrjb.reshape(oldshape)
            distdict['rjb'] = minrjb
        
        if 'rx' in methods:
            # normalize by sum of quad weights
            meanrxdist = meanrxdist/totweight
            meanrxdist = meanrxdist.reshape(oldshape)
            distdict['rx'] = meanrxdist
            
        if 'rrup' in methods:
            minrrup = minrrup.reshape(oldshape)
            distdict['rrup'] = minrrup
    
    else:
        if 'rjb' in methods:
            warnings.warn('No fault; Replacing rjb with repi')
            distdict['rjb'] = distdict['repi']
        if 'rrup' in methods:
            warnings.warn('No fault; Replacing rrup with rhypo')
            distdict['rrup'] = distdict['rhypo']
        if 'rx' in methods:
            warnings.warn('No fault; Replacing rx with repi')
            distdict['rx'] = distdict['repi']
    
    return distdict

def distance_sq_to_segment(p0, p1):
    """
    Calculate the distance^2 from the origin to a segment defined by two vectors
    :param p0: numpy array Nx3 (ECEF)
    :param p1: numpy array Nx3 (ECEF)
    :returns:
        The squared distance from the origin to a segment.
    """
    # /*
    #  * This algorithm is from (Vince's) CS1 class.
    #  * It returns the distance^2 from the origin to a segment defined
    #  * by two vectors
    #  */
    
    dist = np.zeros_like(p1[:,0])
    # /* Are the two points equal? */
    idx_equal = (p0[:,0] == p1[:,0]) & (p0[:,1] == p1[:,1]) & (p0[:,2] == p1[:,2])
    dist[idx_equal] = np.sqrt(p0[idx_equal,0]**2+p0[idx_equal,1]**2+p0[idx_equal,2]**2)
    
    v = p1 - p0
    
    # /*
    #  * C1 = c1/|v| is the projection of the origin O on line (P0,P1).
    #  * If C1 is negative, then O is outside the segment and
    #  * closer to the P0 side.
    #  * If C1 is positive and >V then O is on the other side.
    #  * If C1 is positive and <V then O is inside.
    #  */
    
    c1 = -1 * np.sum(p0*v,axis=1)
    idx_neg = c1 <= 0
    dist[idx_neg] = p0[idx_neg,0]**2+p0[idx_neg,1]**2+p0[idx_neg,2]**2
    
    c2 = np.sum(v*v,axis=1)
    idx_less_c1 = c2 <= c1
    dist[idx_less_c1] = p1[idx_less_c1,0]**2+p1[idx_less_c1,1]**2+p1[idx_less_c1,2]**2
    
    idx_other = np.logical_not(idx_neg | idx_equal | idx_less_c1)

    nr,nc = p0.shape
    t1 = c1/c2
    t1.shape = (nr, 1)
    t2 = np.tile(t1, 3)
    tmp= p0 + (v*t1)
    dist[idx_other] = tmp[idx_other, 0]**2 + tmp[idx_other, 1]**2 + tmp[idx_other, 2]**2
    
    return dist

#call this once per quad
def calc_rupture_distance(P0, P1, P2, P3, points):
    """
    Calculate the shortest distance from a set of points to a rupture surface.
    :param P0: 
        Point object, representing the first top-edge vertex of a fault quadrilateral.
    :param P1: 
        Point object, representing the second top-edge vertex of a fault quadrilateral.
    :param P2: 
        Point object, representing the first bottom-edge vertex of a fault quadrilateral.
    :param P3: 
        Point object, representing the second bottom-edge vertex of a fault quadrilateral.
    :param points: 
        Numpy array Nx3 of points (ECEF) to calculate distance from.
    :returns:
        Array of size N of distances (in km) from input points to rupture surface.
    """
    # Convert to ecef
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)
    
    # Make a unit vector normal to the plane
    normalVector = (p1 - p0).cross(p2 - p0).norm()

    dist = np.ones_like(points[:, 0])*np.nan

    p0d = p0.getArray() - points
    p1d = p1.getArray() - points
    p2d = p2.getArray() - points
    p3d = p3.getArray() - points
    
    # Create 4 planes with normals pointing outside rectangle
    n0 = (p1-p0).cross(normalVector).getArray()
    n1 = (p2-p1).cross(normalVector).getArray()
    n2 = (p3-p2).cross(normalVector).getArray()
    n3 = (p0-p3).cross(normalVector).getArray()
    
    sgn0 = np.signbit(np.sum(n0*p0d, axis = 1))
    sgn1 = np.signbit(np.sum(n1*p1d, axis = 1))
    sgn2 = np.signbit(np.sum(n2*p2d, axis = 1))
    sgn3 = np.signbit(np.sum(n3*p3d, axis = 1))
    
    inside_idx = (sgn0 == sgn1) & (sgn1 == sgn2) & (sgn2 == sgn3)
    dist[inside_idx] = np.power(np.abs(
        np.sum(p0d[inside_idx,:] * normalVector.getArray(), axis=1)), 2)

    outside_idx = np.logical_not(inside_idx)
    s0 = distance_sq_to_segment(p0d, p1d)
    s1 = distance_sq_to_segment(p1d, p2d)
    s2 = distance_sq_to_segment(p2d, p3d)
    s3 = distance_sq_to_segment(p3d, p0d)

    smin = np.minimum(np.minimum(s0,s1),np.minimum(s2, s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist)/1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0],1)
    if np.any(np.isnan(dist)):
        raise ShakeMapException("Could not calculate some distances!")
    dist = np.fliplr(dist)
    return dist

def calc_ry0_distance(P0, P1, lat, lon, dep):
    """Calculate Ry0 distance.

    Compute the minimum distance between sites (lat, lon, dep) and the great
    circle arcs perpendicular to the average strike direction of the
    fault trace and passing through the end-points of the trace.
    
    :param P0:
      Point object, representing the first top-edge vertex of a fault quadrilateral.
    :param P1:
      Point object, representing the second top-edge vertex of a fault quadrilateral.
    :param lat: 
      Numpy array of latitude. 
    :param lon: 
      Numpy array of longitude. 
    :param dep: 
      Numpy array of depths (km). 
    :returns:
      Array of size lon.shape of distances (in km) from input points to rupture surface.
    """
    
    # Strike
    surfaceP0 = point.Point(P0.longitude, P0.latitude, 0.0)
    surfaceP1 = point.Point(P1.longitude, P1.latitude, 0.0)
    strike = P0.azimuth(P1)
    dst1 = geodetic.distance_to_arc(P0.longitude,P0.latitude,
                                    (strike + 90.) % 360, lon, lat)
    dst2 = geodetic.distance_to_arc(P1.longitude,P1.latitude,
                                    (strike + 90.) % 360, lon, lat)
    # Get the shortest distance from the two lines
    idx = np.sign(dst1) == np.sign(dst2)
    dst = np.zeros_like(dst1)
    dst[idx] = np.fmin(np.abs(dst1[idx]), np.abs(dst2[idx]))
    return dst

def calc_u_i(P0, P1, lat, lon, dep):
    """Calculate u_i distance. See Spudich and Chiou OFR 2015-1028. 
    
    :param P0:
      Point object, representing the first top-edge vertex of a fault quadrilateral.
    :param P1:
      Point object, representing the second top-edge vertex of a fault quadrilateral.
    :param lat:
      A numpy array of latitude.
    :param lon:
      A numpy array of longitude.
    :param dep:
      A numpy array of depths (km).
    :returns:
      Array of size mesh.lons of distances (in km) from input points to rupture surface.
    """
    pass

def calc_rx_distance(P0, P1, lat, lon):
    """
    The shortest horizontal distance from a set of points to a line defined by
    extending the fault trace (or the top edge of the rupture) to
    infinity in both directions. Values on the hanging-wall are
    positive and those on the foot-wall are negative.
    :param P0: 
        Point object, representing the first top-edge vertex of a fault quadrilateral.
    :param P1: 
        Point object, representing the second top-edge vertex of a fault quadrilateral.
    :param lat: 
        A numpy array of latitudes. 
    :param lon: 
        A numpy array of longitudes. 
    :returns:
        Array of size N of distances (in km) from input points to rupture surface.
    """
    # Project to Cartesian space:
    west = min(P0.x, P1.x)
    east = max(P0.x, P1.x)
    south = min(P0.y, P1.y)
    north = max(P0.y, P1.y)
    proj = get_orthographic_projection(west, east, north, south)
    
    # projected coordinates are in km
    p0x, p0y = proj(P0.x, P0.y)
    p1x, p1y = proj(P1.x, P1.y)
    
    vhat = Vector(p1y - p0y, -(p1x - p0x), 0)
    vhat = vhat.norm()
    
    # Convert sites to Cartesian
    sx, sy = proj(lon, lat)
    sx1d = np.reshape(sx, (-1,))
    sy1d = np.reshape(sy, (-1,))
    
    r = np.zeros([len(sx1d), 2])
    r[:,0] = sx1d - p0x
    r[:,1] = sy1d - p0y
    dist = np.sum(vhat.getArray()[0:2]*r, axis = 1)
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0], 1)
    dist = np.fliplr(dist)
    return dist

def get_top_edge(lat, lon, dep):
    """
    Determine which edge of a quadrilateral rupture surface is the top.
    :param lat: Sequence of 4 or 5 latitudes defining vertices of rupture surface.
    :param lon: Sequence of 4 or 5 longitudes defining vertices of rupture surface.
    :param dep: Sequence of 4 or 5 depths defining vertices of rupture surface.
    :returns:
        (P0,P1) where P0 and P1 are Vector objects in ECEF coordinates indicating 
        the vertices of the top edge.
    """
    lat = np.array(lat)
    lon = np.array(lon)
    dep = np.array(dep)
    p1 = None
    p2 = None
    if sum(np.diff(dep)) == 0:
        p1 = Vector.fromPoint(point.Point(lon[0], lat[0], dep[0]))
        p2 = Vector.fromPoint(point.Point(lon[1], lat[1], dep[1]))
    else:
        dep2 = dep[0:4]
        dd = np.diff(dep2)
        idx = dd == 0
        idx = np.append(idx,False)
        p1idx = dep2[idx].argmin()+1
        p2idx = p1idx + 1
        p1 = Vector.fromPoint(point.Point(lon[p1idx], lat[p1idx],0.0))
        p2 = Vector.fromPoint(point.Point(lon[p2idx], lat[p2idx],0.0))
    return (p1, p2)

def minimize(a,b):
    """
    Get minimum values from two numpy arrays, ignoring sign in comparison 
    but preserving sign in result.
    :param a:
      Numpy array
    :param b:
      Numpy array
    :returns:
      Numpy array
    """
    d = np.empty_like(a)
    aidx = np.abs(a) < np.abs(b)
    bidx = np.logical_not(aidx)
    d[aidx] = a[aidx]
    d[bidx] = b[bidx]
    return d
    
