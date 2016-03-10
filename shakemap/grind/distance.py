#!/usr/bin/env python

#stdlib imports
import struct
from datetime import datetime
import copy

#third party imports
from .ecef import latlon2ecef
from .vector import Vector
from openquake.hazardlib.geo import point
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo.utils import get_orthographic_projection
import numpy as np
import matplotlib.pyplot as plt

#local imports
from shakemap.utils.exception import ShakeMapException

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
    
def get_distance(method,mesh,quadlist=None,mypoint=None):
    """
    Calculate distance using any one of a number of distance measures. 
    One of quadlist OR mypoint must be specified.
    :param method:
       One of: 'rjb','rx','rrup','ry0','rcdbp','repi','rhypo'
    :param mesh:
       A Mesh object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
    :param quadlist:
       optional list of quadrilaterals (see Fault.py)
    :param mypoint:
       optional Point object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py)
    :returns:
       numpy array of distances, size of mesh.lons
    :raises ShakeMapException:
       if a fault distance method is called without quadlist
    :raises NotImplementedError:
       for unknown distance measures or ones not yet implemented.
    """
    oldshape = mesh.lons.shape
    if len(oldshape) == 2:
        newshape = (oldshape[0]*oldshape[1],1)
    else:
        newshape = (oldshape[0],1)
    if (method == 'rjb') or (method == 'rrup') or (method == 'rx') or (method == 'ry0'):
        if quadlist is None:
            raise ShakeMapException('Cannot calculate rupture distance %s without a list of quadrilaterals' % method)
        mindist = np.ones(newshape,dtype=mesh.lons.dtype)*1e16
        x,y,z = latlon2ecef(mesh.lats,mesh.lons,mesh.depths)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        points = np.hstack((x,y,z))

    if method == 'rjb':
        for quad in quadlist:
            P0,P1,P2,P3 = quad
            P0 = copy.deepcopy(P0)
            P1 = copy.deepcopy(P1)
            P2 = copy.deepcopy(P2)
            P3 = copy.deepcopy(P3)
            P0.depth = 0.0
            P1.depth = 0.0
            P2.depth = 0.0
            P3.depth = 0.0
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            rjbdist = calc_rupture_distance(p0,p1,p2,p3,points)
            mindist = np.minimum(mindist,rjbdist)
        mindist = mindist.reshape(oldshape)
        return mindist
    elif method == 'rx':
        # This is currently just a bandaid to get resonable and 'graceful'
        # results for multiple segment ruptures. The 'correct' wayt to do
        # this is to use the unpublished equations for GC2 (2nd version of
        # Paul Spudich's generalized coordinate system). Rx appears to be
        # defined the same as T in GC2, but gives slightly different values
        # in the NGA W2 flatfile. 
        meanrxdist = np.zeros(newshape, dtype=mesh.lons.dtype)
        totweight = np.zeros(newshape, dtype=mesh.lons.dtype)
        for quad in quadlist:
            # Points on top edge of quad
            P0, P1, P2, P3 = quad
            
            # Project these two points into a Cartesian space
            west = min(P0.x, P1.x)
            east = max(P0.x, P1.x)
            south = min(P0.y, P1.y)
            north = max(P0.y, P1.y)
            proj = get_orthographic_projection(west, east, north, south)
            # projected coordinates are in km
            p0x, p0y = proj(P0.x, P0.y)
            p1x, p1y = proj(P1.x, P1.y)
            
            # Convert from m to km
            PP0 = Vector(p0x*1000, p0y*1000, 0.0)
            PP1 = Vector(p1x*1000, p1y*1000, 0.0)
            
            # Project sites
            ppointx, ppointy = proj(mesh.lons.flatten(), mesh.lats.flatten())
            ppoints = np.zeros((len(ppointx), 3))
            
            # Convert coordinates to meters and put in numpy array
            ppoints[:, 0] = ppointx*1000 
            ppoints[:, 1] = ppointy*1000
            
            # Compute distance to "segment"
            # Hacky... use Rjb function and put bottom points below top trace.
            S0 = copy.deepcopy(P0)
            S1 = copy.deepcopy(P1)
            S2 = copy.deepcopy(P1)
            S3 = copy.deepcopy(P0)
            S0.depth = 0.0
            S1.depth = 0.0
            S2.depth = 1.0
            S3.depth = 1.0
            s0 = Vector.fromPoint(S0)
            s1 = Vector.fromPoint(S1)
            s2 = Vector.fromPoint(S2)
            s3 = Vector.fromPoint(S3)
            topdist = calc_rupture_distance(s0, s1, s2, s3, points)
            
            # Weight of Rx is inverse squared distance
            dweight = 1.0/(topdist**2)
            totweight = totweight + dweight
            rxdist = calc_rx_distance(PP0, PP1, ppoints)
            meanrxdist = meanrxdist + rxdist*dweight
        meanrxdist = meanrxdist/totweight
        meanrxdist = meanrxdist.reshape(oldshape)
        return meanrxdist
    
    elif method == 'rrup':
        for quad in quadlist:
            P0,P1,P2,P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            rrupdist = calc_rupture_distance(p0,p1,p2,p3,points)
            mindist = np.minimum(mindist,rrupdist)
        mindist = mindist.reshape(oldshape)
        return mindist
    elif method == 'ry0':
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
        ry0dist = calc_ry0_distance(FP0, LP1, mesh)
        return ry0dist
    elif method == 'rcdbp':
        raise NotImplementedError('rcdbp distance measure is not implemented yet')
    elif method == 'repi':
        if mypoint is None:
            raise ShakeMapException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(mypoint.longitude,mypoint.latitude,0.0)
        mindist = newpoint.distance_to_mesh(mesh)
        mindist = mindist.reshape(oldshape)
        return mindist
    elif method == 'rhypo':
        if mypoint is None:
            raise ShakeMapException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(mypoint.longitude,mypoint.latitude,mypoint.depth)
        mindist = newpoint.distance_to_mesh(mesh)
        mindist = mindist.reshape(oldshape)
        return mindist
    else:
        raise NotImplementedError('"%s" distance measure is not valid or is not implemented yet' % method)

def get_distance2(methods, mesh, quadlist = None, mypoint = None):
    """
    Calculate distance using any one of a number of distance measures. 
    One of quadlist OR mypoint must be specified.
    :param methods:
       List of strings (or just a string) of distances to compute; can include: 'rjb', 'rx', 'rrup', 'ry0', 'repi', 'rhypo'
    :param mesh:
       A Mesh object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
    :param quadlist:
       optional list of quadrilaterals (see Fault.py)
    :param mypoint:
       optional Point object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py)
    :returns:
       dictionary of numpy array of distances, size of mesh.lons
    :raises ShakeMapException:
       if a fault distance method is called without quadlist
    :raises NotImplementedError:
       for unknown distance measures or ones not yet implemented.
    """
    if not isinstance(methods, list):
        methods = [methods]
    
    oldshape = mesh.lons.shape
    if len(oldshape) == 2:
        newshape = (oldshape[0]*oldshape[1],1)
    else:
        newshape = (oldshape[0],1)
    
    if quadlist is None:
        raise ShakeMapException('Cannot calculate rupture distance %s without a list of quadrilaterals' % method)
    
    if any([a == 'rrup' for a in methods]) or \
       any([a == 'rx' for a in methods]) or \
       any([a == 'rjb' for a in methods]):
       x,y,z = latlon2ecef(mesh.lats,mesh.lons,mesh.depths)
       x.shape = newshape
       y.shape = newshape
       z.shape = newshape
       points = np.hstack((x,y,z))
    
    
    # --------------------------------------------------------
    # Loop over quadlist for those distances that require loop    
    # --------------------------------------------------------
    if any([a == 'rrup' for a in methods]):
        minrrup = np.ones(newshape, dtype = mesh.lons.dtype)*1e16
    if any([a == 'rjb' for a in methods]):
        minrjb = np.ones(newshape, dtype=mesh.lons.dtype)*1e16
    if any([a == 'rx' for a in methods]):
        totweight = np.zeros(newshape, dtype=mesh.lons.dtype)
        meanrxdist = np.zeros(newshape, dtype=mesh.lons.dtype)
    
    
    for quad in quadlist:
        P0, P1, P2, P3 = quad
        
        if any([a == 'rrup' for a in methods]):
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            rrupdist = calc_rupture_distance(p0, p1, p2, p3, points)
            minrrup = np.minimum(minrrup, rrupdist)
        
        if any([a == 'rjb' for a in methods]):
            P0.depth = 0.0
            P1.depth = 0.0
            P2.depth = 0.0
            P3.depth = 0.0
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            rjbdist = calc_rupture_distance(p0, p1, p2, p3, points)
            minrjb = np.minimum(minrjb, rjbdist)
        
        if any([a == 'rx' for a in methods]):
            # This is currently just a bandaid to get resonable and 'graceful'
            # results for multiple segment ruptures. The 'correct' wayt to do
            # this is to use the unpublished equations for GC2 (2nd version of
            # Paul Spudich's generalized coordinate system). Rx appears to be
            # defined the same as T in GC2, but gives slightly different values
            # in the NGA W2 flatfile. 
            
            # Project top two points into a Cartesian space
            west = min(P0.x, P1.x)
            east = max(P0.x, P1.x)
            south = min(P0.y, P1.y)
            north = max(P0.y, P1.y)
            proj = get_orthographic_projection(west, east, north, south)
            # projected coordinates are in km
            p0x, p0y = proj(P0.x, P0.y)
            p1x, p1y = proj(P1.x, P1.y)
            
            # Convert from m to km
            PP0 = Vector(p0x*1000, p0y*1000, 0.0)
            PP1 = Vector(p1x*1000, p1y*1000, 0.0)
            
            # Project sites
            ppointx, ppointy = proj(mesh.lons.flatten(), mesh.lats.flatten())
            ppoints = np.zeros((len(ppointx), 3))
            
            # Convert coordinates to meters and put in numpy array
            ppoints[:, 0] = ppointx*1000 
            ppoints[:, 1] = ppointy*1000
            
            # Compute distance to "segment"
            # Hacky... use Rjb function and put bottom points below top trace.
            S0 = copy.deepcopy(P0)
            S1 = copy.deepcopy(P1)
            S2 = copy.deepcopy(P1)
            S3 = copy.deepcopy(P0)
            S0.depth = 0.0
            S1.depth = 0.0
            S2.depth = 1.0
            S3.depth = 1.0
            s0 = Vector.fromPoint(S0)
            s1 = Vector.fromPoint(S1)
            s2 = Vector.fromPoint(S2)
            s3 = Vector.fromPoint(S3)
            topdist = calc_rupture_distance(s0, s1, s2, s3, points)
            
            # Weight of Rx is inverse squared distance
            dweight = 1.0/(topdist**2)
            totweight = totweight + dweight
            rxdist = calc_rx_distance(PP0, PP1, ppoints)
            meanrxdist = meanrxdist + rxdist*dweight
                
    # Collect distances from loop into a dict
    distdict = dict()
    if any([a == 'rjb' for a in methods]):
        minrjb = minrjb.reshape(oldshape)
        distdict['rjb'] = minrjb
    
    if any([a == 'rx' for a in methods]):
        # normalize by sum of quad weights
        meanrxdist = meanrxdist/totweight
        meanrxdist = meanrxdist.reshape(oldshape)
        distdict['rx'] = meanrxdist
    
    if any([a == 'rrup' for a in methods]):
        minrrup = minrrup.reshape(oldshape)
        distdict['rrup'] = minrrup
    
    # -------------------------------------------------------
    # Remaining distances that do not require loop over quads
    # -------------------------------------------------------
    if any([a == 'ry0' for a in methods]):
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
        ry0dist = calc_ry0_distance(FP0, LP1, mesh)
        distdict['ry0'] = ry0dist
    
    if any([a == 'repi' for a in methods]):
        if mypoint is None:
            raise ShakeMapException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(mypoint.longitude, mypoint.latitude, 0.0)
        repidist = newpoint.distance_to_mesh(mesh)
        repidist = repidist.reshape(oldshape)
        distdict['repi'] = repidist
    
    if any([a == 'rhypo' for a in methods]):
        if mypoint is None:
            raise ShakeMapException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(mypoint.longitude, mypoint.latitude, mypoint.depth)
        rhypodist = newpoint.distance_to_mesh(mesh)
        rhypodist = rhypdist.reshape(oldshape)
        distdict['rhypo'] = rhypodist

#    if 
#        raise NotImplementedError('"%s" distance measure is not valid or is not implemented yet' % method)
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
    t1.shape = (nr,1)
    t2 = np.tile(t1,3)
    tmp= p0+(v*t1)
    dist[idx_other] = tmp[idx_other,0]**2 + tmp[idx_other,1]**2 + tmp[idx_other,2]**2
    
    return dist

#call this once per quad
def calc_rupture_distance(P0,P1,P2,P3,points):
    """
    Calculate the shortest distance from a set of points to a rupture surface.
    :param P0: Vector object (in ECEF coordinates) defining the first vertex of a fault plane.
    :param P1: Vector object (in ECEF coordinates) defining the second vertex of a fault plane.
    :param P2: Vector object (in ECEF coordinates) defining the third vertex of a fault plane.
    :param P3: Vector object (in ECEF coordinates) defining the fourth vertex of a fault plane.
    :param points: numpy array Nx3 of points (ECEF) to calculate distance from.
    :returns:
        Array of size N of distances (in km) from input points to rupture surface.
    """
    #Make a unit vector normal to the plane
    normalVector = (P1-P0).cross(P2-P0).norm()

    dist = np.ones_like(points[:,0])*np.nan

    p0d = P0.getArray()-points
    p1d = P1.getArray()-points
    p2d = P2.getArray()-points
    p3d = P3.getArray()-points
    
    #Create 4 planes with normals pointing outside rectangle
    n0 = (P1-P0).cross(normalVector).getArray()
    n1 = (P2-P1).cross(normalVector).getArray()
    n2 = (P3-P2).cross(normalVector).getArray()
    n3 = (P0-P3).cross(normalVector).getArray()
    
    sgn0 = np.signbit(np.sum(n0*p0d,axis=1))
    sgn1 = np.signbit(np.sum(n1*p1d,axis=1))
    sgn2 = np.signbit(np.sum(n2*p2d,axis=1))
    sgn3 = np.signbit(np.sum(n3*p3d,axis=1))
    
    inside_idx = (sgn0 == sgn1) & (sgn1 == sgn2) & (sgn2 == sgn3)
    dist[inside_idx] = np.power(np.abs(np.sum(p0d[inside_idx,:] * normalVector.getArray(),axis=1)),2)

    outside_idx = np.logical_not(inside_idx)
    s0 = distance_sq_to_segment(p0d, p1d)
    s1 = distance_sq_to_segment(p1d, p2d)
    s2 = distance_sq_to_segment(p2d, p3d)
    s3 = distance_sq_to_segment(p3d, p0d)

    smin = np.minimum(np.minimum(s0,s1),np.minimum(s2,s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist)/1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0],1)
    if np.any(np.isnan(dist)):
        raise ShakeMapException("Could not calculate some distances!")
    dist = np.fliplr(dist)
    return dist

def calc_ry0_distance(P0,P1,mesh):
    """Calculate Ry0 distance.

    Compute the minimum distance between each point of a mesh and the great
    circle arcs perpendicular to the average strike direction of the
    fault trace and passing through the end-points of the trace.
    
    :param P0:
      Point object, representing the first top-edge vertex of a fault quadrilateral.
    :param P1:
      Point object, representing the second top-edge vertex of a fault quadrilateral.
    :returns:
      Array of size mesh.lons of distances (in km) from input points to rupture surface.
    """
    #get the mean strike vector
    surfaceP0 = point.Point(P0.longitude,P0.latitude,0.0)
    surfaceP1 = point.Point(P1.longitude,P1.latitude,0.0)
    strike = P0.azimuth(P1)
    dst1 = geodetic.distance_to_arc(P0.longitude,P0.latitude,(strike + 90.) % 360,mesh.lons, mesh.lats)
    dst2 = geodetic.distance_to_arc(P1.longitude,P1.latitude,(strike + 90.) % 360,mesh.lons, mesh.lats)
    # Get the shortest distance from the two lines
    idx = np.sign(dst1) == np.sign(dst2)
    dst = np.zeros_like(dst1)
    dst[idx] = np.fmin(np.abs(dst1[idx]), np.abs(dst2[idx]))
    return dst

def calc_rx_distance(P0,P1,points):
    """
    The shortest horizontal distance from a set of points to a line defined by
    extending the fault trace (or the top edge of the rupture) to
    infinity in both directions. Values on the hanging-wall are
    positive and those on the foot-wall are negative.
    :param P0: Vector object (in ECEF coordinates) defining the first vertex of the top edge of rupture surface.
    :param P1: Vector object (in ECEF coordinates) defining the second vertex of the top edge of rupture surface.
    :param points: numpy array Nx3 of points (ECEF) to calculate distance from.
    :returns:
        Array of size N of distances (in km) from input points to rupture surface.
    """
    x1 = P0.x
    x2 = P1.x
    y1 = P0.y
    y2 = P1.y
    vhat = Vector(y2-y1,-(x2-x1),0)
    vhat = vhat.norm()
    r = np.zeros_like(points[:,0:2])
    r[:,0] = points[:,0] - x1
    r[:,1] = points[:,1] - y1
    dist = np.sum(vhat.getArray()[0:2]*r,axis=1)/1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0],1)
    dist = np.fliplr(dist)
    return dist

def get_top_edge(lat,lon,dep):
    """
    Determine which edge of a quadrilateral rupture surface is the top.
    :param lat: Sequence of 4 or 5 latitudes defining vertices of rupture surface.
    :param lon: Sequence of 4 or 5 longitudes defining vertices of rupture surface.
    :param dep: Sequence of 4 or 5 depths defining vertices of rupture surface.
    :returns:
        (P0,P1) where P0 and P1 are Vector objects in ECEF coordinates indicating the vertices of the top edge.
    """
    lat = np.array(lat)
    lon = np.array(lon)
    dep = np.array(dep)
    p1 = None
    p2 = None
    if sum(np.diff(dep)) == 0:
        p1 = Vector.fromPoint(point.Point(lon[0],lat[0],dep[0]))
        p2 = Vector.fromPoint(point.Point(lon[1],lat[1],dep[1]))
    else:
        dep2 = dep[0:4]
        dd = np.diff(dep2)
        idx = dd == 0
        idx = np.append(idx,False)
        p1idx = dep2[idx].argmin()+1
        p2idx = p1idx + 1
        p1 = Vector.fromPoint(point.Point(lon[p1idx],lat[p1idx],0.0))
        p2 = Vector.fromPoint(point.Point(lon[p2idx],lat[p2idx],0.0))
    return (p1,p2)

