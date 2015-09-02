#!/usr/bin/env python

#stdlib imports
import struct
from datetime import datetime

#third party imports
from ecef import latlon2ecef
from vector import Vector
from openquake.hazardlib.geo import point
import numpy as np
import matplotlib.pyplot as plt
from neicio.cmdoutput import getCommandOutput

class DistanceException(Exception):
    """
    Used to indicate errors specific to the various methods in this distance.py.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def minimize(a,b):
    """
    Get minimum values from two numpy arrays, ignoring sign in comparison but preserving sign in result.
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
    
def getDistance(method,mesh,quadlist=None,point=None):
    """
    Calculate distance using any one of a number of distance measures. One of quadlist OR point must be specified.
    :param method:
       One of: 'rjb','rx','rrup','ry0','rcdbp','epi','hypo'
    :param mesh:
       A Mesh object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
    :param quadlist:
       optional list of quadrilaterals (see Fault.py)
    :param point:
       optional Point object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py)
    :returns:
       numpy array of distances, size of mesh.lons
    :raises DistanceException:
       if a fault distance method is called without quadlist
    :raises NotImplementedError:
       for unknown distance measures or ones not yet implemented.
    """
    if method == 'rjb':
        raise NotImplementedError('rjb distance measure is not implemented yet')
    if method == 'rx':
        if quadlist is None:
            raise DistanceException('Cannot calculate rx distance without a list of quadrilaterals')
        mindist = np.zeros_like(mesh.lons)
        for quad in quadlist:
            P0,P1 = quad[0:2]
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            x,y,z = ecef.latlon2ecef(mesh.lats,mesh.lons,mesh.depths)
            points = np.vstack(x.flatten(),y.flatten(),z.flatten())
            rxdist = calcRxDistance(p0,p1,points)
            mindist = minimize(mindist,rxdist)
        return mindist
    if method == 'rrup':
        if quadlist is None:
            raise DistanceException('Cannot calculate rupture distance without a list of quadrilaterals')
        mindist = np.zeros_like(mesh.lons)
        for quad in quadlist:
            P0,P1,P2,P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            x,y,z = ecef.latlon2ecef(mesh.lats,mesh.lons,mesh.depths)
            points = np.vstack(x.flatten(),y.flatten(),z.flatten())
            rrupdist = calcRuptureDistance(p0,p1,p2,p3,points)
            mindist = np.minimum(mindist,rrupdist)
        return mindist
    if method == 'ry0':
        raise NotImplementedError('ry0 distance measure is not implemented yet')
    if method == 'rcdpp':
        raise NotImplementedError('rcdbp distance measure is not implemented yet')
    if method == 'epi':
        if point is None:
            raise DistanceException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(point.latitude,point.longitude,0.0)
        return newpoint.distance_to_mesh(mesh) 
    if method == 'hypo':
        if point is None:
            raise DistanceException('Cannot calculate epicentral distance without a point object')
        newpoint = point.Point(point.latitude,point.longitude,point.depth)
        return newpoint.distance_to_mesh(mesh)
    else:
        raise NotImplementedError('"%s" distance measure is not valid or is not implemented yet' % method)

def dist2segment(p0, p1):
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
def calcRuptureDistance(P0,P1,P2,P3,points):
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
    s0 = dist2segment(p0d,p1d)
    s1 = dist2segment(p1d,p2d)
    s2 = dist2segment(p2d,p3d)
    s3 = dist2segment(p3d,p0d)

    smin = np.minimum(np.minimum(s0,s1),np.minimum(s2,s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist)/1000.0

    if np.any(np.isnan(dist)):
        raise Exception,"Could not calculate some distances!"
    return dist

def calcRxDistance(P0,P1,points):
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
    dist = np.sum(vhat.getArray()[0:2]*r,axis=1)
    return dist

def getTopEdge(lat,lon,dep):
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
        p1 = Vector.fromPoint(point.Point(lat[0],lon[0],dep[0]))
        p2 = Vector.fromPoint(point.Point(lat[1],lon[1],dep[1]))
    else:
        dep2 = dep[0:4]
        dd = np.diff(dep2)
        idx = dd == 0
        idx = np.append(idx,False)
        p1idx = dep2[idx].argmin()+1
        p2idx = p1idx + 1
        p1 = Vector.fromPoint(point.Point(lat[p1idx],lon[p1idx],0.0))
        p2 = Vector.fromPoint(point.Point(lat[p2idx],lon[p2idx],0.0))
    return (p1,p2)

if __name__ == '__main__':
    flat = [28.427,27.986,27.469,27.956]
    flon = [84.614,86.179,85.880,84.461]
    fdep = [20.0,20.0,13.0,13.0]
    lat0,lat1,lat2,lat3 = flat
    lon0,lon1,lon2,lon3 = flon
    dep0,dep1,dep2,dep3 = fdep
    
    # lat0,lat1,lat2,lat3 = [28.0,28.0,26.0,26.0]
    # lon0,lon1,lon2,lon3 = [84.0,86.0,86.0,84.0]
    # dep0,dep1,dep2,dep3 = [10.0,10.0,10.0,10.0]

    P0 = Vector.fromPoint(point.Point(lat0,lon0,dep0))
    P1 = Vector.fromPoint(point.Point(lat1,lon1,dep1))
    P2 = Vector.fromPoint(point.Point(lat2,lon2,dep2))
    P3 = Vector.fromPoint(point.Point(lat3,lon3,dep3))

    lons = np.arange(81.5,87.5,0.0083)
    lats = np.arange(25.5,31.0,0.0083)
    lonmat,latmat = np.meshgrid(lons,lats)
    depmat = np.zeros_like(lonmat)

    # lons = np.arange(83.0,88.0,1.0)
    # lats = np.arange(25.0,30.0,1.0)
    # lonmat,latmat = np.meshgrid(lons,lats)
    # depmat = np.zeros_like(lonmat)

    # lonmat = np.array([[25.0]])
    # latmat = np.array([[85.0]])
    # depmat = np.array([[0.0]])
    
    x,y,z = latlon2ecef(lonmat,latmat,depmat)
    nr,nc = x.shape
    x.shape = (nr*nc,1)
    y.shape = (nr*nc,1)
    z.shape = (nr*nc,1)
    points = np.hstack((x,y,z))
    t1 = datetime.now()
    dist = calcRuptureDistance(P0,P1,P2,P3,points)
    t2 = datetime.now()
    dt = t2-t1
    sec = dt.seconds + float(dt.microseconds)/1e6
    npoints = nr*nc
    print 'RRupt: %.4f seconds for %i point grid' % (sec,npoints)
    dist.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.imshow(dist,interpolation='nearest')
        plt.colorbar()
        plt.savefig('dist.png')
    else:
        print dist

    #test this against Bruce's C program
    nr2,nc2 = points.shape
    f = open('points.bin','wb')
    for i in range(0,nr2):
        xt,yt,zt = points[i,:]
        xb = struct.pack('d',xt)
        yb = struct.pack('d',yt)
        zb = struct.pack('d',zt)
        f.write(xb)
        f.write(yb)
        f.write(zb)
    f.close()

    #now write out the quads
    f = open('fault.txt','wt')
    f.write('5\n')
    f.write('%.11f %.11f %.11f\n' % tuple(P0.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P1.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P2.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P3.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P0.getArray()))
    f.close()

    cmd = '/Users/mhearne/src/python/ShakeMap/src/contour/dist_rrupt fault.txt < points.bin > outpoints.bin'
    res,stdout,stderr = getCommandOutput(cmd)
    f = open('outpoints.bin','rb')
    dlist = []
    buf = f.read(8)
    while buf:
        d = struct.unpack('d',buf)
        dlist.append(d[0])
        buf = f.read(8)
    f.close()
    d2 = np.array(dlist)
    d2.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(d2,interpolation='nearest')
        plt.colorbar()
        plt.savefig('dist_c.png')
    else:
        print dist
    
    P0,P1 = getTopEdge(flat,flon,fdep)
    rxdist = calcRxDistance(P0,P1,points)
    rxdist.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(rxdist,interpolation='nearest')
        plt.colorbar()
        plt.savefig('distrx.png')
    else:
        print rxdist

    #test against Bruce's C implementation
    f = open('points.txt','wt')
    for i in range(0,nr2):
        xt,yt = (points[i,0],points[i,1])
        f.write('%.11f %.11f\n' %(xt,yt))
    f.close()
    cmd = '/Users/mhearne/src/python/ShakeMap/src/contour/dist_rx %.11f %.11f %.11f %.11f < points.txt > outpoints.txt' % (P0.x,P0.y,P1.x,P1.y)
    res,stdout,stderr = getCommandOutput(cmd)
    rxdist2 = []
    for line in open('outpoints.txt','rt').readlines():
        rxdist2.append(float(line.strip()))
    rxdist2 = np.array(rxdist2)
    rxdist2.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(rxdist2,interpolation='nearest')
        plt.colorbar()
        plt.savefig('distrx_c.png')
    else:
        print rxdist2
