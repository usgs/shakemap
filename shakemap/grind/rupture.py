#!/usr/bin/env python

# stdlib modules
import copy
from abc import ABC
from abc import abstractmethod
import json

# third party imports
import numpy as np
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection
from openquake.hazardlib.gsim import base

from ..utils.ecef import latlon2ecef
from ..utils.ecef import ecef2latlon
from ..utils.vector import Vector

# local imports
from shakemap.utils.exception import ShakeMapException


#-------------------------------------------------------------------------------
# CONSTANTS

# Depth tolerance in km (for determining if top and bottom edges are horizontal)
DEPTH_TOL = 0.05

# Maximum ratio of distance off of the plane (relative to edge length) for the 
# 4th point to be before being considered non-co-planar and adjusted to actually
# be on the plane?
OFFPLANE_TOLERANCE = 0.05

RAKEDICT = {'SS': 0.0, 'NM': -90.0, 'RS': 90.0, 'ALL': None}
#-------------------------------------------------------------------------------



class Rupture(ABC):
    """
    Abstract base class for ruptures.

    Note on terminology:

        - There are three Ruptuer subclasses: PointRupture, QuadRupture, and
          EdgeRupture.
        - PointRupture represents the rupture as a point source.
        - QuadRupture and EdgeRupture are two different finite source 
          representations.
        - A finite rupture is composed of segments. For QuadRupture, a segment
          is a quadrilaterial; for an EdgeRupture, a segment is a line
          connecting two points.
        - Segments are grouped with a common "group index".
        - Segments within a group must be continuous.
        - The QuadRupture class requires that each segment is a quadrilateral
          with horizontal top and obttom edges.
        - The EdgeRupture class allows for arbitrarily complex top and bottom
          edge specification. 

    """

#    @abstractmethod
#    def writeGeoJson(self, origin):
#        """
#        Write the rupture/origin info to a GeoJson file.
#
#        Args:
#            origin (Origin): Instance of ShakeMap Origin class. 
#
#        """
#        pass

    @abstractmethod
    def getLength(self):
        """
        Returns:
            float: Rupture length in km.
        """
        pass

    @abstractmethod
    def getWidth(self):
        """
        Returns:
            float: Rupture width in km.
        """
        pass

    @abstractmethod
    def getArea(self):
        """
        Returns:
            float: Rupture area in square km.
        """
        pass

    @abstractmethod
    def getStrike(self):
        """
        Return strike angle. If rupture consists of multiple quadrilaterals, the
        average strike angle, weighted by quad length, is returned.
        Note: for ruptures with quads where the strike angle changes by 180 deg
        due to reverses in dip direction are problematic and not handeled well
        by this algorithm.

        Returns:
            float: Strike angle in degrees.
        
        """
        pass

    @abstractmethod
    def getDip(self):
        pass

    @abstractmethod
    def getDepthToTop(self):
        """
        Returns:
           float: Average dip in degrees.

        """
        pass

    @classmethod
    @abstractmethod
    def fromJson(cls, d):
        """
        Args:
           d (dict): Rupture GeoJSON dictionary.

        Returns:
           Rupture subclass instance

        """
        pass


    def getRuptureContext(self, gmpelist, origin):
        """
        Returns:
            An Openquake 
        `RuptureContext <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.RuptureContext>`__.

        Args:
            gmpelist (list): List of hazardlib GMPE objects.
            origin (Origin): Instance of ShakeMap Origin class. 

        Returns:
            RuptureContext object with all known parameters filled in.

        """
        # rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'

        rx = base.RuptureContext()
        rx.mag = origin.mag
        rx.strike = self.getStrike()
        rx.dip = self.getDip()
        rx.ztor = self.getDepthToTop()
        rx.width = self.getWidth()
#            rup.strike = DEFAULT_STRIKE
#            rup.dip = DEFAULT_DIP
#            rup.ztor = DEFAULT_ZTOR
#            rup.width = DEFAULT_WIDTH

        if hasattr(origin, 'rake'):
            rx.rake = origin.rake
        elif hasattr(origin, 'mech'):
            mech = origin.mech
            rx.rake = RAKEDICT[mech]
        else:
            rx.rake = RAKEDICT['ALL']

        rx.hypo_lat = origin.lat
        rx.hypo_lon = origin.lon
        rx.hypo_depth = origin.depth

        return rx
    


class EdgeRupture(Rupture):
    """
    Rupture class that representst the rupture surface by specifying the top
    edge and the bottom edges. These edges do not need to be horizontal. The
    freedom to allow for non-horizontal edges (as compared to QuadRupture)
    comes at the cost of slower distance calculations. This is because the 
    rupture must be discretized and then the distances are compued in a brute
    force fashion based on this mesh, which can be quite large. 
    """


    def __init__(self, toplons, toplats, topdeps, botlons, botlats, botdeps,
                 group_index = None, reference = ''):
        """
        Constructor for EdgeRupture class. 

        Args:
            toplons (ndarray): Array of top edge longitudes.
            toplats (ndarray): Array of top edge latitudes. 
            topdeps (ndarray): Array of top edge depths (km).
            botlons (ndarray): Array of bot edge longitudes.
            botlats (ndarray): Array of bot edge latitudes. 
            botdeps (ndarray): Array of bot edge depths (km).
            group_index (ndarray): Optional array of group index. 
                If None, then assume only single group. 
            reference (str): Citable reference for rupture. 
 
        Returns: 
            EdgeRupture instance.

        """
        self._toplons = toplons
        self._toplats = toplats
        self._topdeps = topdeps
        self._botlons = botlons
        self._botlats = botlats
        self._botdeps = botdeps
        if group_index is not None:
            self._group_index = group_index
        else:
            self._group_index = np.zeros_like(toplons)
        self._reference = reference
        self._computeStikeDip()
        


    def getLength(self):
        """
        Compute length of rupture (km). For EdgeRupture, we compute the length
        as the length of the top edge projected to the surface.

        Returns:
            float: Rupture length in km.
        """
        lons = self._toplons
        lats = self._toplats
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        rlength = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(lons[ind], lats[ind])
                P1 = Point(lons[ind + 1], lats[ind + 1])
                dist = P0.distance(P1)
                rlength = rlength + dist
        return rlength

    def getWidth(self):
        """
        Compute average rupture width in km. For EdgeRupture, we compute this as
        the rupture area divided by teh rupture length. 

        Returns:
            float: Rupture width in km.
        """
        area = self.getArea()
        length = self.getLength()
        return area/length

    def getArea(self):
        """
        Compute the rupture area. For EdgeRupture, we compute this by grouping
        the traces into "quadrilaterals" for which the verticies may not be
        co-planar. We then break up the quadrilaterals into triangles for which
        we can compute area. 

        Returns:
            float: Rupture area in square km.
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        area = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                p0 = latlon2ecef(self._toplats[ind],
                                 self._toplons[ind],
                                 self._topdeps[ind])
                p1 = latlon2ecef(self._toplats[ind+1],
                                 self._toplons[ind+1],
                                 self._topdeps[ind+1])
                p2 = latlon2ecef(self._botlats[ind+1],
                                 self._botlons[ind+1],
                                 self._botdeps[ind+1])
                p3 = latlon2ecef(self._botlats[ind],
                                 self._botlons[ind],
                                 self._botdeps[ind])
                a = np.sqrt((p1[0] - p0[0])**2 +
                            (p1[1] - p0[1])**2 +
                            (p1[2] - p0[2])**2)
                b = np.sqrt((p2[0] - p0[0])**2 +
                            (p2[1] - p0[1])**2 +
                            (p2[2] - p0[2])**2)
                c = np.sqrt((p2[0] - p1[0])**2 +
                            (p2[1] - p1[1])**2 +
                            (p2[2] - p1[2])**2)
                s = (a + b + c)/2
                A1 = np.sqrt(s*(s - a)*(s - b)*(s - c))
                a = np.sqrt((p0[0] - p3[0])**2 +
                            (p0[1] - p3[1])**2 +
                            (p0[2] - p3[2])**2)
                b = np.sqrt((p2[0] - p3[0])**2 +
                            (p2[1] - p3[1])**2 +
                            (p2[2] - p3[2])**2)
                c = np.sqrt((p0[0] - p2[0])**2 +
                            (p0[1] - p2[1])**2 +
                            (p0[2] - p2[2])**2)
                s = (a + b + c)/2
                A2 = np.sqrt(s*(s - a)*(s - b)*(s - c))
                area = area + (A1 + A2)/1000/1000
        return area


    def getStrike(self):
        """
        Return representative strike for this rupture. Note that strike
        can vary along the rupture.

        Returns:
            float: Strike angle in degrees.
        """
        return self._strike

    def getDip(self):
        """
        Representative dip for this rupture. Note that dip
        can vary along the rupture.

        Returns:
            float: dip angle in degrees.
        """
        return self._dip

    def _computeStikeDip(self):
        """
        Loop over all triangles and get the average normal, north, and up vectors
        in ECEF. Use these to compute a representative strike and dip. 
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        norm_vec = Vector(0, 0, 0)
        north_vec = Vector(0, 0, 0)
        up_vec = Vector(0, 0, 0)
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(self._toplons[ind],
                           self._toplats[ind],
                           self._topdeps[ind])
                P1 = Point(self._toplons[ind+1],
                           self._toplats[ind+1],
                           self._topdeps[ind+1])
                P2 = Point(self._botlons[ind+1],
                           self._botlats[ind+1],
                           self._botdeps[ind+1])
                P3 = Point(self._botlons[ind],
                           self._botlats[ind],
                           self._botdeps[ind])
                P1up = Point(self._toplons[ind+1],
                             self._toplats[ind+1],
                             self._topdeps[ind+1]-1.0)
                P1N = Point(self._toplons[ind+1],
                            self._toplats[ind+1]+0.001,
                            self._topdeps[ind+1])
                P3up = Point(self._botlons[ind],
                             self._botlats[ind],
                             self._botdeps[ind]-1.0)
                P3N = Point(self._botlons[ind],
                            self._botlats[ind]+0.001,
                            self._botdeps[ind])
                p0 = Vector.fromPoint(P0)
                p1 = Vector.fromPoint(P1)
                p2 = Vector.fromPoint(P2)
                p3 = Vector.fromPoint(P3)
                p1up = Vector.fromPoint(P1up)
                p1N = Vector.fromPoint(P1N)
                p3up = Vector.fromPoint(P3up)
                p3N = Vector.fromPoint(P3N)

                # Sides
                s01 = p1 - p0
                s02 = p2 - p0
                s03 = p3 - p0
                s21 = p1 - p2
                s23 = p3 - p2
                
                # First triangle
                t1norm = (s02.cross(s01)).norm()
                a = s01.mag()
                b = s02.mag()
                c = s21.mag()
                s = (a + b + c)/2
                A1 = np.sqrt(s*(s - a)*(s - b)*(s - c))/1000

                # Second triangle
                t2norm = (s03.cross(s02)).norm()
                a = s03.mag()
                b = s23.mag()
                c = s02.mag()
                s = (a + b + c)/2
                A2 = np.sqrt(s*(s - a)*(s - b)*(s - c))/1000

                # Up and North
                p1up = (p1up - p1).norm()
                p3up = (p3up - p3).norm()
                p1N = (p1N - p1).norm()
                p3N = (p3N - p3).norm()

                # Combine
                norm_vec = norm_vec + A1*t1norm + A2*t2norm
                north_vec = north_vec + A1*p1N + A2*p3N
                up_vec = up_vec + A1*p1up + A2*p3up

        norm_vec = norm_vec.norm()
        north_vec = north_vec.norm()
        up_vec = up_vec.norm()

        # Do I need to flip the vector because it is pointing down (i.e.,
        # right-hand rule is violated)?
        flip = np.sign(up_vec.dot(norm_vec))
        norm_vec = flip*norm_vec

        # Angle between up_vec and norm_vec is dip
        self._dip = np.arcsin(up_vec.cross(norm_vec).mag())*180/np.pi

        # Normal vector projected to horizontal plane
        nvph = (norm_vec - up_vec.dot(norm_vec)*up_vec).norm()

        # Dip direction is angle between nvph and north; strike is orthogonal.
        cp = nvph.cross(north_vec)
        sign = np.sign(cp.dot(up_vec))
        dp = nvph.dot(north_vec)
        strike = np.arctan2(sign*cp.mag(), dp)*180/np.pi - 90
        if strike < -180:
            strike = strike + 360
        self._strike = strike


    def getDepthToTop(self):
        """
        Returns:
            float: Depth to top of rupture in km.
        """
        return np.min(self._topdeps)

    @classmethod
    def fromJson(cls, d):
        """
        Class method for constructing an EdgeRupture from a GeoJSON dictionary.

        Args: 
            d (dict): Rupture GeoJSON dictionary.

        Returns: 
            EdgeRupture instance.

        """
        polys = d['features'][0]['geometry']['coordinates'][0]
        n_polygons = len(polys)
        toplons = np.empty(shape=(0, 0))
        toplats = np.empty(shape=(0, 0))
        topdeps = np.empty(shape=(0, 0))
        botlons = np.empty(shape=(0, 0))
        botlats = np.empty(shape=(0, 0))
        botdeps = np.empty(shape=(0, 0))
        g_ind = 0
        group_index = []
        for i in range(n_polygons):
            p = polys[i]
            n_points = len(p)
            n_pairs = int((n_points - 1)/2)
            
            p_lons = [pt[0] for pt in p][0:-1]
            p_lats = [pt[1] for pt in p][0:-1]
            p_depths = [pt[2] for pt in p][0:-1]

            tlon = np.array(p_lons[0:n_pairs])
            blon = np.array(p_lons[(n_pairs):])[::-1]
            tlat = np.array(p_lats[0:n_pairs])
            blat = np.array(p_lats[(n_pairs):])[::-1]
            tdep = np.array(p_depths[0:n_pairs])
            bdep = np.array(p_depths[(n_pairs):])[::-1]
            
            toplons = np.append(toplons, tlon)
            toplats = np.append(toplats, tlat)
            topdeps = np.append(topdeps, tdep)
            botlons = np.append(botlons, blon)
            botlats = np.append(botlats, blat)
            botdeps = np.append(botdeps, bdep)

            group_index.extend([g_ind]*n_pairs)
            g_ind = g_ind + 1

        reference = d['features'][0]['properties']['reference']

        return cls(toplons, toplats, topdeps, botlons, botlats, botdeps,
                   group_index = group_index, reference = reference)

    def writeGeoJson(self):
        pass

    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals. Unlike QuadRupture, these
        quadrilaterals are not restricted to be coplanar or have 
        horizontal top/bottom edges.
        
        Return:
            list: List of quadrilaterals, where each quadrilateral is
                a list of OQ Points.
        """
        ugroup = np.unique(self._group_index)
        ngroup = len(ugroup)
        qlist = []
        for i in range(ngroup):
            ind = np.where(self._group_index == ugroup[i])[0]
            nq = len(ind) - 1
            for j in range(nq):
                P0 = Point(self._toplons[j],
                           self._toplats[j],
                           self._topdeps[j])
                P1 = Point(self._toplons[j+1],
                           self._toplats[j+1],
                           self._topdeps[j+1])
                P2 = Point(self._botlons[j+1],
                           self._botlats[j+1],
                           self._botdeps[j+1])
                P3 = Point(self._botlons[j],
                           self._botlats[j],
                           self._botdeps[j])
                qlist.append([P0, P1, P2, P3])

        return qlist




class QuadRupture(Rupture):
    """
    Rupture class that represents the rupture surface as a combination of 
    quadrilaterals. Each quadrilateral must have horizontal top and bottom
    edges and must be coplanar. These restrictions make the computation of
    rupture distances more efficient. The number of points in the top edges
    must match the number of points in the bottom edge. 
    """

    def __init__(self, lon, lat, depth, reference = ''):
        """
        Constructor for QuadRupture class.

        Args:
            lon (array): Sequence of rupture longitude vertices in clockwise 
                order.
            lat (array): Sequence of rupture latitude vertices in clockwise 
                order.
            depth (array): Sequence of rupture depth vertices in clockwise order.
            reference (str): String citeable reference for Rupture.

        """
        self._lon = lon
        self._lat = lat
        self._depth = depth
        self._reference = reference
        self._setQuadrilaterals()

    def getLength(self):
        """
        Compute length of rupture based on top edge in km.

        Returns:
            float: Length of rupture (km).

        """
        flength = 0
        for quad in self._quadrilaterals:
            flength = flength + get_quad_length(quad)
        return flength

    def getWidth(self):
        """
        Compute average rupture width (km) for all quadrilaterals defined for
        the rupture.

        Returns:
            float: Average width in km of all rupture quadrilaterals.
        """
        wsum = 0.0
        for quad in self._quadrilaterals:
            wsum = wsum + get_quad_width(quad)
        mwidth = (wsum / len(self._quadrilaterals)) / 1000.0
        return mwidth

    def getArea(self):
        """
        Compute area of rupture. 
        
        Returns: 
            float: Rupture area in square km. 

        """
        asum = 0.0
        for quad in self._quadrilaterals:
            w = get_quad_width(quad)
            l = get_quad_length(quad)
            asum = asum + w*l
        return asum
        

    @classmethod
    def fromJson(cls, d):
        """
        Create a QuadRupture instance from a GeoJSON dictionary.

        Args:
           d (dict): Rupture GeoJSON dictionary.

        Returns:
            QuadRupture instance.

        """

        #-----------------------------------------------------------------------
        # Assemble the vertices as the constructor needs them:
        #     For each group, where a group consists of N connected quads:
        #       1) N+1 vertices along the top edge
        #       2) N+1 vertices along the bottom edge
        #       3) First vertex repeated
        #       4) A nan (to separate groups)
        #-----------------------------------------------------------------------


        polys = d['features'][0]['geometry']['coordinates'][0]
        n_polygons = len(polys)
        lon = []
        lat = []
        dep = []
        for i in range(n_polygons):
            p = polys[i]
            p_lons = [pt[0] for pt in p][0:-1]
            p_lats = [pt[1] for pt in p][0:-1]
            p_depths = [pt[2] for pt in p][0:-1]
            lon = lon + p_lons + [np.nan]
            lat = lat + p_lats + [np.nan]
            dep = dep + p_depths + [np.nan]

            
        return cls(lon, lat, dep, d['features'][0]['properties']['reference'])

    @classmethod
    def fromTrace(cls, xp0, yp0, xp1, yp1, zp, widths, dips, strike=None,
                  reference=None):
        """
        Create a QuadRupture instance from a set of vertices that define the top
        of the rupture, and an array of widths/dips.

        Each rupture quadrilaterial is defined by specifying the latitude,
        longitude, and depth of the two vertices on the top edges, which must
        have the dame depths. The other verticies are then constructed from
        the top edges and the width and dip of the quadrilateral. 

        Args:
            xp0 (array): Array or list of longitudes (floats) of p0.
            yp0 (array): Array or list of latitudes (floats) of p0.
            xp1 (array): Array or list of longitudes (floats) of p1.
            yp1 (array): Array or list of latitudes (floats) of p1.
            zp (array): Array or list of depths for each of the top of rupture
                rectangles (km).
            widths (array): Array of widths for each of rectangle (km).
            dips (array): Array of dips for each of rectangle (degrees).
            strike (array): If None then strike is computed from verticies of
                top edge of each quadrilateral. If a scalar, then all
                quadrilaterals are constructed assuming this strike direction.
                If an array with the same length as the trace coordinates then
                it specifies the strike for each quadrilateral.
            reference (str): String explaining where the rupture definition came
                from (publication style reference, etc.).

        Returns:
            QuadRupture instance.

        """
        if len(xp0) == len(yp0) == len(xp1) == len(
                yp1) == len(zp) == len(dips) == len(widths):
            pass
        else:
            raise ShakeMapException(
                'Number of xp0,yp0,xp1,yp1,zp,widths,dips points must be equal.')
        if strike is None:
            pass
        else:
            if (len(xp0) == len(strike)) | (len(strike) == 1):
                pass
            else:
                raise ShakeMapException(
                    'Strike must be None, scalar, or same length as trace coordinates.')

        # Convert dips to radians
        dips = np.radians(dips)

        # Ensure that all input sequences are numpy arrays
        xp0 = np.array(xp0, dtype='d')
        xp1 = np.array(xp1, dtype='d')
        yp0 = np.array(yp0, dtype='d')
        yp1 = np.array(yp1, dtype='d')
        zp = np.array(zp, dtype='d')
        widths = np.array(widths, dtype='d')
        dips = np.array(dips, dtype='d')

        # Get a projection object
        west = np.min((xp0.min(), xp1.min()))
        east = np.max((xp0.max(), xp1.max()))
        south = np.min((yp0.min(), yp1.min()))
        north = np.max((yp0.max(), yp1.max()))

        # Projected coordinates are in km
        proj = get_orthographic_projection(west, east, north, south)
        xp2 = np.zeros_like(xp0)
        xp3 = np.zeros_like(xp0)
        yp2 = np.zeros_like(xp0)
        yp3 = np.zeros_like(xp0)
        zpdown = np.zeros_like(zp)
        for i in range(0, len(xp0)):
            # Project the top edge coordinates
            p0x, p0y = proj(xp0[i], yp0[i])
            p1x, p1y = proj(xp1[i], yp1[i])

            # Get the rotation angle defined by these two points
            if strike is None:
                dx = p1x - p0x
                dy = p1y - p0y
                theta = np.arctan2(dx, dy)  # theta is angle from north
            elif len(strike) == 1:
                theta = np.radians(strike[0])
            else:
                theta = np.radians(strike[i])

            R = np.array([[np.cos(theta), -np.sin(theta)],
                          [np.sin(theta), np.cos(theta)]])

            # Rotate the top edge points into a new coordinate system (vertical
            # line)
            p0 = np.array([p0x, p0y])
            p1 = np.array([p1x, p1y])
            p0p = np.dot(R, p0)
            p1p = np.dot(R, p1)

            # Get right side coordinates in project, rotated system
            dz = np.sin(dips[i]) * widths[i]
            dx = np.cos(dips[i]) * widths[i]
            p3xp = p0p[0] + dx
            p3yp = p0p[1]
            p2xp = p1p[0] + dx
            p2yp = p1p[1]

            # Get right side coordinates in un-rotated projected system
            p3p = np.array([p3xp, p3yp])
            p2p = np.array([p2xp, p2yp])
            Rback = np.array([[np.cos(-theta), -np.sin(-theta)],
                              [np.sin(-theta), np.cos(-theta)]])
            p3 = np.dot(Rback, p3p)
            p2 = np.dot(Rback, p2p)
            p3x = np.array([p3[0]])
            p3y = np.array([p3[1]])
            p2x = np.array([p2[0]])
            p2y = np.array([p2[1]])

            # project lower edge points back to lat/lon coordinates
            lon3, lat3 = proj(p3x, p3y, reverse=True)
            lon2, lat2 = proj(p2x, p2y, reverse=True)

            xp2[i] = lon2
            xp3[i] = lon3
            yp2[i] = lat2
            yp3[i] = lat3
            zpdown[i] = zp[i] + dz

        #-----------------------------------------------------------------------
        # Assemble the vertices as the constructor needs them:
        #     For each group, where a group consists of N connected quads:
        #       1) N+1 vertices along the top edge
        #       2) N+1 vertices along the bottom edge
        #       3) First vertex repeated
        #       4) A nan (to separate groups)
        #-----------------------------------------------------------------------
        nrects = len(zp)
        anan = np.ones_like(xp0) * np.nan
        lon = np.array(list(zip(xp0, xp1, xp2, xp3, anan))
                       ).reshape((nrects, 5)).flatten(order='C')
        lat = np.array(list(zip(yp0, yp1, yp2, yp3, anan))
                       ).reshape((nrects, 5)).flatten(order='C')

        # we need an array of depths, but we need to double each zp and zpdown
        # element we have
        dep = []
        for i in range(0, nrects):
            dep += [zp[i], zp[i], zpdown[i], zpdown[i],  np.nan]
        dep = np.array(dep)


        return cls(lon, lat, dep, reference)

    def writeTextFile(self, rupturefile):
        """
        Write rupture data to rupture file format as defined in ShakeMap 
        Software Guide.

        Args:
            Rupturefile (str): Filename of output data file OR file-like object.

        """
        if not hasattr(rupturefile, 'read'):
            f = open(rupturefile, 'wt')
        else:
            f = rupturefile  # just a reference to the input file-like object
        f.write('#%s\n' % self._reference)
        for quad in self.getQuadrilaterals():
            P0, P1, P2, P3 = quad
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write('%.4f %.4f %.4f\n' % (P1.latitude, P1.longitude, P1.depth))
            f.write('%.4f %.4f %.4f\n' % (P2.latitude, P2.longitude, P2.depth))
            f.write('%.4f %.4f %.4f\n' % (P3.latitude, P3.longitude, P3.depth))
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write(u'>\n')
        if not hasattr(rupturefile, 'read'):
            f.close()


    @classmethod
    def fromVertices(cls,
                     xp0, yp0, zp0, xp1, yp1, zp1,
                     xp2, yp2, zp2, xp3, yp3, zp3,
                     reference=None):
        """
        Create a QuadDrupture instance from the vector of vertices that fully
        define the quadrilaterals. The points p0, ..., p3 are labeled below for
        a trapezoid:

        ::

              p0--------p1
             /          |
            /           |
           p3-----------p2

        All of the following vector arguments must have the same length.

        Args: 
            xp0 (array): Array or list of longitudes (floats) of p0.
            yp0 (array): Array or list of latitudes (floats) of p0.
            zp0 (array): Array or list of depths (floats) of p0.
            xp1 (array): Array or list of longitudes (floats) of p1.
            yp1 (array): Array or list of latitudes (floats) of p1.
            zp1 (array): Array or list of depths (floats) of p1.
            xp2 (array): Array or list of longitudes (floats) of p2.
            yp2 (array): Array or list of latitudes (floats) of p2.
            zp2 (array): Array or list of depths (floats) of p2.
            xp3 (array): Array or list of longitudes (floats) of p3.
            yp3 (array): Array or list of latitudes (floats) of p3.
            zp3 (array): Array or list of depths (floats) of p3.
            reference (str): String explaining where the rupture definition came
                from (publication style reference, etc.)

        Returns:
            QuadRupture object, where the rupture is modeled as a series of
                trapezoids.

        """
        if len(xp0) == len(yp0) == len(zp0) == len(xp1) == len(yp1) == len(zp1) == \
           len(xp2) == len(yp2) == len(zp2) == len(xp3) == len(yp3) == len(zp3):
            pass
        else:
            raise ShakeMapException('All vectors specifying quadrilateral '\
                                    'vertices must have the same length.')

        nq = len(xp0)

        xp0 = np.array(xp0, dtype='d')
        yp0 = np.array(yp0, dtype='d')
        zp0 = np.array(zp0, dtype='d')
        xp1 = np.array(xp1, dtype='d')
        yp1 = np.array(yp1, dtype='d')
        zp1 = np.array(zp1, dtype='d')
        xp2 = np.array(xp2, dtype='d')
        yp2 = np.array(yp2, dtype='d')
        zp2 = np.array(zp2, dtype='d')
        xp3 = np.array(xp3, dtype='d')
        yp3 = np.array(yp3, dtype='d')
        zp3 = np.array(zp3, dtype='d')

        #-----------------------------------------------------------------------
        # Assemble the vertices as the constructor needs them:
        #     For each group, where a group consists of N connected quads:
        #       1) N+1 vertices along the top edge
        #       2) N+1 vertices along the bottom edge
        #       3) First vertex repeated
        #       4) A nan (to separate groups)
        #-----------------------------------------------------------------------

        anan = np.ones_like(xp0) * np.nan
        lon = np.array(list(zip(xp0, xp1, xp2, xp3, anan))
                       ).reshape((nq, 5)).flatten(order='C')
        lat = np.array(list(zip(yp0, yp1, yp2, yp3, anan))
                       ).reshape((nq, 5)).flatten(order='C')
        dep = np.array(list(zip(zp0, zp1, zp2, zp3, anan))
                       ).reshape((nq, 5)).flatten(order='C')

        return cls(lon, lat, dep, reference)


    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals.

        Returns:
            list: List of quadrilaterals where each quad is a tuple of four
                `Point <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py>`__
                objects.
        """
        return copy.deepcopy(self._quadrilaterals)

    def getStrike(self):
        """
        Return strike angle. If rupture consists of multiple quadrilaterals, the
        average strike angle, weighted by quad length, is returned.
        Note: for ruptures with quads where the strike angle changes by 180 deg
        due to reverses in dip direction are problematic and not handeled well
        by this algorithm.

        Returns:
            float: Strike angle in degrees.

        """
        nq = len(self._quadrilaterals)
        strikes = np.zeros(nq)
        lengths = np.zeros(nq)
        for i in range(nq):
            P0 = self._quadrilaterals[i][0]
            P1 = self._quadrilaterals[i][1]
            strikes[i] = P0.azimuth(P1)
            lengths[i] = get_quad_length(self._quadrilaterals[i])
        x = np.sin(np.radians(strikes))
        y = np.cos(np.radians(strikes))
        xbar = np.sum(x * lengths) / np.sum(lengths)
        ybar = np.sum(y * lengths) / np.sum(lengths)
        return np.degrees(np.arctan2(xbar, ybar))

    def getDepthToTop(self):
        """
        Determine shallowest vertex of entire rupture.

        :returns:
            Shallowest depth of all vertices (float).
        """
        mindep = 9999999
        for quad in self._quadrilaterals:
            P0, P1, P2, P3 = quad
            depths = np.array([P0.depth, P1.depth, P2.depth, P3.depth])
            if np.min(depths) < mindep:
                mindep = np.min(depths)
        return mindep

    def getDip(self):
        """
        Return average dip of all quadrilaterals in the rupture.

        Returns:
           float: Average dip in degrees.

        """
        dipsum = 0.0
        for quad in self._quadrilaterals:
            N = get_quad_normal(quad)
            V = get_vertical_vector(quad)
            dipsum = dipsum + np.degrees(np.arccos(Vector.dot(N, V)))
        dip = dipsum / len(self._quadrilaterals)
        return dip


    def getIndividualWidths(self):
        """
        Return an array of rupture widths (km), one for each quadrilateral
        defined for the rupture.

        Returns:
            Array of quad widths in km of all rupture quadrilaterals.
        """
        nquad = self.getNumQuads()
        widths = np.zeros(nquad)
        for i in range(nquad):
            q = self._quadrilaterals[i]
            widths[i] = get_quad_width(q) / 1000.0
        return widths

    def getIndividualTopLengths(self):
        """
        Return an array of rupture lengths along top edge (km),
        one for each quadrilateral defined for the rupture.

        :returns:
            Array of lengths in km of top edge of quadrilaterals.
        """
        nquad = self.getNumQuads()
        lengths = np.zeros(nquad)
        for i in range(nquad):
            P0, P1, P2, P3 = self._quadrilaterals[i]
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            lengths[i] = (p1 - p0).mag() / 1000.0
        return lengths





    @staticmethod
    def _fixStrikeDirection(quad):
        P0, P1, P2, P3 = quad
        eps = 1e-6
        p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
        p1 = Vector.fromPoint(P1)
        p2 = Vector.fromPoint(P2)
        p1p0 = p1 - p0
        p2p0 = p2 - p0
        qnv = Vector.cross(p2p0, p1p0).norm()
        tmp = p0 + qnv
        tmplat, tmplon, tmpz = ecef2latlon(tmp.x, tmp.y, tmp.z)
        if (tmpz - P0.depth) < eps: # If True then do nothing
            fixed = quad
        else:
            newP0 = copy.deepcopy(P1)
            newP1 = copy.deepcopy(P0)
            newP2 = copy.deepcopy(P3)
            newP3 = copy.deepcopy(P2)
            fixed = [newP0, newP1, newP2, newP3]
        return fixed
        

    def _setQuadrilaterals(self):
        """
        Create internal list of N quadrilaterals. Reverses quad if dip direction is
        incorrect.
        """

        # Make sure arrays are numpy arrays. 
        self._lon = np.array(self._lon)
        self._lat = np.array(self._lat)
        self._depth = np.array(self._depth)

        # Find the nans, which tells is where the separate polygons/groups are
        group_ends = np.where(np.isnan(self._lon))[0]
        n_groups = len(group_ends)

        # Check that arrays are the same length
        if len(self._lon) != len(self._lat) != len(self._depth):
            raise IndexError(
                'Length of input lon, lat, depth arrays must be equal')

        # Construct quads
        group_start = 0

        self._quadrilaterals = []
        self._group_index = []
        groupind = 0
        for i in range(n_groups):
            lonseg = self._lon[group_start:group_ends[i]]
            latseg = self._lat[group_start:group_ends[i]]
            depthseg = self._depth[group_start:group_ends[i]]
            
            # Each group can have many contiguous quadrilaterals defined in it
            # separations (nans) between segments mean that segments are not
            # contiguous.

            npoints = len(lonseg)
            nquads = int((npoints - 4) / 2) + 1
            quad_start = 0
            quad_end = -1
            for j in range(nquads):
                P0 = Point(lonseg[quad_start],
                           latseg[quad_start],
                           depthseg[quad_start])
                P1 = Point(lonseg[quad_start + 1],
                           latseg[quad_start + 1],
                           depthseg[quad_start + 1])
                P2 = Point(lonseg[quad_end - 1],
                           latseg[quad_end - 1],
                           depthseg[quad_end - 1])
                P3 = Point(lonseg[quad_end],
                           latseg[quad_end],
                           depthseg[quad_end])
                quad = [P0, P1, P2, P3]

                # Enforce plane by moving P2 -- already close because of check
                # in read_rupture_file/is_quadrupture_class/is_quad
                
                dummy, fixed_quad = is_quad(quad)

                # Reverse quad if necessary
                fixed_quad = self._fixStrikeDirection(fixed_quad)

                self._quadrilaterals.append(fixed_quad)
                
                quad_start = quad_start + 1
                quad_end = quad_end - 1

            group_start = group_ends[i] + 1
            self._group_index.extend([groupind] * nquads)
            groupind = groupind + 1

    def _getGroupIndex(self):
        """
        Return a list of segment indexes.

        Returns:
            list: Segment indexes; lenght equals the number of quadrilaterals.
        """
        return copy.deepcopy(self._group_index)

    def getLats(self):
        """
        Return a copy of the array of latitudes for the rupture verticies.

        Returns:
            array: Numpy array of latitude values.
        """
        return self._lat.copy()

    def getLons(self):
        """
        Return a copy of the array of longitudes for the rupture verticies.

        Returns:
            array: Numpy array of latitude values.
        """
        return self._lon.copy()

    def getDeps(self):
        """
        Return a copy of the array of depths for the rupture verticies.

        Returns:
            array: Numpy array of latitude values.
        """
        return self._depth.copy()

    def getNumGroups(self):
        """
        Return a count of the number of rupture groups.

        Returns:
            int:Rnumber of rupture groups. 

        """
        return len(np.unique(self._group_index))

    def getNumQuads(self):
        """
        Return a count of the number of rupture quadrilaterals.

        Returns:
            int: Number of rupture quadrilaterals.
        """
        return len(self._quadrilaterals)

    def getRuptureAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X, Y, Z (lon,lat,depth)
        coordinates. Rupture groups are separated by numpy.NaN values.

        Returns:
            tuple: 3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) 
                coordinates.
        """
        return (np.array(self._lon), np.array(self._lat), np.array(self._depth))

    def getRuptureAsMesh(self):
        """
        Return rupture segments as a OQ-Hazardlib Mesh object.

        Returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """
        rupture = Mesh(self._lon, self._lat, self._depth)
        return rupture


def read_rupture_file(file):
    """
    This is a module-level function to read in a rupture file. This allows for
    the ShakeMap 3 text file specification or the ShakeMap 4 JSON rupture format.
    The ShakeMap 3 (".txt" extension) only supports QuadRupture style rupture
    representation and so this method will always return a QuadRupture instance. 
    The ShakeMap 4 JSON format supports QuadRupture and EdgeRupture
    represenations and so this method detects the rupture class and returns the
    appropriate Rupture subclass instance.

    Args:
        file (srt): Path to rupture file.
    
    Returns:
        Rupture subclass instance. 

    """
    try:
        #-----------------------------------------------------------------------
        # First, try to read as a json file
        #-----------------------------------------------------------------------
        if isinstance(file, str):
            with open(file) as f:
                d = json.load(f)
        else:
            d = json.loads(str(file))

        rupt = json_to_rupture(d)
        
    except json.JSONDecodeError:
        #-----------------------------------------------------------------------
        # Reading as json failed, so hopefully it is a ShakeMap 3 text file
        #-----------------------------------------------------------------------
        try:
            d = text_to_json(file)
            rupt = json_to_rupture(d)
        except:
            raise Exception("Unknown rupture file format.")
    return rupt


def json_to_rupture(d):
    """
    Method returns either a QuadRupture or EdgeRupture object based on a 
    GeoJSON dictionary. 

    Args: 
        d (dict): Rupture GeoJSON dictionary.

    Returns:
        a Rupture subclass.

    """
    validate_json(d)

    # Is this a QuadRupture or an EdgeRupture?
    valid_quads = is_quadrupture_class(d)

    if valid_quads is True:
        rupt = QuadRupture.fromJson(d)
    else:
        rupt = EdgeRupture.fromJson(d)

    return rupt


def text_to_json(file):
    """
    Read in old ShakeMap 3 textfile rupture format and convert to GeoJSON. 

    Args:
        rupturefile (srt): Path to rupture file OR file-like object in GMT
            psxy format, where:

                * Rupture vertices are space separated lat, lon, depth triplets 
                  on a single line.
                * Rupture groups are separated by lines containing ">"
                * Rupture groups must be closed.
                * Verticies within a rupture group must start along the top edge
                  and move in the strike direction then move to the bottom edge
                  and move back in the opposite direction.

    Returns:
        dict: GeoJSON rupture dictionary.

    """

    #---------------------------------------------------------------------------
    # First read in the data
    #---------------------------------------------------------------------------
    x = []
    y = []
    z = []
    isFile = False
    if isinstance(file, str):
        isFile = True
        file = open(file, 'rt')
        lines = file.readlines()
    else:
        lines = file.readlines()
    reference = ''
    for line in lines:
        sline = line.strip()
        if sline.startswith('#'):
            reference += sline
            continue
        if sline.startswith('>'):
            if len(x):  # start of new line segment
                x.append(np.nan)
                y.append(np.nan)
                z.append(np.nan)
                continue
            else:  # start of file
                continue
        if not len(sline.strip()):
            continue
        parts = sline.split()
        if len(parts) < 3:
            raise ShakeMapException(
                'Rupture file %s has no depth values.' % file)
        y.append(float(parts[0]))
        x.append(float(parts[1]))
        z.append(float(parts[2]))
    if isFile:
        file.close()

    # Construct GeoJSON dictionary, note that some things like eventid and
    # metadata are available in the old file so these will be empty.
    # We could add an optional argument to include an Origin object, which
    # would be used to fill in these values.

    coords = []
    poly = []
    for lon, lat, dep in zip(x, y, z):
        if np.isnan(lon):
            coords.append(poly)
            poly = []
        else:
            poly.append([lon, lat, dep])
    if poly != []:
        coords.append(poly)

    d = {
        "type": "FeatureCollection",
        "metadata": {
            "magnitude": None,
            "eventtime": "",
            "eventid": "",
            "title": ""
        },
        "features":[
            {
                "type": "Feature",
                "properties": {
                    "rupture type": "rupture extent",
		    "reference": reference
                },
                "geometry": {
	            "type": "MultiPolygon",
	            "coordinates":[coords]
                }
            }
        ]
    }
    return d


def validate_json(d):
    """
    Check that the JSON format is acceptable. This is only for requirements that
    are common to both QuadRupture and EdgeRupture.

    Args:
        d (dict): Rupture JSON dictionary.
    """
    if d['type'] != 'FeatureCollection':
        raise Exception('JSON file is not a \"FeatureColleciton\".')

    if 'eventid' not in d['metadata'].keys():
        raise Exception('\"eventid\" not in metadata.')

    if len(d['features']) != 1:
        raise Exception('JSON file should contain excactly one feature.')

    f = d['features'][0]

    if 'reference' not in f['properties'].keys():
        raise Exception('Feature property dictionary should contain '\
                        '\"referencey\" key.')

    if f['type'] != 'Feature':
        raise Exception('Feature type should be \"Feature\".')

    geom = f['geometry']

    if geom['type'] != 'MultiPolygon':
        raise Exception('Geometry type should be \"MultiPolygon\".')

    if 'coordinates' not in geom.keys():
        raise Exception('Geometry dictionary should contain \"coordinates\" '\
                        'key.')

    polygons = geom['coordinates'][0]

    n_polygons = len(polygons)
    for i in range(n_polygons):
        p = polygons[i]
        n_points = len(p)
        if n_points % 2 == 0:
            raise Exception('Number of points in polyon must be odd.')

        if p[0] != p[-1]:
            raise Exception('First and last points in polygon must be '\
                            'identical.')

        n_pairs = int((n_points - 1)/2)
        for j in range(n_pairs):
            #-------------------------------------------------------------------
            # Points are paired and in each pair the top is first, as in:
            #
            #      _.-P1-._
            #   P0'        'P2---P3
            #   |                  \
            #   P7---P6----P5-------P4
            #
            # Pairs: P0-P7, P1-P6, P2-P5, P3-P4
            #-------------------------------------------------------------------
            top_depth = p[j][2]
            bot_depth = p[-(j+2)][2]
            if top_depth > bot_depth:
                raise Exception('Top points must be ordered before bottom points.')


def is_quadrupture_class(d):
    """
    Check if JSON file fulfills QuadRupture class criteria:
    
        - Are top and bottom edges horizontal?
        - Are the four points in each quad coplanar?

    Args:
        d (dict): Rupture JSON dictionary.

    Returns:
        bool: Can the rupture be represented in the QuadRupture class?
    """
    isQuad = True

    f = d['features'][0]
    geom = f['geometry']
    polygons = geom['coordinates'][0]
    n_polygons = len(polygons)
    for i in range(n_polygons):
        p = polygons[i]
        n_points = len(p)
        n_pairs = int((n_points - 1)/2)

        # Within each polygon, top and bottom edges must be horizontal
        depths = [pt[2] for pt in p]
        tops = np.array(depths[0:n_pairs])
        if not np.isclose(tops[0], tops, rtol = 0, atol = DEPTH_TOL).all():
            isQuad = False
        bots = np.array(depths[(n_pairs):-1])
        if not np.isclose(bots[0], bots, rtol = 0, atol = DEPTH_TOL).all():
            isQuad = False

        n_quads = n_pairs - 1
        for j in range(n_quads):
            # Four points of each quad should be co-planar within a tolerance
            quad = [Point(p[j][0], p[j][1], p[j][2]),
                    Point(p[j+1][0], p[j+1][1], p[j+1][2]),
                    Point(p[-(j+3)][0], p[-(j+3)][1], p[-(j+3)][2]),
                    Point(p[-(j+2)][0], p[-(j+2)][1], p[-(j+2)][2])]

            test = is_quad(quad)
            if test[0] is False:
                isQuad = False

    return isQuad


def is_quad(q):
    """
    Checks that an individual quad is coplanar. 

    Args: 
        q (list): A quadrilateral; list of four OQ Points.

    Returns:
        tuple: Bool for whether or not the points are planar within tolerance;
            and also the corrected quad where p2 is adjusted to be on the same
            plane as the other points.
    """
    P0, P1, P2, P3 = q

    # Convert points to ECEF
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Unit vector along top edge
    v0 = (p1 - p0).norm()

    # Distance along bottom edge
    d = (p3 - p2).mag()

    # New location for p2 by extending from p3 the same distance and
    # direction that p1 is from p0:
    new_p2 = p3 + v0*d

    # How far off of the plane is the origin p2?
    planepoints = [p0, p1, p2]
    dist = get_distance_to_plane(planepoints, p2)

    # Is it close enough?
    if dist / d > OFFPLANE_TOLERANCE:
        on_plane = False
    else:
        on_plane = True

    # Fixed quad
    fquad = [p0.toPoint(),
             p1.toPoint(),
             new_p2.toPoint(),
             p3.toPoint()]

    return (on_plane, fquad)


def get_quad_width(q):
    """
    Return width of an individual planar trapezoid, where the p0-p1 distance
    represents the long side.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Width of planar trapezoid.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    AB = p0 - p1
    AC = p0 - p3
    t1 = (AB.cross(AC).cross(AB)).norm()
    width = t1.dot(AC)

    return width

def get_quad_mesh(q, dx):
    """
    Create a mesh from a quadrilateal. 

    Args:
        q (list): A quadrilateral; list of four points. 
        dx (float):  Target dx in km; used to get nx and ny of mesh, but mesh
            snaps to edges of rupture so actual dx/dy will not actually equal this
        value in general.

    Returns:
        dict: Mesh dictionary, which includes numpy arrays:

        - llx: lower left x coordinate in ECEF coords.
        - lly: lower left y coordinate in ECEF coords.
        - llz: lower left z coordinate in ECEF coords.
        - ulx: upper left x coordinate in ECEF coords.
        - etc.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)
    # Get nx based on length of top edge, minimum allowed is 2
    toplen_km = get_quad_length(q)
    nx = int(np.max([round(toplen_km / dx, 0) + 1, 2]))

    # Get array of points along top and bottom edges
    xfac = np.linspace(0, 1, nx)
    topp = [p0 + (p1 - p0) * a for a in xfac]
    botp = [p3 + (p2 - p3) * a for a in xfac]

    # Get ny based on mean length of vectors connecting top and bottom points
    ylen_km = np.ones(nx)
    for i in range(nx):
        ylen_km[i] = (topp[i] - botp[i]).mag() / 1000
    ny = int(np.max([round(np.mean(ylen_km) / dx, 0) + 1, 2]))
    yfac = np.linspace(0, 1, ny)

    # Build mesh: dict of ny by nx arrays (x, y, z):
    mesh = {'x': np.zeros([ny, nx]), 'y': np.zeros(
        [ny, nx]), 'z': np.zeros([ny, nx])}
    for i in range(nx):
        mpts = [topp[i] + (botp[i] - topp[i]) * a for a in yfac]
        mesh['x'][:, i] = [a.x for a in mpts]
        mesh['y'][:, i] = [a.y for a in mpts]
        mesh['z'][:, i] = [a.z for a in mpts]

    # Make arrays of pixel corners
    mesh['llx'] = mesh['x'][1:, 0:-1]
    mesh['lrx'] = mesh['x'][1:, 1:]
    mesh['ulx'] = mesh['x'][0:-1, 0:-1]
    mesh['urx'] = mesh['x'][0:-1, 1:]
    mesh['lly'] = mesh['y'][1:, 0:-1]
    mesh['lry'] = mesh['y'][1:, 1:]
    mesh['uly'] = mesh['y'][0:-1, 0:-1]
    mesh['ury'] = mesh['y'][0:-1, 1:]
    mesh['llz'] = mesh['z'][1:, 0:-1]
    mesh['lrz'] = mesh['z'][1:, 1:]
    mesh['ulz'] = mesh['z'][0:-1, 0:-1]
    mesh['urz'] = mesh['z'][0:-1, 1:]
    mesh['cpx'] = np.zeros_like(mesh['llx'])
    mesh['cpy'] = np.zeros_like(mesh['llx'])
    mesh['cpz'] = np.zeros_like(mesh['llx'])

    # i and j are indices over subruptures
    ni, nj = mesh['llx'].shape
    for i in range(0, ni):
        for j in range(0, nj):
            # Rupture corner points
            pp0 = Vector(
                mesh['ulx'][i, j], mesh['uly'][i, j], mesh['ulz'][i, j])
            pp1 = Vector(
                mesh['urx'][i, j], mesh['ury'][i, j], mesh['urz'][i, j])
            pp2 = Vector(
                mesh['lrx'][i, j], mesh['lry'][i, j], mesh['lrz'][i, j])
            pp3 = Vector(
                mesh['llx'][i, j], mesh['lly'][i, j], mesh['llz'][i, j])
            # Find center of quad
            mp0 = pp0 + (pp1 - pp0) * 0.5
            mp1 = pp3 + (pp2 - pp3) * 0.5
            cp = mp0 + (mp1 - mp0) * 0.5
            mesh['cpx'][i, j] = cp.x
            mesh['cpy'][i, j] = cp.y
            mesh['cpz'][i, j] = cp.z
    return mesh


def get_local_unit_slip_vector(strike, dip, rake):
    """
    Compute the components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (deg) between rupture plane and the horizontal plane
            normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike) + \
        np.cos(rake) * np.sin(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike) + \
        np.cos(rake) * np.cos(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_DS(strike, dip, rake):
    """
    Compute the DIP SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.sin(rake) * np.cos(dip) * np.cos(strike)
    sy = np.sin(rake) * np.cos(dip) * np.sin(strike)
    sz = np.sin(rake) * np.sin(dip)
    return Vector(sx, sy, sz)


def get_local_unit_slip_vector_SS(strike, dip, rake):
    """
    Compute the STRIKE SLIP components of a unit slip vector.

    Args:
        strike (float): Clockwise angle (deg) from north of the line at the
            intersection of the rupture plane and the horizontal plane.
        dip (float): Angle (degrees) between rupture plane and the horizontal
            plane normal to the strike (0-90 using right hand rule).
        rake (float): Direction of motion of the hanging wall relative to the
            foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in 'local' N-S, E-W, U-D coordinates.

    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.cos(rake) * np.sin(strike)
    sy = np.cos(rake) * np.cos(strike)
    sz = 0.0
    return Vector(sx, sy, sz)

def reverse_quad(q):
    """
    Reverse the verticies of a quad in the sense that the strike direction
    is flipped. 

    Args:
        q (list): A quadrilateral; list of four points.

    Returns: 
        list: Reversed quadrilateral.

    """
    return [q[1], q[0], q[3], q[2]]

def get_quad_slip(q, rake):
    """
    Compute the unit slip vector in ECEF space for a quad and rake angle.

    Args:
        q (list): A quadrilateral; list of four points.
        rake (float): Direction of motion of the hanging wall relative to
        the foot wall, as measured by the angle (deg) from the strike vector.

    Returns:
        Vector: Unit slip vector in ECEF space.

    """
    P0, P1, P2 = q[0:3]
    strike = P0.azimuth(P1)
    dip = get_quad_dip(q)
    s1_local = get_local_unit_slip_vector(strike, dip, rake)
    s0_local = Vector(0, 0, 0)
    qlats = [a.latitude for a in q]
    qlons = [a.longitude for a in q]
    proj = get_orthographic_projection(
        np.min(qlons), np.max(qlons), np.min(qlats), np.max(qlats))
    s1_ll = proj(np.array([s1_local.x]), np.array([s1_local.y]), reverse=True)
    s0_ll = proj(np.array([s0_local.x]), np.array([s0_local.y]), reverse=True)
    s1_ecef = Vector.fromTuple(latlon2ecef(s1_ll[1], s1_ll[0], s1_local.z))
    s0_ecef = Vector.fromTuple(latlon2ecef(s0_ll[1], s0_ll[0], s0_local.z))
    slp_ecef = (s1_ecef - s0_ecef).norm()
    return slp_ecef


def get_quad_length(q):
    """
    Length of top eduge of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Length of quadrilateral in km.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    qlength = (p1 - p0).mag() / 1000
    return qlength


def get_quad_dip(q):
    """
    Dip of a quadrilateral.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        float: Dip in degrees.

    """
    N = get_quad_normal(q)
    V = get_vertical_vector(q)
    dip = np.degrees(np.arccos(Vector.dot(N, V)))
    return dip


def get_quad_normal(q):
    """
    Compute the unit normal vector for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized normal vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p3 = Vector.fromPoint(P3)
    v1 = p1 - p0
    v2 = p3 - p0
    vn = Vector.cross(v2, v1).norm()
    return vn


def get_quad_strike_vector(q):
    """
    Compute the unit vector pointing in the direction of strike for a
    quadrilateral in ECEF coordinates. Top edge assumed to be horizontal.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing in strike direction in ECEF coords.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    v1 = (p1 - p0).norm()
    return v1


def get_quad_down_dip_vector(q):
    """
    Compute the unit vector pointing down dip for a quadrilateral in
    ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: The unit vector pointing down dip in ECEF coords.

    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    p0p1 = p1 - p0
    qnv = get_quad_normal(q)
    ddv = Vector.cross(p0p1, qnv).norm()
    return ddv


def get_vertical_vector(q):
    """
    Compute the vertical unit vector for a quadrilateral
    in ECEF coordinates.

    Args:
        q (list): A quadrilateral; list of four points.

    Returns:
        Vector: Normalized vertical vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    P0_up = copy.deepcopy(P0)
    P0_up.depth = P0_up.depth - 1.0
    p0 = Vector.fromPoint(P0)   # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P0_up)
    v1 = (p1 - p0).norm()
    return v1


def get_distance_to_plane(planepoints, otherpoint):
    """
    Calculate a point's distance to a plane.  Used to figure out if a
    quadrilateral points are all co-planar.

    Args:
        planepoints (list): List of three points (from Vector class) defining a
            plane.
        otherpoint (Vector): 4th Vector to compare to points defining the plane.
        
    Returns:
        float: Distance (in meters) from otherpoint to plane.

    """
    # from
    # https://en.wikipedia.org/wiki/Plane_(geometry)#Describing_a_plane_through_three_points
    p0, p1, p2 = planepoints
    x1, y1, z1 = p0.getArray()
    x2, y2, z2 = p1.getArray()
    x3, y3, z3 = p2.getArray()
    D = np.linalg.det(np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]))
    if D != 0:
        d = -1
        at = np.linalg.det(np.array([[1, y1, z1], [1, y2, z2], [1, y3, z3]]))
        bt = np.linalg.det(np.array([[x1, 1, z1], [x2, 1, z2], [x3, 1, z3]]))
        ct = np.linalg.det(np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
        a = (-d / D) * at
        b = (-d / D) * bt
        c = (-d / D) * ct

        numer = np.abs(a * otherpoint.x +
                       b * otherpoint.y +
                       c * otherpoint.z + d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer / denom
    else:
        dist = 0
    return dist



