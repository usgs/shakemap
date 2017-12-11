#!/usr/bin/env python

# stdlib modules
import copy

# third party imports
import numpy as np
from openquake.hazardlib.geo.mesh import Mesh
from openquake.hazardlib.geo.point import Point
from openquake.hazardlib.geo.utils import get_orthographic_projection

from impactutils.vectorutils.ecef import latlon2ecef
from impactutils.vectorutils.ecef import ecef2latlon
from impactutils.vectorutils.vector import Vector
from impactutils.time.ancient_time import HistoricTime
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture.base import Rupture
from shakelib.rupture import utils
from shakelib.rupture import gc2


class QuadRupture(Rupture):
    """
    Rupture class that represents the rupture surface as a combination of
    quadrilaterals. Each quadrilateral must have horizontal top and bottom
    edges and must be coplanar. These restrictions make the computation of
    rupture distances more efficient. The number of points in the top edges
    must match the number of points in the bottom edge.
    """

    def __init__(self, d, origin):
        """
        Create a QuadRupture instance from a GeoJSON dictionary and an Origin.

        Args:
           d (dict): Rupture GeoJSON dictionary.
           origin (Origin): Reference to a ShakeMap Origin object.

        Returns:
            QuadRupture instance.

        """

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

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v

        self._geojson = d
        self._lon = lon
        self._lat = lat
        self._depth = dep
        self._origin = origin
        self._reference = d['metadata']['reference']
        self._setQuadrilaterals()

    def getDepthAtPoint(self, lat, lon):
        SMALL_DISTANCE = 2e-03  # 2 meters
        depth = np.nan

        tmp = self.computeRjb(np.array([lon]), np.array([lat]), np.array([0]))
        if tmp > SMALL_DISTANCE:
            return depth

        i = 0
        imin = -1
        dmin = 9999999999999999
        for quad in self.getQuadrilaterals():
            pX = Vector.fromPoint(Point(lon, lat, 0))
            points = np.reshape(np.array([pX.x, pX.y, pX.z]), (1, 3))
            rjb = utils._quad_distance(quad, points, horizontal=True)
            if rjb[0][0] < dmin:
                dmin = rjb[0][0]
                imin = i
            i += 1

        quad = self._quadrilaterals[imin]
        P0, P1, P2, P3 = quad
        # project the quad and the point in question to orthographic defined by
        # quad
        xmin = np.min([P0.x, P1.x, P2.x, P3.x])
        xmax = np.max([P0.x, P1.x, P2.x, P3.x])
        ymin = np.min([P0.y, P1.y, P2.y, P3.y])
        ymax = np.max([P0.y, P1.y, P2.y, P3.y])
        proj = get_orthographic_projection(xmin, xmax, ymax, ymin)

        # project each vertex of quad (at 0 depth)
        s0x, s0y = proj(P0.x, P0.y)
        s1x, s1y = proj(P1.x, P1.y)
        s2x, s2y = proj(P2.x, P2.y)
        s3x, s3y = proj(P3.x, P3.y)
        sxx, sxy = proj(lon, lat)

        # turn these to vectors
        s0 = Vector(s0x, s0y, 0)
        s1 = Vector(s1x, s1y, 0)
        s3 = Vector(s3x, s3y, 0)
        sx = Vector(sxx, sxy, 0)

        # Compute vector from s0 to s1
        s0s1 = s1 - s0
        # Compute the vector from s0 to s3
        s0s3 = s3 - s0
        # Compute the vector from s0 to sx
        s0sx = sx - s0

        # cross products
        s0normal = s0s3.cross(s0s1)
        dd = s0s1.cross(s0normal)
        # normalize dd (down dip direction)
        ddn = dd.norm()
        # dot product
        sxdd = ddn.dot(s0sx)

        # get width of quad
        w = utils.get_quad_width(quad)

        # Get weights for top and bottom edge depths
        N = utils.get_quad_normal(quad)
        V = utils.get_vertical_vector(quad)
        dip = np.degrees(np.arccos(Vector.dot(N, V)))
        ws = (w * np.cos(np.radians(dip)))
        wtt = (ws - sxdd) / ws
        wtb = sxdd / ws

        # Compute the depth of of the plane at Px:
        depth = wtt * P0.z + wtb * P3.z * 1000

        return depth

    def getLength(self):
        """
        Compute length of rupture based on top edge in km.

        Returns:
            float: Length of rupture (km).

        """
        flength = 0
        for quad in self._quadrilaterals:
            flength = flength + utils.get_quad_length(quad)
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
            wsum = wsum + utils.get_quad_width(quad)
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
            width = utils.get_quad_width(quad)
            length = utils.get_quad_length(quad)
            asum = asum + width * length
        return asum

    @classmethod
    def fromTrace(cls, xp0, yp0, xp1, yp1, zp, widths, dips, origin,
                  strike=None, group_index=None, reference=""):
        """
        Create a QuadRupture instance from a set of vertices that define the
        top of the rupture, and an array of widths/dips.

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
            origin (Origin): Reference to a ShakeMap origin object.
            strike (array): If None then strike is computed from verticies of
                top edge of each quadrilateral. If a scalar, then all
                quadrilaterals are constructed assuming this strike direction.
                If an array with the same length as the trace coordinates then
                it specifies the strike for each quadrilateral.
            group_index (list): List of integers to indicate group index. If
                None then each quadrilateral is assumed to be in a different
                group since there is no guarantee that any of them are
                continuous.
            reference (str): String explaining where the rupture definition
                came from (publication style reference, etc.).

        Returns:
            QuadRupture instance.

        """
        if len(xp0) == len(yp0) == len(xp1) == len(
                yp1) == len(zp) == len(dips) == len(widths):
            pass
        else:
            raise ShakeLibException(
                'Number of xp0,yp0,xp1,yp1,zp,widths,dips points must be '
                'equal.')
        if strike is None:
            pass
        else:
            if (len(xp0) == len(strike)) | (len(strike) == 1):
                pass
            else:
                raise ShakeLibException(
                    'Strike must be None, scalar, or same length as '
                    'trace coordinates.')

        if group_index is None:
            group_index = np.array(range(len(xp0)))

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

        # ---------------------------------------------------------------------
        # Create GeoJSON object
        # ---------------------------------------------------------------------

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate(
                    [xp0[ind[0]].reshape((1,)),
                     xp1[ind], xp2[ind][::-1],
                     xp3[ind][::-1][-1].reshape((1,)),
                     xp0[ind[0]].reshape((1,))
                     ])
            lats = np.concatenate(
                    [yp0[ind[0]].reshape((1,)),
                     yp1[ind],
                     yp2[ind][::-1],
                     yp3[ind][::-1][-1].reshape((1,)),
                     yp0[ind[0]].reshape((1,))
                     ])
            deps = np.concatenate(
                    [zp[ind[0]].reshape((1,)),
                     zp[ind],
                     zpdown[ind][::-1],
                     zpdown[ind][::-1][-1].reshape((1,)),
                     zp[ind[0]].reshape((1,))])

            poly = []
            for lon, lat, dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)

        d = {"type": "FeatureCollection",
             "metadata": {
                 "reference": reference
             },
             "features": [{
                 "type": "Feature",
                 "properties": {
                     "rupture type": "rupture extent"
                 },
                 "geometry": {
                     "type": "MultiPolygon",
                     "coordinates": [coords]
                 }
             }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v

        return cls(d, origin)

    def writeTextFile(self, rupturefile):
        """
        Write rupture data to rupture file format as defined in ShakeMap
        Software Guide.

        Note that this currently treats each quadrilateral as a separate
        polygon. This needs to be udpated.

        Args:
            rupturefile (str): Filename of output data file OR file-like
                object.

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
                     origin,
                     group_index=None,
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
            origin (Origin): Reference to a ShakeMap Origin object.
            group_index (list): List of integers to indicate group index. If
                None then each quadrilateral is assumed to be in a different
                group since there is no guarantee that any of them are
                continuous.
            reference (str): String explaining where the rupture definition
                came from (publication style reference, etc.)

        Returns:
            QuadRupture object, where the rupture is modeled as a series of
                trapezoids.

        """
        if len(xp0) == len(yp0) == len(zp0) == len(xp1) == len(yp1) == \
           len(zp1) == len(xp2) == len(yp2) == len(zp2) == len(xp3) == \
           len(yp3) == len(zp3):
            pass
        else:
            raise ShakeLibException('All vectors specifying quadrilateral '
                                    'vertices must have the same length.')

        nq = len(xp0)
        if group_index is not None:
            if len(group_index) != nq:
                raise Exception(
                    "group_index must have same length as vertices.")
        else:
            group_index = np.array(range(nq))

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

        # ---------------------------------------------------------------------
        # Create GeoJSON object
        # ---------------------------------------------------------------------

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate(
                    [xp0[ind[0]].reshape((1,)),
                     xp1[ind],
                     xp2[ind][::-1],
                     xp3[ind][::-1][-1].reshape((1,)),
                     xp0[ind[0]].reshape((1,))
                     ])
            lats = np.concatenate(
                    [yp0[ind[0]].reshape((1,)),
                     yp1[ind],
                     yp2[ind][::-1],
                     yp3[ind][::-1][-1].reshape((1,)),
                     yp0[ind[0]].reshape((1,))
                     ])
            deps = np.concatenate(
                    [zp0[ind[0]].reshape((1,)),
                     zp1[ind],
                     zp2[ind][::-1],
                     zp3[ind][::-1][-1].reshape((1,)),
                     zp0[ind[0]].reshape((1,))
                     ])

            poly = []
            for lon, lat, dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)

        d = {"type": "FeatureCollection",
             "metadata": {
                 "reference": reference
             },
             "features": [{
                 "type": "Feature",
                 "properties": {
                     "rupture type": "rupture extent"
                 },
                 "geometry": {
                     "type": "MultiPolygon",
                     "coordinates": [coords]
                 }
             }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v
        if hasattr(origin, 'id'):
            d['metadata']['eventid'] = origin.id

        return cls(d, origin)

    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals.

        Returns:
            list: List of quadrilaterals where each quad is a tuple of four
                `Point <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py>`__
                objects.
        """  # noqa
        return copy.deepcopy(self._quadrilaterals)

    def getStrike(self):
        """
        Return strike angle. If rupture consists of multiple quadrilaterals,
        the average strike angle, weighted by quad length, is returned.
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
            lengths[i] = utils.get_quad_length(self._quadrilaterals[i])
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
            N = utils.get_quad_normal(quad)
            V = utils.get_vertical_vector(quad)
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
            widths[i] = utils.get_quad_width(q) / 1000.0
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
        if (tmpz - P0.depth) < eps:  # If True then do nothing
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
        Create internal list of N quadrilaterals. Reverses quad if dip
        direction is incorrect.
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

                dummy, fixed_quad = utils.is_quad(quad)

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
        Return a list of segment group indexes.

        Returns:
            list: Segment group indexes; length equals the number of
                quadrilaterals.
        """
        return copy.deepcopy(self._group_index)

    @property
    def lats(self):
        """
        Return an array of latitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans.
        """
        lats = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_lats = []
            bot_lats = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].latitude]
                    bot0 = [quads[q_ind[j]][3].latitude]
                    top_lats = top_lats + top0
                    bot_lats = bot_lats + bot0
                top_lats = top_lats + [quads[q_ind[j]][1].latitude]
                bot_lats = bot_lats + [quads[q_ind[j]][2].latitude]
            lats = lats + top_lats + bot_lats[::-1] + top0 + [np.nan]

        return np.array(lats)

    @property
    def lons(self):
        """
        Return an array of longitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop longitude values; disconnected
                segments are separated by nans.
        """
        lons = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_lons = []
            bot_lons = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].longitude]
                    bot0 = [quads[q_ind[j]][3].longitude]
                    top_lons = top_lons + top0
                    bot_lons = bot_lons + bot0
                top_lons = top_lons + [quads[q_ind[j]][1].longitude]
                bot_lons = bot_lons + [quads[q_ind[j]][2].longitude]
            lons = lons + top_lons + bot_lons[::-1] + top0 + [np.nan]
        return np.array(lons)

    @property
    def depths(self):
        """
        Return an array of depths for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop depths; disconnected
                segments are separated by nans.
        """
        deps = []
        quads = self.getQuadrilaterals()
        groups = self._getGroupIndex()
        u_groups = np.unique(groups)
        ng = len(u_groups)
        for i in range(ng):
            q_ind = np.where(groups == u_groups[i])[0]
            nq = len(q_ind)
            top_deps = []
            bot_deps = []
            for j in range(nq):
                if j == 0:
                    top0 = [quads[q_ind[j]][0].depth]
                    bot0 = [quads[q_ind[j]][3].depth]
                    top_deps = top_deps + top0
                    bot_deps = bot_deps + bot0
                top_deps = top_deps + [quads[q_ind[j]][1].depth]
                bot_deps = bot_deps + [quads[q_ind[j]][2].depth]
            deps = deps + top_deps + bot_deps[::-1] + top0 + [np.nan]

        return np.array(deps)

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
        return (np.array(self._lon),
                np.array(self._lat),
                np.array(self._depth))

    def getRuptureAsMesh(self):
        """
        Return rupture segments as a OQ-Hazardlib Mesh object.

        Returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """  # noqa
        rupture = Mesh(self._lon, self._lat, self._depth)
        return rupture

    def computeRjb(self, lon, lat, depth):
        """
        Method for computing Joyner-Boore distance.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Joyner-Boore distance (km).

        """

        # ---------------------------------------------------------------------
        # Sort out sites
        # ---------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        minrjb = np.ones(newshape, dtype=lon.dtype) * 1e16
        quads = self.getQuadrilaterals()

        for i in range(len(quads)):
            P0, P1, P2, P3 = quads[i]
            S0 = copy.deepcopy(P0)
            S1 = copy.deepcopy(P1)
            S2 = copy.deepcopy(P2)
            S3 = copy.deepcopy(P3)
            S0.depth = 0.0
            S1.depth = 0.0
            S2.depth = 0.0
            S3.depth = 0.0
            squad = [S0, S1, S2, S3]
            rjbdist = utils._quad_distance(squad, sites_ecef, horizontal=True)
            minrjb = np.minimum(minrjb, rjbdist)

        minrjb = minrjb.reshape(oldshape)
        return minrjb

    def computeRrup(self, lon, lat, depth):
        """
        Method for computing rupture distance.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Rupture distance (km).

        """

        # ---------------------------------------------------------------------
        # Sort out sites
        # ---------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        minrrup = np.ones(newshape, dtype=lon.dtype) * 1e16
        quads = self.getQuadrilaterals()

        for i in range(len(quads)):
            rrupdist = utils._quad_distance(quads[i], sites_ecef)
            minrrup = np.minimum(minrrup, rrupdist)

        minrrup = minrrup.reshape(oldshape)
        return minrrup

    def computeGC2(self, lon, lat, depth):
        """
        Method for computing version 2 of the Generalized Coordinate system
        (GC2) by Spudich and Chiou OFR 2015-1028.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            dict: Dictionary with keys for each of the GC2-related distances,
                which include 'rx', 'ry', 'ry0', 'U', 'T'.
        """
        # This just hands off to the module-level method
        dict = gc2._computeGC2(self, lon, lat, depth)
        return dict
