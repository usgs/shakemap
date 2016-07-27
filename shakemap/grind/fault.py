#!/usr/bin/env python

# stdlib modules
import io
import sys
import copy

# third party imports
import numpy as np
from openquake.hazardlib.geo import mesh
from openquake.hazardlib.geo import point, utils
from openquake.hazardlib.geo.utils import get_orthographic_projection
from .ecef import latlon2ecef
from .ecef import ecef2latlon
from .vector import Vector
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# local imports
from shakemap.utils.exception import ShakeMapException

# CONSTANTS
# what is the maximum ratio of distance out of the plane defined by 3 points a
# 4th point can be before being considered non-co-planar?
OFFPLANE_TOLERANCE = 0.05


class Fault(object):
    """
    Class to handle fault files of various types and output fault data in 
    various ways.
    """

    def __init__(self, lon, lat, depth, reference):
        """
        Constructor for Fault.
        :param lon:
            sequence of fault longitude vertices in clockwise order.
        :param lat:
            sequence of fault latitude vertices in clockwise order.
        :param depth:
            sequence of fault depth vertices in clockwise order.
        :param reference:
            string citeable reference for Fault.
        """
        self._lon = lon
        self._lat = lat
        self._depth = depth
        self._reference = reference
        self._validate()
        self.setQuadrilaterals()

    @classmethod
    def fromTrace(cls, xp0, yp0, xp1, yp1, zp, widths, dips, strike=None,
                  reference=None):
        """
        Create a fault object from a set of vertices that define the top of the
        fault, and an array of widths/dips.

        These top of rupture points are defined by specifying the x and y
        coordinates of each of the two vertices, and then specifying an array of
        depths,widths, and dips for each rectangle.
        :param xp0:
          Array of longitude coordinates for the first (top of rupture) vertex of
          each rectangle (decimal degrees).
        :param yp0:
          Array of latitude coordinates for the first (top of rupture) vertex of
          each rectangle (decimal degrees).
        :param xp1:
          Array of longitude coordinates for the second (top of rupture) vertex of
          each rectangle (decimal degrees).
        :param yp1:
          Array of latitude coordinates for the second (top of rupture) vertex of
          each rectangle (decimal degrees).
        :param zp:
          Array of depths for each of the top of rupture rectangles (km).
        :param widths:
          Array of widths for each of rectangle (km).
        :param dips:
          Array of dips for each of rectangle (degrees).
        :param strike:
          If None then strike is computed from verticies of top edge of each
          quadrilateral. If a scalar, then all quadrilaterals are constructed
          assuming this strike direction. If a vector with the same length as
          the trace coordinates then it specifies the strike for each quadrilateral.
        :param reference:
          String explaining where the fault definition came from (publication style
          reference, etc.)
        :returns:
          Fault object, where the fault is modeled as a series of rectangles.
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

        # convert dips to radians
        dips = np.radians(dips)

        # ensure that all input sequences are numpy arrays
        xp0 = np.array(xp0, dtype='d')
        xp1 = np.array(xp1, dtype='d')
        yp0 = np.array(yp0, dtype='d')
        yp1 = np.array(yp1, dtype='d')
        zp = np.array(zp, dtype='d')
        widths = np.array(widths, dtype='d')
        dips = np.array(dips, dtype='d')

        # get a projection object
        west = np.min((xp0.min(), xp1.min()))
        east = np.max((xp0.max(), xp1.max()))
        south = np.min((yp0.min(), yp1.min()))
        north = np.max((yp0.max(), yp1.max()))

        # projected coordinates are in km
        proj = get_orthographic_projection(west, east, north, south)
        surfaces = []
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
                theta = np.radians(strike)
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

            # Get right side coordinates in project,rotated system
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

        # assemble the vertices as the Fault constructor needs them...
        # which is: for each rectangle, there should be the four corners, the
        # first corner repeated, and then a nan.
        nrects = len(zp)
        anan = np.ones_like(xp0) * np.nan
        lon = np.array(list(zip(xp0, xp1, xp2, xp3, xp0, anan))
                       ).reshape((nrects, 6)).flatten(order='C')
        lat = np.array(list(zip(yp0, yp1, yp2, yp3, yp0, anan))
                       ).reshape((nrects, 6)).flatten(order='C')

        # we need an array of depths, but we need to double each zp and zpdown
        # element we have
        dep = []
        for i in range(0, nrects):
            dep += [zp[i], zp[i], zpdown[i], zpdown[i], zp[i], np.nan]
        dep = np.array(dep)

        # take the nans off the end of each array
        lon = lon[0:-1]
        lat = lat[0:-1]
        dep = dep[0:-1]

        return cls(lon, lat, dep, reference)

    def writeFaultFile(self, faultfile):
        """
        Write fault data to fault file format as defined in ShakeMap Software Guide.
        :param faultfile:
          Filename of output data file OR file-like object.
        """
        if not hasattr(faultfile, 'read'):
            f = open(faultfile, 'wt')
        else:
            f = faultfile  # just a reference to the input file-like object
        f.write('#%s\n' % self._reference)
        for quad in self.getQuadrilaterals():
            P0, P1, P2, P3 = quad
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write('%.4f %.4f %.4f\n' % (P1.latitude, P1.longitude, P1.depth))
            f.write('%.4f %.4f %.4f\n' % (P2.latitude, P2.longitude, P2.depth))
            f.write('%.4f %.4f %.4f\n' % (P3.latitude, P3.longitude, P3.depth))
            f.write('%.4f %.4f %.4f\n' % (P0.latitude, P0.longitude, P0.depth))
            f.write(u'>\n')
        if not hasattr(faultfile, 'read'):
            f.close()

    @classmethod
    def readFaultFile(cls, faultfile):
        """
        Read fault file format as defined in ShakeMap Software Guide.
        :param faultfile:
            Path to fault file OR file-like object in GMT psxy format, where
             * Fault vertices are space separated lat,lon,depth triplets on a single line.
             * Fault segments are separated by lines containing ">"
             * Fault segments must be closed.
             * Fault segments must be all clockwise or all counter-clockwise.
        :returns:
           Fault object.
        :raises ShakeMapException:
            when any of above conditions are not met.
        """
        x = []
        y = []
        z = []
        isFile = False
        if isinstance(faultfile, str):
            isFile = True
            faultfile = open(faultfile, 'rt')
            faultlines = faultfile.readlines()
        else:
            faultlines = faultfile.readlines()
        reference = ''
        for line in faultlines:
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
                    'Finite fault file %s has no depth values.' %
                    faultfile)
            y.append(float(parts[0]))
            x.append(float(parts[1]))
            z.append(float(parts[2]))
        if isFile:
            faultfile.close()
        if np.isnan(x[-1]):
            x = x[0:-1]
            y = y[0:-1]
            z = z[0:-1]

        return cls(x, y, z, reference)

    @classmethod
    def fromVertices(cls,
                     xp0, yp0, zp0, xp1, yp1, zp1,
                     xp2, yp2, zp2, xp3, yp3, zp3,
                     reference=None):
        """
        Create a fault object from the vector of vertices that fully define the
        quadrilaterals. The points p0, ..., p3 are labeled below for a trapezoid:

        p0--------p1
        |          \
        |           \
        p3----------p2

        All of the following vector arguments must have the same length.
        :param xp0:
          Vector of longitudes of p0.
        :param yp0:
          Vector of latitudes of p0.
        :param zp0:
          Vector of depths of p0 (positive down).
        :param xp1:
          Vector of longitudes of p1.
        :param yp1:
          Vector of latitudes of p1.
        :param zp1:
          Vector of depths of p1 (positive down).
        :param xp2:
          Vector of longitudes of p2.
        :param yp2:
          Vector of latitudes of p2.
        :param zp2:
          Vector of depths of p2 (positive down).
        :param xp3:
          Vector of longitudes of p3.
        :param yp3:
          Vector of latitudes of p3.
        :param zp3:
          Vector of depths of p3 (positive down).
        :param reference:
          String explaining where the fault definition came from (publication style
          reference, etc.)
        :returns:
          Fault object, where the fault is modeled as a series of trapezoids.
        """
        if len(xp0) == len(yp0) == len(zp0) == len(xp1) == len(yp1) == len(zp1) == \
           len(xp2) == len(yp2) == len(zp2) == len(xp3) == len(yp3) == len(zp3):
            pass
        else:
            raise ShakeMapException('All vectors specifying quadrilateral '
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

        # assemble the vertices as the Fault constructor needs them...
        # which is: for each rectangle, there should be the four corners, the
        # first corner repeated, and then a nan.
        anan = np.ones_like(xp0) * np.nan
        lon = np.array(list(zip(xp0, xp1, xp2, xp3, xp0, anan))
                       ).reshape((nq, 6)).flatten(order='C')
        lat = np.array(list(zip(yp0, yp1, yp2, yp3, yp0, anan))
                       ).reshape((nq, 6)).flatten(order='C')
        dep = np.array(list(zip(zp0, zp1, zp2, zp3, zp0, anan))
                       ).reshape((nq, 6)).flatten(order='C')

        return cls(lon, lat, dep, reference)

    def multiplyFaultLength(self, factor):
        # What does this do and why do we want to do it????
        for i in range(0, len(self._quadrilaterals)):
            quad = self._quadrilaterals[i]
            P0, P1, P2, P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            dtop = p1 - p0
            toplen = dtop.mag()
            topnorm = dtop.norm()
            dbottom = p2 - p3
            bottomlen = dbottom.mag()
            bottomnorm = dbottom.norm()
            newtoplen = toplen * factor
            newbottomlen = bottomlen * factor
            newp1 = p0 + topnorm * newtoplen
            newp2 = p3 + bottomnorm * newbottomlen
            newP1 = newp1.toPoint()
            newP2 = newp2.toPoint()
            newP1.depth = P0.depth
            newP2.depth = P3.depth
            self._quadrilaterals[i][1] = newP1
            self._quadrilaterals[i][2] = newP2

    def getFaultLength(self):
        """
        Compute lenght of fault based on top edge.
        :returns:
            Length of fault.
        """
        flength = 0
        for quad in self._quadrilaterals:
            flength = flength + get_quad_length(quad)
        return flength

    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals.
        :returns:
            each quad is a tuple of four Point objects
            (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py)
        """
        return copy.deepcopy(self._quadrilaterals)

    def getQuadWidth(self, p0, p1, p3):
        """
        Return width of an individual planar trapezoid, where the p0-p1 distance
        represents the long side.
        :param p0:
            ECEF x,y,z point representing the first vertex of a quadrilateral.
        :param p1:
            ECEF x,y,z point representing the second vertex of a quadrilateral.
        :param p3:
            ECEF x,y,z point representing the fourth vertex of a quadrilateral.
        :returns:
           width of planar trapezoid.
        """
        AB = p0 - p1
        AC = p0 - p3
        t1 = (AB.cross(AC).cross(AB)).norm()
        width = t1.dot(AC)

        return width

    def getStrike(self):
        """
        Return strike angle. If fault consists of multiple quadrilaterals, the
        average strike angle, weighted by quad length, is returned.
        Note: for faults with quads where the strike angle changes by 180 deg
        due to reverses in dip direction are problematic and not handeled well
        by this algorithm.
        :returns:
          float strike angle
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

    def getTopOfRupture(self):
        """
        Determine shallowest vertex of entire fault.
        :returns:
            float shallowest depth of all vertices.
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
        Return average dip of all quadrilaterals in the fault.
        :returns:
           Average dip in degrees.
        """
        dipsum = 0.0
        for quad in self._quadrilaterals:
            N = get_quad_normal(quad)
            V = get_vertical_vector(quad)
            dipsum = dipsum + np.degrees(np.arccos(Vector.dot(N, V)))
        dip = dipsum / len(self._quadrilaterals)
        return dip

    def getWidth(self):
        """
        Return the average fault width (km) for all quadrilaterals defined for the fault.
        :returns:
            Average width in km of all fault quadrilaterals.
        """
        wsum = 0.0
        for quad in self._quadrilaterals:
            P0, P1, P2, P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p3 = Vector.fromPoint(P3)
            wsum += self.getQuadWidth(p0, p1, p3)
        mwidth = (wsum / len(self._quadrilaterals)) / 1000.0
        return mwidth

    def getIndividualWidths(self):
        """
        Return an array of fault widths (km), one for each quadrilateral defined for the fault.
        :returns:
            Array of quad widths in km of all fault quadrilaterals.
        """
        nquad = self.getNumQuads()
        widths = np.zeros(nquad)
        for i in range(nquad):
            P0, P1, P2, P3 = self._quadrilaterals[i]
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p3 = Vector.fromPoint(P3)
            widths[i] = self.getQuadWidth(p0, p1, p3) / 1000.0
        return widths

    def getIndividualTopLengths(self):
        """
        Return an array of fault lengths along top edge (km),
        one for each quadrilateral defined for the fault.
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

    def _getTrapMeanLength(self, p0, p1, p2, p3):
        """
        Return the sqrt of the area of a quadrilateral (used for QA of fault plane).
        :param p0:
            ECEF x,y,z point representing the first vertex of a quadrilateral.
        :param p1:
            ECEF x,y,z point representing the second vertex of a quadrilateral.
        :param p2:
            ECEF x,y,z point representing the third vertex of a quadrilateral.
        :param p3:
            ECEF x,y,z point representing the fourth vertex of a quadrilateral.
        :returns:
            square root of trapezoid area.
        """
        # area of a trapezoid: A = (a+b)/2 * h
        # (https://en.wikipedia.org/wiki/Trapezoid)
        h = self.getQuadWidth(p0, p1, p3)
        a = (p1 - p0).mag()
        b = (p2 - p3).mag()
        A = ((a + b) / 2.0) * h
        length = np.sqrt(A)
        return length

    def getDistanceToPlane(self, planepoints, otherpoint):
        """
        Calculate a point's distance to a plane.  Used to figure out if a quadrilateral points are all co-planar.
        :param planepoints:
            List of three points (Vector objects) defining a plane.
        :param otherpoint:
            4th Vector to compare to points defining the plane
        :returns:
            Distance (in meters) from otherpoint to plane.
        """
        # from
        # https://en.wikipedia.org/wiki/Plane_(geometry)#Describing_a_plane_through_three_points
        p0, p1, p2 = planepoints
        x1, y1, z1 = p0.getArray()
        x2, y2, z2 = p1.getArray()
        x3, y3, z3 = p2.getArray()
        D = np.linalg.det(np.array([[x1, y1, z1], [x2, y2, z2], [x3, y3, z3]]))
        d = -1
        at = np.linalg.det(np.array([[1, y1, z1], [1, y2, z2], [1, y3, z3]]))
        bt = np.linalg.det(np.array([[x1, 1, z1], [x2, 1, z2], [x3, 1, z3]]))
        ct = np.linalg.det(np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
        a = (-d / D) * at
        b = (-d / D) * bt
        c = (-d / D) * ct

        numer = np.abs(
            a *
            otherpoint.x +
            b *
            otherpoint.y +
            c *
            otherpoint.z +
            d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer / denom
        return dist

    def isPointToRight(self, P0, P1, P2):
        eps = 1e-6
        p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
        p1 = Vector.fromPoint(P1)
        p2 = Vector.fromPoint(P2)
        p1p0 = p1 - p0
        p2p0 = p2 - p0
        qnv = Vector.cross(p2p0, p1p0).norm()
        tmp = p0 + qnv
        tmplat, tmplon, tmpz = ecef2latlon(tmp.x, tmp.y, tmp.z)
        if (tmpz - P0.depth) < eps:
            return True
        return False

    def reverseQuad(self, P0, P1, P2, P3):
        newP0 = copy.deepcopy(P1)
        newP1 = copy.deepcopy(P0)
        newP2 = copy.deepcopy(P3)
        newP3 = copy.deepcopy(P2)
        if not self.isPointToRight(newP0, newP1, newP2):
            raise ShakeMapException(
                'Third vertex of quadrilateral must be to the right of the second vertex')
        return (newP0, newP1, newP2, newP3)

    def validateQuad(self, P0, P1, P2, P3):
        """
        Validate and fix* a given quadrilateral (*currently "fix" means check third vertex for co-planarity
        with other three points, and force it to be co-planar if it's not wildly out of the plane.()
        :param P0:
            First vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P1:
            Second vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P2:
            Third vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P3:
            Fourth vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :returns:
           Tuple of (potentially) modified vertices.
        :raises ShakeMapException:
           * if top and bottom edges are not parallel to surface
           * if dip angle is not dipping to the right relative to strike (defined by first two vertices)
           * if all 4 points are not reasonably co-planar (P2 is more than 5% of mean length of trapezoid out of plane)
        """
        # TODO: Someday fix the rule about dip angle being clockwise and 0-90 degrees
        # In theory, you could flip the quadrilateral by 180 degrees and it
        # would be ok.

        # Are the top and bottom edges both parallel to the surface?
        topDepthsEqual = P0.depth == P1.depth
        bottomDepthsEqual = P2.depth == P3.depth
        if not topDepthsEqual or not bottomDepthsEqual:
            raise ShakeMapException(
                'Top and bottom edges of fault quadrilateral must be parallel to the surface')
        # Is top edge defined by first two vertices?
        if P1.depth > P2.depth:
            raise ShakeMapException(
                'Top edge of a quadrilateral must be defined by the first two vertices')
        # Is dip angle clockwise and btw 0-90 degrees?
        if not self.isPointToRight(P0, P1, P2):
            P0, P1, P2, P3 = self.reverseQuad(P0, P1, P2, P3)
            print('Reversing quad where dip not between 0 and 90 degrees.')
        # Are all 4 points (reasonably) co-planar?
        # Translate vertices to ECEF
        p0 = Vector.fromPoint(P0)
        p1 = Vector.fromPoint(P1)
        p2 = Vector.fromPoint(P2)
        p3 = Vector.fromPoint(P3)
        # Calculate normalized vector along top edge
        v0 = (p1 - p0).norm()
        # Calculate distance btw p3 and p2
        d = (p3 - p2).mag()
        # get the new P2 value
        v1 = v0 * d
        newp2 = p3 + v1
        planepoints = [p0, p1, p2]
        dnormal = self.getDistanceToPlane(planepoints, p2)
        geometricMean = self._getTrapMeanLength(p0, p1, newp2, p3)
        if dnormal / geometricMean > OFFPLANE_TOLERANCE:
            raise ShakeMapException(
                'Points in quadrilateral are not co-planar')
        newP0 = p0.toPoint()
        newP1 = p1.toPoint()
        newP2 = newp2.toPoint()
        newP3 = p3.toPoint()
        return [newP0, newP1, newP2, newP3]

    def setQuadrilaterals(self):
        """
        Create internal list of N quadrilaterals.
        """
        # Fault QA rules
        # 1) Fault must consist of 1 or more quadrilaterals, where each quad top/bottom edges are
        #   parallel to the surface
        # 2) The strike angle of each quadrilateral is defined by the first two vertices of that quad
        # 3) The dip angle is defined by segments 2 and 3, or 1 and 4.  This angle must be clockwise with respect to the
        #   strike angle, and between 0 and 90 degrees.
        # 4) The top edge of each quad must be defined by the first two vertices of that quad.
        # 5) 4 points of quadrilateral must be co-planar
        self._lon = np.array(self._lon)
        self._lat = np.array(self._lat)
        self._depth = np.array(self._depth)
        inan = np.isnan(self._lon)
        numnans = len(self._lon[inan])
        numsegments = numnans + 1
        # requirements:
        # 1) Coordinate arrays must be same length
        # 2) Polygons must be quadrilaterals
        # 3) Quads must be closed
        # 4) Quads must be planar
        if len(self._lon) != len(self._lat) != len(self._depth):
            raise IndexError(
                'Length of input lon,lat,depth arrays must be equal')

        # Addition: self._segment_index
        #   It is also convenient to create a list of indexes for which 'trace'
        #   each quad is on (i.e, grouped). This information is required for GC2
        #   calculations. Also, in some rare situations, we need to be able to
        #   overwrite the default values.

        istart = 0
        endpoints = list(np.where(np.isnan(self._lon))[0])
        endpoints.append(len(self._lon))
        self._quadrilaterals = []
        self._segment_index = []
        segind = 0
        for iend in endpoints:
            lonseg = self._lon[istart:iend][0:-1]  # remove closing points
            latseg = self._lat[istart:iend][0:-1]
            depthseg = self._depth[istart:iend][0:-1]
            # each segment can have many contiguous quadrilaterals defined in it
            # separations (nans) between segments mean that segments are not
            # contiguous.
            npoints = len(lonseg)
            nquads = int((npoints - 4) / 2) + 1
            startidx = 0
            endidx = -1
            for i in range(0, nquads):
                topLeft = point.Point(lonseg[startidx], latseg[startidx],
                                      depthseg[startidx])
                topRight = point.Point(
                    lonseg[
                        startidx + 1],
                    latseg[
                        startidx + 1],
                    depthseg[
                        startidx + 1])
                bottomRight = point.Point(
                    lonseg[
                        endidx - 1],
                    latseg[
                        endidx - 1],
                    depthseg[
                        endidx - 1])
                bottomLeft = point.Point(lonseg[endidx], latseg[endidx],
                                         depthseg[endidx])
                surface = self.validateQuad(topLeft, topRight, bottomRight,
                                            bottomLeft)
                self._quadrilaterals.append(surface)
                startidx += 1
                endidx -= 1
            istart = iend + 1
            self._segment_index.extend([segind] * nquads)
            segind = segind + 1

    def getSegmentIndex(self):
        """
        Return a list of segment indexes.
        :returns:
            List of segment indexes; lenght equals the number of quadrilaterals.
        """
        return copy.deepcopy(self._segment_index)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            if 'xlim3d' not in list(ax.properties.keys()):
                raise ShakeMapException(
                    'Non-3d axes object passed to plot() method.')
        for quad in self._quadrilaterals:
            P0, P1, P2, P3 = quad
            ax.plot([P0.longitude], [P0.latitude], [-P0.depth], 'B.')
            ax.text([P0.longitude], [P0.latitude], [-P0.depth], 'P0')
            ax.plot([P1.longitude], [P1.latitude], [-P1.depth], 'b.')
            ax.text([P1.longitude], [P1.latitude], [-P1.depth], 'P1')
            ax.plot([P2.longitude], [P2.latitude], [-P2.depth], 'b.')
            ax.text([P2.longitude], [P2.latitude], [-P2.depth], 'P2')
            ax.plot([P3.longitude], [P3.latitude], [-P3.depth], 'b.')
            ax.text([P3.longitude], [P3.latitude], [-P3.depth], 'P3')

    def getLats(self):
        """
        Return a copy of the array of latitude points for the fault.
        :returns:
            Numpy array of latitude values.
        """
        return self._lat.copy()

    def getLons(self):
        """
        Return a copy of the array of longitude points for the fault.
        :returns:
            Numpy array of latitude values.
        """
        return self._lon.copy()

    def getReference(self):
        """
        Return whatever reference information was contained in fault file.
        :returns:
           string citeable reference
        """
        return copy.deepcopy(self._reference)

    def getNumSegments(self):
        """
        Return a count of the number of fault segments.
        :returns:
            number of fault segments
        """
        return len(np.where(np.isnan(self._lon))[0]) + 1

    def getNumQuads(self):
        """
        Return a count of the number of fault quadrilaterals.
        :returns:
            number of fault quadrilaterals.
        """
        return len(self._quadrilaterals)

    def getFaultAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.  Fault segments are separated by
        numpy.NaN values.
        :returns:
            3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.
        """
        return (np.array(self._x), np.array(self._y), np.array(self._z))

    def getFaultAsMesh(self):
        """
        Return fault segments as a OQ-Hazardlib Mesh object.
        :returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """
        fault = mesh.Mesh(self._x, self._y, self._z)
        return fault

    def _validate(self):
        """
        Ensure that all segments are closed.
        :raises ShakeMapException:
            if unclosed segments exist.
        """
        # TODO - implement ccw algorithm...
        # http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        if len(
            self._lon) != len(
            self._lat) or len(
            self._lon) != len(
                self._depth):
            raise ShakeMapException("Fault coordinates don't match")
        inan = np.where(np.isnan(self._lon))[0]
        if not len(inan):
            return
        if not np.isnan(self._lon[inan[-1]]):
            inan = list(inan).append(len(self._lon))
        istart = 0
        for i in range(0, len(inan)):
            iend = inan[i] - 1
            x1 = self._lon[istart]
            x2 = self._lon[iend]
            y1 = self._lat[istart]
            y2 = self._lat[iend]
            z1 = self._depth[istart]
            z2 = self._depth[iend]
            if x1 != x2 or y1 != y2 or z1 != z2:
                raise ShakeMapException(
                    'Unclosed segments exist in fault file.')
            istart = inan[i] + 1


def get_quad_mesh(q, dx):
    """
    Length of top eduge of a quadrilateral.
    :param q:
        A quadrilateral.
    :param dx:
        Target dx in km; used to get nx and ny of mesh, but mesh snaps
        to edges of fault so actual dx/dy will not actually equal this
        value in general.
    :returns:
        Mesh dictionary, which includes numpy arrays:
           llx: lower left x coordinate in ECEF coords.
           lly: lower left y coordinate in ECEF coords.
           llz: lower left z coordinate in ECEF coords.
           ulx: upper left x coordinate in ECEF coords.
           etc.
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

    # i and j are indices over subfaults
    ni, nj = mesh['llx'].shape
    for i in range(0, ni):
        for j in range(0, nj):
            # Fault corner points
            pp0 = Vector(
                mesh['ulx'][
                    i, j], mesh['uly'][
                    i, j], mesh['ulz'][
                    i, j])
            pp1 = Vector(
                mesh['urx'][
                    i, j], mesh['ury'][
                    i, j], mesh['urz'][
                    i, j])
            pp2 = Vector(
                mesh['lrx'][
                    i, j], mesh['lry'][
                    i, j], mesh['lrz'][
                    i, j])
            pp3 = Vector(
                mesh['llx'][
                    i, j], mesh['lly'][
                    i, j], mesh['llz'][
                    i, j])
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
    :param strike:
        Clockwise angle (deg) from north of the line at the intersection
        of the fault plane and the horizontal plane.
    :param dip:
        Angle (deg) between fault plane and the horizontal plane normal
        to the strike (0-90 using right hand rule).
    :param rake:
        Direction of motion of the hanging wall relative to the
        foot wall, as measured by the angle (deg) from the strike vector.
    :returns:
        Unit slip vector in 'local' N-S, E-W, U-D coordinates.
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
    :param strike:
        Clockwise angle (deg) from north of the line at the intersection
        of the fault plane and the horizontal plane.
    :param dip:
        Angle (deg) between fault plane and the horizontal plane normal
        to the strike (0-90 using right hand rule).
    :param rake:
        Direction of motion of the hanging wall relative to the
        foot wall, as measured by the angle (deg) from the strike vector.
    :returns:
        Unit slip vector in 'local' N-S, E-W, U-D coordinates.
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
    :param strike:
        Clockwise angle (deg) from north of the line at the intersection
        of the fault plane and the horizontal plane.
    :param dip:
        Angle (deg) between fault plane and the horizontal plane normal
        to the strike (0-90 using right hand rule).
    :param rake:
        Direction of motion of the hanging wall relative to the
        foot wall, as measured by the angle (deg) from the strike vector.
    :returns:
        Unit slip vector in 'local' N-S, E-W, U-D coordinates.
    """
    strike = np.radians(strike)
    dip = np.radians(dip)
    rake = np.radians(rake)
    sx = np.cos(rake) * np.sin(strike)
    sy = np.cos(rake) * np.cos(strike)
    sz = 0.0
    return Vector(sx, sy, sz)


def get_quad_slip(q, rake):
    """
    Compute the unit slip vector in ECEF space for a quad and rake angle.
    :param q:
        A quadrilateral.
    :param rake:
        Direction of motion of the hanging wall relative to the
        foot wall, as measured by the angle (deg) from the strike vector.
    :returns:
        Unit slip vector in ECEF space.
    """
    P0, P1, P2 = q[0:3]
    strike = P0.azimuth(P1)
    dip = getQuadDip(q)
    s1_local = getLocalUnitSlipVector(strike, dip, rake)
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
    :param q:
        A quadrilateral.
    :returns:
        Length in km.
    """
    P0, P1, P2, P3 = q
    p0 = Vector.fromPoint(P0)  # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P1)
    qlength = (p1 - p0).mag() / 1000
    return qlength


def get_quad_dip(q):
    """
    Dip of a quadrilateral.
    :param q:
        A quadrilateral.
    :returns:
        Dip in degrees.
    """
    N = get_quad_normal(q)
    V = get_vertical_vector(q)
    dip = np.degrees(np.arccos(Vector.dot(N, V)))
    return dip


def get_quad_normal(q):
    """
    Compute the unit normal vector for a quadrilateral in
    ECEF coordinates.
    :param q:
        A quadrilateral.
    :returns:
        Normalized normal vector for the quadrilateral in ECEF coords.
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
    :param q:
        A quadrilateral.
    :returns:
        The unit vector pointing in strike direction in ECEF coords.
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
    :param q:
        A quadrilateral.
    :returns:
        The unit vector pointing down dip in ECEF coords.
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
    :param q:
        A quadrilateral.
    :returns:
        Normalized vertical vector for the quadrilateral in ECEF coords.
    """
    P0, P1, P2, P3 = q
    P0_up = copy.deepcopy(P0)
    P0_up.depth = P0_up.depth - 1.0
    p0 = Vector.fromPoint(P0)   # fromPoint converts to ECEF
    p1 = Vector.fromPoint(P0_up)
    v1 = (p1 - p0).norm()
    return v1
