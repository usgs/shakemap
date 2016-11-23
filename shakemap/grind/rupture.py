#!/usr/bin/env python

# stdlib modules
import copy
from abc import ABC
from abc import abstractmethod
import json

# third party imports
import numpy as np
from openquake.hazardlib.geo import mesh
from openquake.hazardlib.geo import point
from openquake.hazardlib.geo.utils import get_orthographic_projection
from ..utils.ecef import latlon2ecef
from ..utils.ecef import ecef2latlon
from ..utils.vector import Vector

# local imports
from shakemap.utils.exception import ShakeMapException

# CONSTANTS
# what is the maximum ratio of distance out of the plane defined by 3 points a
# 4th point can be before being considered non-co-planar?
OFFPLANE_TOLERANCE = 0.05

class Rupture(ABC):
    """
    Abstract base class for ruptures.
    """

    @abstractmethod
    def getRuptureLength(self):
        pass

    @abstractmethod
    def getStrike(self):
        pass

    @abstractmethod
    def getTopOfRupture(self):
        pass

    @abstractmethod
    def getDip(self):
        pass

    @abstractmethod
    def getWidth(self):
        pass

    ## Does this need to be an abstract method?? 
    def getRuptureContext(self, gmpelist):
        """
        Return an Openquake 
        `RuptureContext <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.RuptureContext>`__.

        If Origin does not contain a Rupture, then strike, dip, ztor, and width
        will be filled with default values. Rake may not be known, or may be
        estimated from a focal mechanism.

        Args:
            gmpelist (list): List of hazardlib GMPE objects.

        Returns:
            RuptureContext object with all known parameters filled in.

        """
        # rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'


        rup = base.RuptureContext()
        rup.mag = self.getEventParam('mag')
        if self._rupture is not None:
            rup.strike = self._rupture.getStrike()
            rup.dip = self._rupture.getDip()
            rup.ztor = self._rupture.getTopOfRupture()
            rup.width = self._rupture.getWidth()
        else:
            rup.strike = DEFAULT_STRIKE
            rup.dip = DEFAULT_DIP
            rup.ztor = DEFAULT_ZTOR
            rup.width = DEFAULT_WIDTH

        if 'rake' in self._event_dict:
            rup.rake = self.getEventParam('rake')
        elif 'mech' in self._event_dict:
            mech = self._event_dict['mech']
            rup.rake = RAKEDICT[mech]
        else:
            rup.rake = RAKEDICT['ALL']

        rup.hypo_lat = self.getEventParam('lat')
        rup.hypo_lon = self.getEventParam('lon')
        rup.hypo_depth = self.getEventParam('depth')

        if rup.rake is None:
            rup.rake = DEFAULT_RAKE

        return rup
    



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
                 segment_index = None, reference = ''):
        """
        Constructor for EdgeRupture class. 

        Args:
            toplons (ndarray): Array of top edge longitudes.
            toplats (ndarray): Array of top edge latitudes. 
            topdeps (ndarray): Array of top edge depths (km).
            botlons (ndarray): Array of bot edge longitudes.
            botlats (ndarray): Array of bot edge latitudes. 
            botdeps (ndarray): Array of bot edge depths (km).
            segment_index (ndarray): Optional array of segment index. 
                If None, then assume only single segment. 
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
        if segment_index is not None:
            self._segment_index = segment_index
        else:
            self._segment_index = np.zeros_like(toplons)
        self._reference = reference
        # TODO: add validation
        #      - Check that arrays are the same length
        #      - Check that if segment index is supplied that it makes sens
        #        (integers, sequential, etc)
        #      - ???


    @classmethod
    def readJsonFile(cls, filename):
        """
        Method for reading a JSON file that specifies the top and bottom edges
        of the rupture. 

        Args: 
            filename (str): Name of JSON file. 

        Returns: 
            EdgeRupture instance.

        """
        with open(filename) as f:
            rj = json.load(f)
        toplats = np.array(rj['toplats'])
        toplons = np.array(rj['toplons'])
        topdeps = np.array(rj['topdeps'])
        botlats = np.array(rj['botlats'])
        botlons = np.array(rj['botlons'])
        botdeps = np.array(rj['botdeps'])
        reference = rj['reference']
        if 'segment index' in rj.keys():
            segment_index = rj['segment index']
        else:
            segment_index = np.zeros_like(toplats)

        return cls(toplons, toplats, topdeps, botlons, botlats, botdeps,
                   segment_index = segment_index, reference = reference)

    def writeRuptureFile(self):
        pass

    def getRuptureLength(self):
        pass

    def getQuadrilaterals(self):
        pass

    def getStrike(self):
        pass

    def getTopOfRupture(self):
        pass

    def getDip(self):
        pass

    def getWidth(self):
        pass

    def getIndividualWidths(self):
        pass

    def getIndividualTopLengths(self):
        pass

    def getLats(self):
        pass

    def getLons(self):
        pass

    def getDeps(self):
        pass

    def getNumSegments(self):
        pass

    def getNumQuads(self):
        pass



class QuadRupture(Rupture):
    """
    Rupture class that represents the rupture surface as a combination of 
    quadrilaterals. Each quadrilateral must have horizontal top and bottom
    edges. This restriction makes the computation of rupture distances 
    more efficient. The number of points in the top edges must match the 
    number of points in the bottom edge. 
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
        self._validate()
        self._setQuadrilaterals()

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

        # assemble the vertices as the constructor needs them...
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
    def readTextFile(cls, rupturefile):
        """
        Read rupture file format as defined in ShakeMap 3.5 Software Guide.

        Args:
            rupturefile (srt): Path to rupture file OR file-like object in GMT
                psxy format, where:

                * Rupture vertices are space separated lat,lon,depth triplets on
                  a single line.
                * Rupture segments are separated by lines containing ">"
                * Rupture segments must be closed.
                * Rupture segments must be all clockwise or all 
                  counter-clockwise.

        Returns:
           QuadRupture object.

        """
        x = []
        y = []
        z = []
        isFile = False
        if isinstance(rupturefile, str):
            isFile = True
            rupturefile = open(rupturefile, 'rt')
            rupturelines = rupturefile.readlines()
        else:
            rupturelines = rupturefile.readlines()
        reference = ''
        for line in rupturelines:
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
                    'Finite rupture file %s has no depth values.' %
                    rupturefile)
            y.append(float(parts[0]))
            x.append(float(parts[1]))
            z.append(float(parts[2]))
        if isFile:
            rupturefile.close()
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

        # assemble the vertices as the constructor needs them...
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

    def getRuptureLength(self):
        """
        Compute lenght of rupture based on top edge in km.

        Returns:
            float: Length of rupture (km).

        """
        flength = 0
        for quad in self._quadrilaterals:
            flength = flength + get_quad_length(quad)
        return flength

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
            float: Strike angle.

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

        :returns:
           Average dip in degrees (float).
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
        Return the average rupture width (km) for all quadrilaterals defined for
        the rupture.

        :returns:
            Average width in km of all rupture quadrilaterals (float).
        """
        wsum = 0.0
        for quad in self._quadrilaterals:
            P0, P1, P2, P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p3 = Vector.fromPoint(P3)
            wsum += get_quad_width(p0, p1, p3)
        mwidth = (wsum / len(self._quadrilaterals)) / 1000.0
        return mwidth

    def getIndividualWidths(self):
        """
        Return an array of rupture widths (km), one for each quadrilateral
        defined for the rupture.

        :returns:
            Array of quad widths in km of all rupture quadrilaterals.
        """
        nquad = self.getNumQuads()
        widths = np.zeros(nquad)
        for i in range(nquad):
            P0, P1, P2, P3 = self._quadrilaterals[i]
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p3 = Vector.fromPoint(P3)
            widths[i] = get_quad_width(p0, p1, p3) / 1000.0
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
    def _getTrapMeanLength(p0, p1, p2, p3):
        """
        Return the sqrt of the area of a quadrilateral (used for QA of rupture 
        plane).

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
        h = get_quad_width(p0, p1, p3)
        a = (p1 - p0).mag()
        b = (p2 - p3).mag()
        A = ((a + b) / 2.0) * h
        length = np.sqrt(A)
        return length

    @staticmethod
    def getDistanceToPlane(planepoints, otherpoint):
        """
        Calculate a point's distance to a plane.  Used to figure out if a
        quadrilateral points are all co-planar.

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
        if D != 0:
            d = -1
            at = np.linalg.det(
                np.array([[1, y1, z1], [1, y2, z2], [1, y3, z3]]))
            bt = np.linalg.det(
                np.array([[x1, 1, z1], [x2, 1, z2], [x3, 1, z3]]))
            ct = np.linalg.det(
                np.array([[x1, y1, 1], [x2, y2, 1], [x3, y3, 1]]))
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

    @staticmethod
    def _isPointToRight(P0, P1, P2):
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

    def _reverseQuad(self, P0, P1, P2, P3):
        newP0 = copy.deepcopy(P1)
        newP1 = copy.deepcopy(P0)
        newP2 = copy.deepcopy(P3)
        newP3 = copy.deepcopy(P2)
        if not self._isPointToRight(newP0, newP1, newP2):
            raise ShakeMapException(
                'Third vertex of quadrilateral must be to the right of the second vertex')
        return (newP0, newP1, newP2, newP3)

    def _validateQuad(self, P0, P1, P2, P3):
        """
        Validate and fix* a given quadrilateral (*currently "fix" means check
        third vertex for co-planarity with other three points, and force it to
        be co-planar if it's not wildly out of the plane.()

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
           * if dip angle is not dipping to the right relative to strike 
             (defined by first two vertices)
           * if all 4 points are not reasonably co-planar (P2 is more than 5% of
             mean length of trapezoid out of plane)
        """
        # TODO: Someday fix the rule about dip angle being clockwise and 0-90 degrees
        # In theory, you could flip the quadrilateral by 180 degrees and it
        # would be ok.

        # Are the top and bottom edges both parallel to the surface?
        topDepthsEqual = np.allclose(P0.depth, P1.depth, atol = 3e-3)
        bottomDepthsEqual = np.allclose(
            P2.depth, P3.depth, atol = 5e-2, rtol = 1e-4)
        if not topDepthsEqual or not bottomDepthsEqual:
            raise ShakeMapException(
                'Top and bottom edges of rupture quadrilateral '\
                'must be parallel to the surface')

        # Is top edge defined by first two vertices?
        if P1.depth > P2.depth:
            raise ShakeMapException(
                'Top edge of a quadrilateral must be defined by '\
                'the first two vertices')

        # Is dip angle clockwise and btw 0-90 degrees?
        if not self._isPointToRight(P0, P1, P2):
            P0, P1, P2, P3 = self._reverseQuad(P0, P1, P2, P3)
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

    def _setQuadrilaterals(self):
        """
        Create internal list of N quadrilaterals.
        """
        # QuadRupture QA rules
        # 1) Rupture must consist of 1 or more quadrilaterals, where each quad
        #    top/bottom edges are parallel to the surface
        # 2) The strike angle of each quadrilateral is defined by the first two
        #    vertices of that quad
        # 3) The dip angle is defined by segments 2 and 3, or 1 and 4.  This
        #    angle must be clockwise with respect to the strike angle, and
        #    between 0 and 90 degrees.
        # 4) The top edge of each quad must be defined by the first two vertices
        #    of that quad.
        # 5) 4 points of quadrilateral must be co-planar
        self._lon = np.array(self._lon)
        self._lat = np.array(self._lat)
        self._depth = np.array(self._depth)
        inan = np.isnan(self._lon)
        numnans = len(self._lon[inan])

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
                    lonseg[startidx + 1],
                    latseg[startidx + 1],
                    depthseg[startidx + 1])
                bottomRight = point.Point(
                    lonseg[endidx - 1],
                    latseg[endidx - 1],
                    depthseg[endidx - 1])
                bottomLeft = point.Point(lonseg[endidx], latseg[endidx],
                                         depthseg[endidx])
                surface = self._validateQuad(topLeft, topRight, bottomRight,
                                             bottomLeft)
                self._quadrilaterals.append(surface)
                startidx += 1
                endidx -= 1
            istart = iend + 1
            self._segment_index.extend([segind] * nquads)
            segind = segind + 1

    def _getSegmentIndex(self):
        """
        Return a list of segment indexes.

        :returns:
            List of segment indexes; lenght equals the number of quadrilaterals.
        """
        return copy.deepcopy(self._segment_index)

    def getLats(self):
        """
        Return a copy of the array of latitudes for the rupture verticies.

        :returns:
            Numpy array of latitude values.
        """
        return self._lat.copy()

    def getLons(self):
        """
        Return a copy of the array of longitudes for the rupture verticies.

        :returns:
            Numpy array of latitude values.
        """
        return self._lon.copy()

    def getDeps(self):
        """
        Return a copy of the array of depths for the rupture verticies.

        :returns:
            Numpy array of latitude values.
        """
        return self._depth.copy()

    def getNumSegments(self):
        """
        Return a count of the number of rupture segments.

        :returns:
            number of rupture segments
        """
        return len(np.where(np.isnan(self._lon))[0]) + 1

    def getNumQuads(self):
        """
        Return a count of the number of rupture quadrilaterals.

        :returns:
            number of rupture quadrilaterals.
        """
        return len(self._quadrilaterals)

    def getRuptureAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X, Y, Z (lon,lat,depth)
        coordinates. Rupture segments are separated by numpy.NaN values.

        :returns:
            3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.
        """
        return (np.array(self._lon), np.array(self._lat), np.array(self._depth))

    def getRuptureAsMesh(self):
        """
        Return rupture segments as a OQ-Hazardlib Mesh object.

        :returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """
        rupture = mesh.Mesh(self._lon, self._lat, self._depth)
        return rupture

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
            raise ShakeMapException("Rupture coordinates don't match")
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
                    'Unclosed segments exist in rupture file.')
            istart = inan[i] + 1


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
        if isinstance(file, str):
            with open(file) as f:
                rj = json.load(f)
        else:
            rj = json.loads(str(file))
        if rj['rupture type'] == 'EdgeRupture':
            rupt = EdgeRupture.readJsonFile(file)
        else:
            rupt = QuadRupture.readJsonFile(file)
    except json.JSONDecodeError:
        try:
            rupt = QuadRupture.readTextFile(file)
        except:
            raise Exception("Unknown rupture file format.")
    return rupt


def get_quad_width(p0, p1, p3):
    """
    Return width of an individual planar trapezoid, where the p0-p1 distance
    represents the long side.

    Args:
        p0 (point): ECEF x,y,z point representing the first vertex of a 
            quadrilateral.
        p1 (point): ECEF x,y,z point representing the second vertex of a 
            quadrilateral.
        p3 (point):  ECEF x,y,z point representing the fourth vertex of a 
            quadrilateral.

    Returns:
        float: Width of planar trapezoid.
    """
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
