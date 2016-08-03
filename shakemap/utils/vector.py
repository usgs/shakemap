#!/usr/bin/env python

import numpy as np
from openquake.hazardlib.geo import point
from ..utils.ecef import latlon2ecef, ecef2latlon


class Vector(object):
    """
    Three-dimensional vector object, stored as three floats of x,y,z.

    Todo: 
        Optimize/vectorize calculations like dot/cross for arrays of
        vectors. 
    """

    def __init__(self, x, y, z):
        """
        Create three dimensional vector object in cartesian space.

        :param x:
            x coordinate (float). 
        :param y:
            y coordinate (float). 
        :param z:
            z coordinate (float). 
        :returns:
            Vector object containing x,y,z coordinates as floats.
        """
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @classmethod
    def fromPoint(cls, oqpoint):
        """
        Class method which allows user to create a Vector from a GEM Hazardlib 
        Point object.  The Point lat, lon, depth values are converted to 
        Earth-Centered-Earth-Fixed (ECEF) cartesian coordinates.

        :param oqpoint:
            Openquake 
            `Point <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py>`))
            object. 
        :returns:
            A Vector object.
        """
        x, y, z = latlon2ecef(
            oqpoint.latitude, oqpoint.longitude, oqpoint.depth)
        return Vector(x, y, z)

    @classmethod
    def fromTuple(cls, a):
        """
        Class method which allows user to create a Vector from an x/y/z tuple.

        :param a:
            an x/y/z tuple.
        :returns:
            a Vector object.
        """
        x, y, z = a
        return Vector(x, y, z)

    def __add__(self, other):
        """
        Add another Vector object to this one (x+x,y+y,z+z).

        :param other:
            Another Vector object
        :returns:
            A third Vector object.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError("Cannot add Vector and %s objects" % type(other))
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        """
        Subtract another Vector object from this one (x+x,y+y,z+z).

        :param other:
            Another Vector object.
        :returns:
            A third Vector object.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "Cannot subtract Vector and %s objects" %
                type(other))
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, length):
        """
        Multiply the Vector by a scalar, changing it's length.

        :param length:
            A scalar number.
        :returns:
            A Vector object.
        :raises TypeError:
            when length is not a number.
        """
        try:
            length = float(length)
        except ValueError:
            raise TypeError(
                "Cannot multiply Vector and %s objects" %
                type(length))
        return Vector(self.x * length, self.y * length, self.z * length)

    def __rmul__(self, length):
        """
        Multiply the Vector by a scalar, changing it's length.

        :param length:
            A scalar number.
        :returns:
            A Vector object.
        :raises TypeError:
            When length is not a number.
        """
        try:
            length = float(length)
        except ValueError:
            raise TypeError(
                "Cannot multiply Vector and %s objects" %
                type(length))
        return Vector(self.x * length, self.y * length, self.z * length)

    def __eq__(self, other):
        """
        Check equality between this Vector and another.

        :param other:
            Another Vector object.
        :returns:
            True or False.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "Cannot compare Vector and %s objects" %
                type(other))
        if other.x == self.x and other.y == self.y and other.z == self.z:
            return True
        return False

    def distance(self, other):
        """
        Calculate distance between this Vector and another.

        :param other:
            Another Vector object.
        :returns:
            float distance between Vectors.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "Cannot calculate distance between Vector and %s objects" %
                type(other))
        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)
                       ** 2 + (self.z - other.z)**2)

    def cross(self, other):
        """
        Calculate cross product between this Vector and another.

        :param other:
            Another Vector object.
        :returns:
            a Vector object.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "Cannot calculate cross product between Vector and %s objects" %
                type(other))
        cp = np.cross(self.getArray(), other.getArray())
        return Vector(cp[0], cp[1], cp[2])

    def dot(self, other):
        """
        Calculate dot product between this Vector and another.

        :param other:
           Another Vector object.
        :returns:
           a float dot product.
        :raises TypeError:
            If other is not a Vector object.
        """
        if not isinstance(other, Vector):
            raise TypeError(
                "Cannot calculate cross product between Vector and %s objects" %
                type(other))
        dp = np.dot(self.getArray(), other.getArray())
        return dp

    def getArray(self):
        """
        :returns:
            3 element Numpy array of [x,y,z]
        """
        return np.array((self.x, self.y, self.z))

    def getTuple(self):
        """
        :returns:
            3 element tuple of (x,y,z)
        """
        return (self.x, self.y, self.z)

    def norm(self):
        """
        :returns:
            Normalized Vector.
        """
        length = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        x = self.x / length
        y = self.y / length
        z = self.z / length
        return Vector(x, y, z)

    def mag(self):
        """
        :returns:
            Length of Vector (float).
        """
        length = np.sqrt(self.x**2 + self.y**2 + self.z**2)
        return length

    def toPoint(self):
        """
        Convert the Vector to a hazardlib Point object, after translating 
        back to lat, lon, depth.

        :returns:
            An Openquake 
            `Point <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py>`__.
        """
        lat, lon, dep = ecef2latlon(self.x, self.y, self.z)
        return point.Point(lon, lat, dep)

    def __repr__(self):
        """
        String representation of Vector.
        """
        return '<x=%.4f,y=%.4f,z=%.4f>' % (self.x, self.y, self.z)
