#!/usr/bin/env python

# stdlib modules
from abc import ABC
from abc import abstractmethod
from abc import abstractproperty
import json

# third party imports
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.gsim import base


from shakelib.plotting.plotrupture import plot_rupture_wire3d
from shakelib.plotting.plotrupture import map_rupture
from shakelib.rupture import constants


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

    def writeGeoJson(self, file):
        """
        Write the rupture to a GeoJson file.

        Args:
            file (str): Name of file.
        """
        with open(file, 'w') as f:
            json.dump(self._geojson, f)

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
        Return strike angle. If rupture consists of multiple quadrilaterals,
        the average strike angle, weighted by quad length, is returned.
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

    @abstractmethod
    def getQuadrilaterals(self):
        """
        Method to return rupture quadrilaterals. Returns None for
        PointRupture.
        """
        pass

    def getReference(self):
        """
        Returns:
           string: Reference info from file.

        """
        return self._reference

    def getOrigin(self):
        """
        Returns:
           Origin object

        """
        return self._origin

    @abstractproperty
    def lats(self):
        pass

    @abstractproperty
    def lons(self):
        pass

    @abstractproperty
    def depths(self):
        pass

    def getRuptureContext(self, gmpelist):
        """
        Returns an Openquake `RuptureContext <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.RuptureContext>`__.

        Args:
            gmpelist (list): List of hazardlib GMPE objects.

        Returns:
            RuptureContext object with all known parameters filled in.

        """  # noqa

        origin = self._origin

        # rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'

        rx = base.RuptureContext()
        rx.mag = origin.mag
        rx.strike = self.getStrike()
        rx.dip = self.getDip()
        rx.ztor = self.getDepthToTop()
        rx.width = self.getWidth()

        if hasattr(origin, 'rake'):
            rx.rake = origin.rake
        elif hasattr(origin, 'mech'):
            mech = origin.mech
            rx.rake = constants.RAKEDICT[mech]
        else:
            rx.rake = constants.RAKEDICT['ALL']

        rx.hypo_lat = origin.lat
        rx.hypo_lon = origin.lon
        rx.hypo_depth = origin.depth

        return rx

    def computeRhyp(self, lon, lat, depth):
        """
        Method for computing hypocentral distance.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Hypocentral distance (km).
        """
        origin = self._origin
        oldshape = lon.shape

        rhyp = geodetic.distance(origin.lon, origin.lat, origin.depth,
                                 lon, lat, depth)
        rhyp = rhyp.reshape(oldshape)
        return rhyp

    def computeRepi(self, lon, lat, depth):
        """
        Method for computing epicentral distance.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
           array: Epicentral distance (km).
        """
        origin = self._origin
        oldshape = lon.shape

        repi = geodetic.distance(origin.lon, origin.lat, 0.0,
                                 lon, lat, depth)
        repi = repi.reshape(oldshape)
        return repi

    @abstractmethod
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
        pass

    @abstractmethod
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
        pass

    @abstractmethod
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
        pass

    def plot3d(self):
        """
        Method for making a quick 3D wireframe plot of rupture.
        """
        plot_rupture_wire3d(self)

    def map(self):
        """
        Method for making a quick map of the fault.
        """
        map_rupture(self)
