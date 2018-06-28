#!/usr/bin/env python

# third party imports
import numpy as np
import logging

from impactutils.time.ancient_time import HistoricTime
from shakelib.rupture.base import Rupture
from shakelib.rupture import constants

from ps2ff.constants import DistType, MagScaling, Mechanism
from ps2ff.interpolate import PS2FF


class PointRupture(Rupture):
    """
    Rupture class for point sources. The purpose is to gracefully handle:

        - Requests for rupture distances when no rupture is available.
        - Provide reasonable default values for rupture parameters.
    """

    def __init__(self, origin, reference="Origin"):
        """
        Constructs a PointRupture instance.

        Args:
            origin (Origin): Reference to a ShakeMap Origin instance.
            reference (str): Citable reference for rupture; in the case of a
                PointRupture, the 'reference' is probably the origin.

        Returns:
            PointRupture instance.
        """
        self._origin = origin
        self._reference = reference

        coords = [origin.lon, origin.lat, origin.depth]

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
                     "type": "Point",
                     "coordinates": coords
                 }
             }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime(constants.TIMEFMT)
            else:
                d['metadata'][k] = v

        self._geojson = d

    def getLength(self):
        """
        Rupture length, which is None for a PointRupture.
        Could potentially put in a default value based on magnitude.
        """
        return None

    def getWidth(self):
        """
        Rupture width.
        Could potentially put in a default value based on magnitude.
        """
        return constants.DEFAULT_WIDTH

    def getArea(self):
        """
        Rupture area, which is None for a PointRupture.
        Could potentially put in a default value based on magnitude.
        """
        return None

    def getStrike(self):
        """
        Strike, which is None.
        Could potentially get from strec or something?
        """
        return constants.DEFAULT_STRIKE

    def getDip(self):
        """
        Dip, which is None.
        Could potentially get from strec or something?
        """
        return constants.DEFAULT_DIP

    def getDepthToTop(self):
        """
        Depth to top of rupture.
        Could get from hypo/magnitude?
        """
        return constants.DEFAULT_ZTOR

    def getQuadrilaterals(self):
        return None

    @property
    def lats(self):
        """
        Returns rupture latitudes, which is just the hypocenter for a
        PointRupture."""
        return self._origin.lat

    @property
    def lons(self):
        """
        Returns rupture longitudes, which is just the hypocenter for a
        PointRupture."""
        return self._origin.lon

    @property
    def depths(self):
        """
        Returns rupture depths, which is just the hypocenter for a
        PointRupture."""
        return self._origin.depth

    def computeRjb(self, lon, lat, depth, var=False):
        """
        Convert point-distances to Joyner-Boore distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).
            var (bool): Also return variance of prediction.

        Returns:
            If var is True then this returns a tuple of two arrays: first, the
                predicted Rjb values, and second an array of the variance of
                those predictions. If var is False then just the first element
                of the tuple is returned.
        """
        return self._computeRdist('Rjb', lon, lat, depth, var)

    def computeRrup(self, lon, lat, depth, var=False):
        """
        Convert point-distances to rupture distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).
            var (bool): Also return variance of prediction.

        Returns:
            If var is True then this returns a tuple of two arrays: first, the
                predicted Rrup, and second an array of the variance of
                those predictions. If var is False then just the first element
                of the tuple is returned.
        """
        return self._computeRdist('Rrup', lon, lat, depth, var)

    def _computeRdist(self, rtype, lon, lat, depth, var):
        """
        Helper function to actually do the approximate fault distance
        computations.

        Args:
            type (str): Either 'Rjb' or 'Rrup'.
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).
            var (bool): Also return variance of prediction.

        Returns:
            If var is True then this returns a tuple of two arrays: first, the
                predicted approximate fault distance values (either Rjb or
                Rrup, depending on the 'type' argument), and second an array of
                the variance of those predictions. If var is False then just
                the first element of the tuple is returned.
        """
        if rtype == 'Rjb':
            dtype = DistType.Rjb
        elif rtype == 'Rrup':
            dtype = DistType.Rrup
        else:
            raise ValueError('Unknown distance type in _computeRdist')

        # ----------------------------
        # Sort out ps2ff parameters
        # ----------------------------
        origin = self._origin
        mech = origin.mech
        if not hasattr(origin, '_tectonic_region'):
            mscale = MagScaling.WC94
            smech = Mechanism.A
            aspect = 1.7
            min_sdepth = 0
            max_sdepth = 20
        elif origin._tectonic_region == 'Active Shallow Crust':
            mscale = MagScaling.HB08
            aspect = 1.7
            min_sdepth = 0
            max_sdepth = 20
            if mech == 'ALL':
                # HB08 doesn't have an 'ALL' mechanism, so use WC94
                mscale = MagScaling.WC94
                smech = Mechanism.A
            elif mech == 'RS':
                smech = Mechanism.R
            elif mech == 'NM':
                smech = Mechanism.N
            elif mech == 'SS':
                smech = Mechanism.SS
        elif origin._tectonic_region == 'Stable Shallow Crust':
            mscale = MagScaling.S14
            aspect = 1.0
            min_sdepth = 0
            max_sdepth = 15
            if mech == 'ALL':
                smech = Mechanism.A
            elif mech == 'RS':
                smech = Mechanism.R
            elif mech == 'NM':
                smech = Mechanism.N
            elif mech == 'SS':
                smech = Mechanism.SS
        else:
            logging.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            mscale = MagScaling.WC94
            smech = Mechanism.A
            aspect = 1.7
            min_sdepth = 0
            max_sdepth = 20

        p2f = PS2FF.fromParams(dist_type=dtype,
                               mag_scaling=mscale,
                               mechanism=smech,
                               AR=aspect,
                               min_seis_depth=min_sdepth,
                               max_seis_depth=max_sdepth)

        repis = np.clip(self.computeRepi(lon, lat, depth), 0.0001, None)
        mags = np.full_like(repis, origin.mag)

        r_hat = p2f.r2r(repis, mags)

        if var is True:
            r_var = p2f.var(repis, mags)
            return (r_hat, r_var)
        else:
            return r_hat

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
        # This just returns defaults of zero, which will hopefully behave
        # gracefully as used in GMPEs but should eventually be updated
        # so that things like the hangingwall terms are unbiased.
        gc2 = {
            "rx": np.zeros_like(lon),
            "ry": np.zeros_like(lon),
            "ry0": np.zeros_like(lon),
            "U": np.zeros_like(lon),
            "T": np.zeros_like(lon)
        }
        return gc2
