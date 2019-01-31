#!/usr/bin/env python

# third party imports
import numpy as np
import scipy.interpolate as spint
import logging

from impactutils.time.ancient_time import HistoricTime
from shakelib.rupture.base import Rupture
from shakelib.rupture import constants

from ps2ff.constants import MagScaling, Mechanism
from ps2ff.run import single_event_adjustment


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

    def computeRjb(self, lon, lat, depth):
        """
        Convert point-distances to Joyner-Boore distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            tuple: A tuple of two arrays: first, the predicted Rjb values,
            and second an array of the variance of those predictions.
        """
        return self._computeRdist('Rjb', lon, lat, depth)

    def computeRrup(self, lon, lat, depth):
        """
        Convert point-distances to rupture distances based on magnitude.

        Args:
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            tuple: A tuple of two arrays: first, the predicted Rrup, and
            second an array of the variance of those predictions.
        """
        return self._computeRdist('Rrup', lon, lat, depth)

    def _computeRdist(self, rtype, lon, lat, depth):
        """
        Helper function to actually do the approximate fault distance
        computations.

        Args:
            type (str): Either 'Rjb' or 'Rrup'.
            lon (array): Numpy array of longitudes.
            lat (array): Numpy array of latitudes.
            depth (array): Numpy array of depths (km; positive down).

        Returns:
            tuple: A tuple of two arrays: first, the predicted approximate
            fault distance values (either Rjb or Rrup, depending on the
            'type' argument), and second an array of the variance of those
            predictions.
        """

        # ----------------------------
        # Sort out ps2ff parameters
        # ----------------------------
        origin = self._origin
        mech = getattr(origin, 'mech', 'ALL')
        if not hasattr(origin, '_tectonic_region'):
            mscale = MagScaling.WC94
            smech = Mechanism.A
            mindip_deg = 10.0
            maxdip_deg = 90.0
            aspect = 1.7
        elif origin._tectonic_region == 'Active Shallow Crust':
            mscale = MagScaling.HB08
            aspect = 1.7
            if mech == 'ALL':
                # HB08 doesn't have an 'ALL' mechanism, so use WC94
                mscale = MagScaling.WC94
                smech = Mechanism.A
                mindip_deg = 10.0
                maxdip_deg = 90.0
            elif mech == 'RS':
                smech = Mechanism.R
                mindip_deg = 35.0
                maxdip_deg = 50.0
            elif mech == 'NM':
                smech = Mechanism.N
                mindip_deg = 40.0
                maxdip_deg = 60.0
            elif mech == 'SS':
                smech = Mechanism.SS
                mindip_deg = 75.0
                maxdip_deg = 90.0
        elif origin._tectonic_region == 'Stable Shallow Crust':
            mscale = MagScaling.S14
            aspect = 1.0
            if mech == 'ALL':
                smech = Mechanism.A
                mindip_deg = 10.0
                maxdip_deg = 90.0
            elif mech == 'RS':
                smech = Mechanism.R
                mindip_deg = 30.0
                maxdip_deg = 60.0
            elif mech == 'NM':
                smech = Mechanism.N
                mindip_deg = 40.0
                maxdip_deg = 60.0
            elif mech == 'SS':
                smech = Mechanism.SS
                mindip_deg = 60.0
                maxdip_deg = 90.0
        else:
            logging.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            mscale = MagScaling.WC94
            smech = Mechanism.A
            aspect = 1.7
            mindip_deg = 10.0
            maxdip_deg = 90.0

        mindip = mindip_deg * np.pi / 180.0
        maxdip = maxdip_deg * np.pi / 180.0

        repis = np.clip(self.computeRepi(lon, lat, depth), 0.0001, None)

        repi, Rjb_hat, Rrup_hat, Rjb_var, Rrup_var = \
            single_event_adjustment(origin.mag, origin.depth, ar=aspect,
                                    mechanism=smech, mag_scaling=mscale,
                                    n_repi=13,
                                    min_repi=np.min(repis) - 1e-5,
                                    max_repi=np.max(repis) + 0.1,
                                    nxny=7, n_theta=19,
                                    n_dip=4, min_dip=mindip, max_dip=maxdip,
                                    n_eps=5, trunc=2)

        if rtype == 'Rjb':
            spline = spint.interp1d(repi, np.vstack((Rjb_hat, Rjb_var)),
                                    kind='linear', copy=False,
                                    assume_sorted=True)
            rv_hat = spline(repis)
        elif rtype == 'Rrup':
            spline = spint.interp1d(repi, np.vstack((Rrup_hat, Rrup_var)),
                                    kind='linear', copy=False,
                                    assume_sorted=True)
            rv_hat = spline(repis)
        else:
            raise ValueError('Unknown distance type in _computeRdist')
        return (rv_hat[0], rv_hat[1])

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
