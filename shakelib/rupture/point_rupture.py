#!/usr/bin/env python

# stdlib modules
import os
import warnings
import re

# third party imports
import numpy as np
import pandas as pd
import scipy.interpolate as spint

from impactutils.time.ancient_time import HistoricTime
from shakelib.rupture.base import Rupture
from shakelib.rupture import utils
from shakelib.rupture import constants


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
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
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
        cdir, tmp = os.path.split(__file__)
        origin = self._origin

        # -------------------
        # Sort out file names
        # -------------------
        mech = origin.mech
        if not hasattr(origin, '_tectonic_region'):
            rf = os.path.join(
                cdir, "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
            vf = os.path.join(
                cdir, "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")
        elif origin._tectonic_region == 'Active Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechA_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechA_ar1p7_seis0_20_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechR_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechR_ar1p7_seis0_20_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechN_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechN_ar1p7_seis0_20_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechSS_ar1p7_seis0_20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_WC94_mechSS_ar1p7_seis0_20_Var.csv")
        elif origin._tectonic_region == 'Stable Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechA_ar1p0_seis0_15_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechR_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechR_ar1p0_seis0_15_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechN_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechN_ar1p0_seis0_15_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechSS_ar1p0_seis0_15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rjb_S14_mechSS_ar1p0_seis0_15_Var.csv")
        else:
            warnings.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            rf = os.path.join(
                cdir, "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
            vf = os.path.join(
                cdir, "ps2ff",
                "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")

        # -----------------
        # Start with ratios
        # -----------------
        repi2rjb_ratios_tbl = pd.read_csv(rf, comment='#')
        r2rrt_cols = repi2rjb_ratios_tbl.columns[1:]
        mag_list = []
        for column in (r2rrt_cols):
            if re.search('R\d+\.*\d*', column):
                magnitude = float(re.findall(
                    'R(\d+\.*\d*)', column)[0])
                mag_list.append(magnitude)
        mag_list = np.array(mag_list)
        dist_list = np.log(np.array(repi2rjb_ratios_tbl['Repi_km']))
        repi2rjb_grid = repi2rjb_ratios_tbl.values[:, 1:]
        repi2rjb_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rjb_grid, kx=1, ky=1)

        def repi2rjb_tbl(repi, M):
            ratio = repi2rjb_obj.ev(np.log(repi), M)
            rjb = repi * ratio
            return rjb

        repis = self.computeRepi(lon, lat, depth)
        mags = np.ones_like(repis) * origin.mag
        rjb_hat = repi2rjb_tbl(repis, mags)

        # -------------------
        # Additional Variance
        # -------------------
        repi2rjbvar_ratios_tbl = pd.read_csv(vf, comment='#')
        repi2rjbvar_grid = repi2rjbvar_ratios_tbl.values[:, 1:]
        repi2rjbvar_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rjbvar_grid, kx=1, ky=1)
        rjbvar = repi2rjbvar_obj.ev(np.log(repis), mags)

        if var is True:
            return (rjb_hat, rjbvar)
        else:
            return rjb_hat

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
                predicted Rjb values, and second an array of the variance of
                those predictions. If var is False then just the first element
                of the tuple is returned.
        """
        cdir, tmp = os.path.split(__file__)
        origin = self._origin

        # -------------------
        # Sort out file names
        # -------------------
        rake = float(origin.rake)
        mech = utils.rake_to_mech(rake)
        if not hasattr(origin, '_tectonic_region'):
            rf = os.path.join(
                cdir, "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
            vf = os.path.join(
                cdir, "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")
        elif origin._tectonic_region == 'Active Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechA_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechA_ar1p7_seis0-20_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechR_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechR_ar1p7_seis0-20_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechN_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechN_ar1p7_seis0-20_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechSS_ar1p7_seis0-20_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_WC94_mechSS_ar1p7_seis0-20_Var.csv")
        elif origin._tectonic_region == 'Stable Shallow Crust':
            if mech == 'ALL':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechA_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechA_ar1p0_seis0-15_Var.csv")
            elif mech == 'RS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechR_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechR_ar1p0_seis0-15_Var.csv")
            elif mech == 'NM':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechN_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechN_ar1p0_seis0-15_Var.csv")
            elif mech == 'SS':
                rf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechSS_ar1p0_seis0-15_Ratios.csv")
                vf = os.path.join(
                    cdir, "ps2ff",
                    "Rrup_S14_mechSS_ar1p0_seis0-15_Var.csv")
        else:
            warnings.warn(
                'Unsupported tectonic region; using coefficients for unknown'
                'tectonic region.')
            rf = os.path.join(
                cdir, "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
            vf = os.path.join(
                cdir, "ps2ff",
                "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")

        # -----------------
        # Start with ratios
        # -----------------
        repi2rrup_ratios_tbl = pd.read_csv(rf, comment='#')
        r2rrt_cols = repi2rrup_ratios_tbl.columns[1:]
        mag_list = []
        for column in (r2rrt_cols):
            if re.search('R\d+\.*\d*', column):
                magnitude = float(re.findall(
                    'R(\d+\.*\d*)', column)[0])
                mag_list.append(magnitude)
        mag_list = np.array(mag_list)
        dist_list = np.log(np.array(repi2rrup_ratios_tbl['Repi_km']))
        repi2rrup_grid = repi2rrup_ratios_tbl.values[:, 1:]
        repi2rrup_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rrup_grid, kx=1, ky=1)

        def repi2rrup_tbl(repi, M):
            ratio = repi2rrup_obj.ev(np.log(repi), M)
            rrup = repi * ratio
            return rrup

        repis = self.computeRepi(lon, lat, depth)
        mags = np.ones_like(repis) * origin.mag
        rrup_hat = repi2rrup_tbl(repis, mags)

        # -------------------
        # Additional Variance
        # -------------------
        repi2rrupvar_ratios_tbl = pd.read_csv(vf, comment='#')
        repi2rrupvar_grid = repi2rrupvar_ratios_tbl.values[:, 1:]
        repi2rrupvar_obj = spint.RectBivariateSpline(
            dist_list, mag_list, repi2rrupvar_grid, kx=1, ky=1)
        rrupvar = repi2rrupvar_obj.ev(np.log(repis), mags)

        if var is True:
            return (rrup_hat, rrupvar)
        else:
            return rrup_hat

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
        # gracefully as used in GMPEs.
        dict = {"rx": np.zeros_like(lon),
                "ry": np.zeros_like(lon),
                "ry0": np.zeros_like(lon),
                "U": np.zeros_like(lon),
                "T": np.zeros_like(lon)
                }
        return dict
