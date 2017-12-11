#!/usr/bin/env python

# stdlib imports

# third party imports
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D
from mapio.geodict import GeoDict
from openquake.hazardlib.gsim.base import SitesContext
import numpy as np

# local imports
from shakelib.utils.exception import ShakeLibException


class Sites(object):
    """
    An object to encapsulate information used to generate a GEM
    `SitesContext <https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py>`__.
    """  # noqa

    def __init__(self, vs30grid, vs30measured_grid=None, backarc=None,
                 defaultVs30=686.0):
        """
        Construct a Sites object.

        Args:
            vs30grid: MapIO Grid2D object containing Vs30 values.
            vs30measured_grid: Boolean array indicating whether Vs30 values
                were measured or derived (i.e., from topographic slope).
            backarc: Boolean array indicating whether site is in the subduction
                `backarc <http://earthquake.usgs.gov/learn/glossary/?term=backarc>`__.
            defaultVs30: Default Vs30 value to use in locations where Vs30Grid
                is not specified.
        """  # noqa
        self._Vs30 = vs30grid
        if backarc is None:
            self._backarc = np.zeros_like(vs30grid.getData(), dtype=bool)
        else:
            self._backarc = backarc
            # Could add a check here that if backarc is provided that it's type
            # is bool and dimensions match vs30grid

        self._defaultVs30 = defaultVs30
        self._vs30measured_grid = vs30measured_grid
        self._GeoDict = vs30grid.getGeoDict().copy()
        self._lons = np.linspace(self._GeoDict.xmin,
                                 self._GeoDict.xmax,
                                 self._GeoDict.nx)
        self._lats = np.linspace(self._GeoDict.ymin,
                                 self._GeoDict.ymax,
                                 self._GeoDict.ny)

    @classmethod
    def _create(cls, geodict, defaultVs30, vs30File, padding, resample):
        if vs30File is not None:
            fgeodict = cls._getFileGeoDict(vs30File)
            if not resample:
                if not padding:
                    # we want something that is within and aligned
                    geodict = fgeodict.getBoundsWithin(geodict)
                else:
                    # we want something that is just aligned, since we're
                    # padding edges
                    geodict = fgeodict.getAligned(geodict)
            vs30grid = cls._load(vs30File, samplegeodict=geodict,
                                 resample=resample, method='linear',
                                 doPadding=padding, padValue=defaultVs30)

        return vs30grid

    @classmethod
    def fromBounds(cls, xmin, xmax, ymin, ymax, dx, dy, defaultVs30=686.0,
                   vs30File=None, vs30measured_grid=None,
                   backarc=None, padding=False, resample=False):
        """
        Create a Sites object by defining a center point, resolution, extent,
        and Vs30 values.

        Args:
            xmin: X coordinate of left edge of bounds.
            xmax: X coordinate of right edge of bounds.
            ymin: Y coordinate of bottom edge of bounds.
            ymax: Y coordinate of top edge of bounds.
            dx: Resolution of desired grid in X direction.
            dy: Resolution of desired grid in Y direction.
            defaultVs30: Default Vs30 value to use if vs30File not specified.
            vs30File: Name of GMT or GDAL format grid file containing Vs30
                values.
            vs30measured_grid: Boolean grid indicating whether Vs30 values were
                measured or derived (i.e., from slope).
            backarc: Boolean array indicating whether site is in the subduction
                `backarc <http://earthquake.usgs.gov/learn/glossary/?term=backarc>`__.
            padding: Boolean indicating whether or not to pad resulting Vs30
                grid out to edges of input bounds. If False, grid will be
                clipped to the extent of the input file.
            resample: Boolean indicating whether or not the grid should be
                resampled.
        """  # noqa
        geodict = GeoDict.createDictFromBox(xmin, xmax, ymin, ymax, dx, dy)
        if vs30File is not None:
            vs30grid = cls._create(geodict, defaultVs30,
                                   vs30File, padding, resample)
        else:
            griddata = np.ones((geodict.ny, geodict.nx),
                               dtype=np.float64) * defaultVs30
            vs30grid = Grid2D(griddata, geodict)
        return cls(vs30grid, vs30measured_grid=vs30measured_grid,
                   backarc=backarc, defaultVs30=defaultVs30)

    @classmethod
    def fromCenter(cls, cx, cy, xspan, yspan, dx, dy, defaultVs30=686.0,
                   vs30File=None, vs30measured_grid=None,
                   backarc=None, padding=False, resample=False):
        """
        Create a Sites object by defining a center point, resolution, extent,
        and Vs30 values.

        Args:
            cx: X coordinate of desired center point.
            cy: Y coordinate of desired center point.
            xspan: Width of desired grid.
            yspan: Height of desired grid.
            dx: Resolution of desired grid in X direction.
            dy: Resolution of desired grid in Y direction.
            defaultVs30: Default Vs30 value to use if vs30File not specified.
            vs30File: Name of GMT or GDAL format grid file containing Vs30
                values.
            vs30measured_grid: Boolean grid indicating whether Vs30 values were
                measured or derived (i.e., from slope).
            backarc: Boolean array indicating whether site is in the subduction
                `backarc <http://earthquake.usgs.gov/learn/glossary/?term=backarc>`__.
            padding: Boolean indicating whether or not to pad resulting Vs30
                grid out to edges of input bounds. If False, grid will be
                clipped to the extent of the input file.
            resample: Boolean indicating whether or not the grid should be
                resampled.
        """  # noqa
        geodict = GeoDict.createDictFromCenter(cx, cy, dx, dy, xspan, yspan)
        if vs30File is not None:
            vs30grid = cls._create(geodict, defaultVs30,
                                   vs30File, padding, resample)
        else:
            griddata = np.ones((geodict.ny, geodict.nx),
                               dtype=np.float64) * defaultVs30
            vs30grid = Grid2D(griddata, geodict)
        return cls(vs30grid, vs30measured_grid=vs30measured_grid,
                   backarc=backarc, defaultVs30=defaultVs30)

    def getSitesContext(self, lldict=None, rock_vs30=None):
        """
        Create a SitesContext object by sampling the current Sites object.

        Args:
            lldict: Either

                - None, in which case the SitesContext for the complete Sites
                  grid is returned, or
                - A location dictionary (elements are 'lats' and 'lons' and
                  each is a numpy array). Each element must have the same
                  shape. In this case the SitesContext for these locaitons is
                  returned.

            rock_vs30: Either

                - None, in which case the SitesContext will reflect the Vs30
                  grid in the Sites instance, or
                - A float for the rock Vs30 value, in which case the
                  SitesContext will be constructed for this constant Vs30
                  value.

        Returns:
            SitesContext object.

        Raises:
            ShakeLibException: When lat/lon input sequences do not share
                dimensionality.

        """  # noqa

        sctx = SitesContext()

        if lldict is not None:
            lats = lldict['lats']
            lons = lldict['lons']
            latshape = lats.shape
            lonshape = lons.shape
            if latshape != lonshape:
                msg = 'Input lat/lon arrays must have the same dimensions'
                raise ShakeLibException(msg)

            if rock_vs30 is not None:
                tmp = self._Vs30.getValue(
                    lats, lons, default=self._defaultVs30)
                sctx.vs30 = np.ones_like(tmp) * rock_vs30
            else:
                sctx.vs30 = self._Vs30.getValue(
                    lats, lons, default=self._defaultVs30)
            sctx.lats = lats
            sctx.lons = lons
        else:
            sctx.lats = self._lats.copy()
            sctx.lons = self._lons.copy()
            if rock_vs30 is not None:
                sctx.vs30 = np.ones_like(self._Vs30.getData()) * rock_vs30
            else:
                sctx.vs30 = self._Vs30.getData().copy()

        sctx = Sites._addDepthParameters(sctx)

        # For ShakeMap purposes, vs30 measured is always Fales
        sctx.vs30measured = np.zeros_like(sctx.vs30, dtype=bool)

        # Backarc should be a numpy array
        if lldict is not None:
            backarcgrid = Grid2D(self._backarc, self._Vs30.getGeoDict())
            sctx.backarc = backarcgrid.getValue(lats, lons, default=False)
        else:
            sctx.backarc = self._backarc.copy()

        return sctx

    def getVs30Grid(self):
        """

        Returns: Grid2D object containing Vs30 values for this Sites object.
        """
        return self._Vs30

    def getNxNy(self):
        """
        Returns: The number of grid points in x and y.
        """
        return self._GeoDict.nx, self._GeoDict.ny

    @staticmethod
    def _load(vs30File, samplegeodict=None, resample=False, method='linear',
              doPadding=False, padValue=np.nan):
        try:
            vs30grid = GMTGrid.load(vs30File,
                                    samplegeodict=samplegeodict,
                                    resample=resample,
                                    method=method,
                                    doPadding=doPadding,
                                    padValue=padValue)
        except Exception as msg1:
            try:
                vs30grid = GDALGrid.load(vs30File,
                                         samplegeodict=samplegeodict,
                                         resample=resample,
                                         method=method,
                                         doPadding=doPadding,
                                         padValue=padValue)
            except Exception as msg2:
                msg = 'Load failure of %s - error messages: "%s"\n "%s"' % (
                    vs30File, str(msg1), str(msg2))
                raise ShakeLibException(msg)

        if vs30grid.getData().dtype != np.float64:
            vs30grid.setData(vs30grid.getData().astype(np.float64))

        return vs30grid

    @staticmethod
    def _getFileGeoDict(fname):
        geodict = None
        try:
            geodict, t = GMTGrid.getFileGeoDict(fname)
        except Exception as msg1:
            try:
                geodict, t = GDALGrid.getFileGeoDict(fname)
            except Exception as msg2:
                msg = 'File geodict failure with %s - error messages: '\
                      '"%s"\n "%s"' % (fname, str(msg1), str(msg2))
                raise ShakeLibException(msg)
        return geodict

    @staticmethod
    def _addDepthParameters(sctx):
        """
        Add the different depth parameters to a sites context from
        Vs30 values.

        Args:
            sctx: A sites context.

        Returns: A sites context with the depth parameters set.
        """
        sctx.z1pt0_cy14_cal = Sites._z1pt0_from_vs30_cy14_cal(sctx.vs30)
        sctx.z1pt0_ask14_cal = Sites._z1pt0_from_vs30_ask14_cal(sctx.vs30)
        sctx.z2pt5_cb14_cal = Sites._z2pt5_from_vs30_cb14_cal(
            sctx.vs30) / 1000.0
        sctx.z1pt0_cy08 = Sites._z1pt0_from_vs30_cy08(sctx.vs30)
        sctx.z2pt5_cb07 = Sites._z2pt5_from_z1pt0_cb07(sctx.z1pt0_cy08)

        return sctx

    @staticmethod
    def _z1pt0_from_vs30_cy14_cal(vs30):
        """
        Compute z1.0 using CY14 relationship for California.

        Args:
            vs30: Numpy array of Vs30 values in m/s.

        Returns: Numpy array of z1.0 in m.
        """
        z1 = np.exp(-(7.15 / 4.0) *
                    np.log((vs30**4.0 + 571.**4) / (1360**4.0 + 571.**4)))
        return z1

    @staticmethod
    def _z1pt0_from_vs30_ask14_cal(vs30):
        """
        Calculate z1.0 using ASK14 relationship for California.

        Args:
            vs30: Numpy array of Vs30 values in m/s.

        Returns: Numpy array of z1.0 in m.

        """
        # ASK14 define units as km, but implemented as m in OQ
        z1 = np.exp(-(7.67 / 4.0) *
                    np.log((vs30**4.0 + 610.**4) / (1360**4.0 + 610.**4)))
        return z1

    @staticmethod
    def _z2pt5_from_vs30_cb14_cal(vs30):
        """
        Calculate z2.5 using CB14 relationship for California.

        Args:
            vs30: Numpy array of Vs30 values in m/s.

        Returns: Numpy array of z2.5 in m. *NOTE*: OQ's CampbellBozorgnia2014
            class expects z2.5 to be in km!
        """
        z2p5 = 1000 * np.exp(7.089 - 1.144 * np.log(vs30))
        return z2p5

    @staticmethod
    def _z1pt0_from_vs30_cy08(vs30):
        """
        Chiou and Youngs (2008) z1.0 equation.

        Args:
            vs30: Numpy array of Vs30 values in m/s.

        Returns: Numpy array of z1.0 in m.
        """
        z1pt0 = np.exp(28.5 - (3.82 / 8.0) * np.log(vs30**8 + 378.7**8))
        return z1pt0

    @staticmethod
    def _z2pt5_from_z1pt0_cb07(z1pt0):
        """
        Equations are from 2007 PEER report by Campbell and Bozorgnia.

        Args:
            z1pt0: Numpy array of z1.0 in m.

        Returns: Numpy array of z2.5 in m.
        """
        z2pt5 = 519.0 + z1pt0 * 3.595
        return z2pt5
