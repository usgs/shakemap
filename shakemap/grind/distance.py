#!/usr/bin/env python

# stdlib imports
import copy
import warnings
import itertools as it
import os

# third party imports
import numpy as np
import pandas as pd
import re
import scipy.interpolate as spint

from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo.utils import get_orthographic_projection
from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.gsim import base

# local imports
from shakemap.utils.exception import ShakeMapException
from shakemap.utils.ecef import latlon2ecef
from shakemap.utils.vector import Vector
from shakemap.grind.rupture import QuadRupture
from shakemap.grind.rupture import EdgeRupture
from shakemap.grind.rupture import PointRupture
from shakemap.grind.rupture import get_quad_length
from shakemap.grind.rupture import reverse_quad


class Distance(object):
    """
    Class for distance calculations. Primary method is 'get_distance'. 
    To gracefully handle multiple segment ruptures, many of the distances
    are based on the Spudich and Chiou (2015) GC2 coordinate system. 


    References: 
        Spudich, Paul and Chiou, Brian, 2015, Strike-parallel and strike-normal
        coordinate system around geometrically complicated rupture traces—Use by
        NGA-West2 and further improvements: U.S. Geological Survey Open-File
        Report 2015-1028, 20 p., http://dx.doi.org/10.3133/ofr20151028.
    """

    def __init__(self, gmpe, lon, lat, dep, origin = None, rupture = None):
        """
        Constructor for Distance class. 

        Args:
            gmpe (GMPE): Concrete subclass of GMPE
                `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
                can be individual instance or list of instances.
            lon (array): A numpy array of site longitudes.
            lat (array): A numpy array of site latitudes.
            dep (array): A numpy array of site depths (km); down is positive. 
            origin (Origin): Shakemap Origin instance. For point-source 
                distances (Repi, Rhyp) the hypocenter is taken from the origin
                instance. The finite-rupture distances (Rrup, Rjb, ...), are
                based on the rupture representation in the source instance if
                available.
            rupture (Rupture): A Shakemap Rupture instance.

        Returns:
            Distance object.
        """
        self._origin = origin
        self._rupture = rupture

        self._distance_context = self._calcDistanceContext(
            gmpe, lat, lon, dep)

    @classmethod
    def fromSites(cls, gmpe, origin, sites, rup = None):
        """
        Convenience class method to construct a Distance object from a sites
        object.

        Args:
            gmpe (GMPE): Concrete subclass of GMPE
                `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
                can be individual instance or list of instances.
            origin (Origin): A Shakeamp Origin object.
            sites (Sites): A Shakemap Sites object.

        Returns:
            Distance object.
        """
        sm_dict = sites._GeoDict
        west = sm_dict.xmin
        east = sm_dict.xmax
        south = sm_dict.ymin
        north = sm_dict.ymax
        nx = sm_dict.nx
        ny = sm_dict.ny
        lats = np.linspace(north, south, ny)
        lons = np.linspace(west, east, nx)
        lon, lat = np.meshgrid(lons, lats)
        dep = np.zeros_like(lon)
        return cls(gmpe, lon, lat, dep, origin, rup)

    def getDistanceContext(self):
        """
        Returns:
            Openquake distance context
            `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html?highlight=distancescontext#openquake.hazardlib.gsim.base.DistancesContext>`__.
        """
        return copy.deepcopy(self._distance_context)

    def getOrigin(self):
        """
        Returns:
            Shakemap Origin object. 
        """
        return copy.deepcopy(self._origin)

    def _calcDistanceContext(self, gmpe, lat, lon, dep):
        """
        Create a DistancesContext object.

        Args:
            gmpe (GMPE): Concrete subclass of GMPE
                (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
                can be individual instance or list of instances.
            lat (array): Numpy array of latitudes.
            lon (array): Numpy array of longitudes.
            dep (array): Numpy array of depths (km).
        
        Returns:
            DistancesContext object with distance grids required by input 
            gmpe(s).
        
        Raises:
            TypeError: If gmpe is not a subclass of GMPE. 
        """
        if not isinstance(gmpe, list):
            gmpe = [gmpe]

        # require rhypo always
        requires = set(['rhypo'])

        for ig in gmpe:
            if not isinstance(ig, GMPE):
                raise TypeError(
                    'getDistanceContext() cannot work with objects of type "%s"' %
                    type(ig))
            requires = requires | ig.REQUIRES_DISTANCES

        context = base.DistancesContext()

        ddict = get_distance(list(requires), lat, lon, dep, self._origin,
                             self._rupture)

        for method in requires:
            (context.__dict__)[method] = ddict[method]

        return context

def get_distance_measures():
    """
    Returns a list of strings specifying the distance measure types 
    (e.g. "repi", "rhypo") available for the "methods" argument to 
    get_diistance().

    Returns:
        A list of strings.
    """

    return ['repi', 'rhypo', 'rjb', 'rrup', 'rx', 'ry', 'ry0', 'U', 'T']

def get_distance(methods, lat, lon, dep, origin = None, rupture = None, 
                 mesh_dx = 0.5):
    """
    Calculate distance using any one of a number of distance measures.
    One of quadlist OR hypo must be specified. The following table gives
    the allowed distance strings and a description of each. 

    +--------+----------------------------------------------------------+
    | String | Description                                              |
    +========+==========================================================+
    | repi   | Distance to epicenter.                                   |
    +--------+----------------------------------------------------------+
    | rhypo  | Distance to hypocenter.                                  |
    +--------+----------------------------------------------------------+
    | rjb    | Joyner-Boore distance; this is closest distance to the   |
    |        | surface projection of the rupture plane.                 |
    +--------+----------------------------------------------------------+
    | rrup   | Rupture distance; closest distance to the rupture plane. |
    +--------+----------------------------------------------------------+
    | rx     | Strike-normal distance; same as GC2 coordiante T.        |
    +--------+----------------------------------------------------------+
    | ry     | Strike-parallel distance; same as GC2 coordiante U, but  |
    |        | with a shift in origin definition. See Spudich and Chiou |
    |        | (2015) http://dx.doi.org/10.3133/ofr20151028.            |
    +--------+----------------------------------------------------------+
    | ry0    | Horizontal distance off the end of the rupture measured  |
    |        | parallel to strike. Can only be zero or positive. We     |
    |        | compute this as a function of GC2 coordinate U.          |
    +--------+----------------------------------------------------------+
    | U      | GC2 coordinate U.                                        |
    +--------+----------------------------------------------------------+
    | T      | GC2 coordinate T.                                        |
    +--------+----------------------------------------------------------+

    Args:
        methods (list): List of strings (or just a string) of distances to compute.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longidues.
        dep (array): A numpy array of depths (km).
        origin (Origin): A ShakeMap Origin instance.
        rupture (Rupture): A ShakeMap Rupture instance.
        mesh_dx (float): Mesh spacing in km. Only used if rupture is an EdgeRupture
            instance (meaning that distances are computed based on meshing the 
            rupture surface). 

    Returns:
       dict: dictionary of numpy arrays of distances, size of lon.shape.

    Note that magnitude-based median rupture distances are computed if a PointRupture
    is provided for the 'rupture' argument. If no rupture is provided and rupture
    distances are requested then Rhypo is substituded for Rrup and Repi for Rjb. 
    """

    # Dictionary for holding/returning the requested distances
    distdict = dict()

    # Coerce methods into list if it isn't
    if not isinstance(methods, list):
        methods = [methods]

    # Check that all requested distances are available
    methods_available = set(get_distance_measures())
    if not set(methods).issubset(methods_available):
        raise NotImplementedError(
            'One or more requested distance method is not '
            'valid or is not implemented yet')

    # Check dimensions of site coordinates
    if (lat.shape != lon.shape) or (lat.shape != dep.shape):
        raise ShakeMapException('lat, lon, and dep must have the same shape.')

    #---------------------------------------------------------------------------
    # Point distances
    #---------------------------------------------------------------------------
    if 'rhypo' in methods:
        if origin is not None:
            distdict['rhypo'] = origin.computeRhyp(lon, lat, dep)
        else:
            raise Exception('Origin required for Rhypo.')

    if 'repi' in methods:
        if origin is not None:
            distdict['repi'] = origin.computeRepi(lon, lat, dep)
        else:
            raise Exception('Origin required for Repi.')

    #---------------------------------------------------------------------------
    # Rupture distances
    #---------------------------------------------------------------------------
    gc2_distances = set(['rx', 'ry', 'ry0', 'U', 'T'])
    if rupture is not None:
        if 'rrup' in methods:
            # Question: would it be better to have all the computeR* methods
            # have the same arguments (but some unused)? Benefit would be 
            # avoiding these annoying if-statements. 
            if isinstance(rupture, QuadRupture):
                distdict['rrup'] = rupture.computeRrup(lon, lat, dep)
            elif isinstance(rupture, EdgeRupture):
                distdict['rrup'] = rupture.computeRrup(lon, lat, dep, mesh_dx)
            elif isinstance(rupture, PointRupture):
                tmp = rupture.computeRrup(lon, lat, dep, origin)
                distdict['rrup'] = tmp[0]
                distdict['rrup_var'] = tmp[1]
        if 'rjb' in methods:
            if isinstance(rupture, QuadRupture):
                distdict['rjb'] = rupture.computeRjb(lon, lat, dep)
            elif isinstance(rupture, EdgeRupture):
                distdict['rjb'] = rupture.computeRjb(lon, lat, dep, mesh_dx)
            elif isinstance(rupture, PointRupture):
                tmp = rupture.computeRjb(lon, lat, dep, origin)
                distdict['rjb'] = tmp[0]
                distdict['rjb_var'] = tmp[1]

        # If any of the GC2-related distances are requested, may as well do all
        if len(set(methods).intersection(gc2_distances)) > 0: 
            distdict.update(rupture.computeGC2(lon, lat, dep))
    else:
        if origin is not None:
            if 'rrup' in methods:
                distdict['rrup'] = origin.computeRhyp(lon, lat, dep)
            if 'rjb' in methods:
                distdict['rjb'] = origin.computeRepi(lon, lat, dep)
        else:
            raise Exception('An Origin or Rupture is required..')

        # If an origin is provided but no rupture AND a GC2 distance is 
        # requested, then a sensible thing to do is to return the GC2
        # default values in a PointRupture
        if len(set(methods).intersection(gc2_distances)) > 0:
            tmp_point_rupture = PointRupture()
            distdict.update(tmp_point_rupture.computeGC2(lon, lat, dep))
    return distdict



