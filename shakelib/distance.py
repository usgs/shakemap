#!/usr/bin/env python

# stdlib imports
import copy

# third party imports
import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.gsim import base

# local imports
from shakelib.utils.exception import ShakeLibException
from shakelib.rupture.edge_rupture import EdgeRupture


class Distance(object):
    """
    Class for distance calculations. Primary method is 'get_distance'. To
    handle multiple segment ruptures, many of the distances are based on the
    Spudich and Chiou (2015) GC2 coordinate system.


    References:
        Spudich, Paul and Chiou, Brian, 2015, Strike-parallel and strike-normal
        coordinate system around geometrically complicated rupture traces—Use
        by NGA-West2 and further improvements: U.S. Geological Survey Open-File
        Report 2015-1028, 20 p., http://dx.doi.org/10.3133/ofr20151028.
    """

    def __init__(self, gmpe, lon, lat, dep, rupture=None):
        """
        Constructor for Distance class.

        Args:
            gmpe (GMPE): Concrete subclass of GMPE
                `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
                can be individual instance or list of instances.
            lon (array): A numpy array of site longitudes.
            lat (array): A numpy array of site latitudes.
            dep (array): A numpy array of site depths (km); down is positive.
            rupture (Rupture): A Shakemap Rupture instance.

        Returns:
            Distance object.
        """  # noqa

        self._rupture = rupture

        self._distance_context = self._calcDistanceContext(
            gmpe, lat, lon, dep)

    @classmethod
    def fromSites(cls, gmpe, sites, rup):
        """
        Convenience class method to construct a Distance object from a sites
        object.

        Args:
            gmpe (GMPE): Concrete subclass of GMPE
                `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
                can be individual instance or list of instances.
            sites (Sites): A Shakemap Sites object.
            rup (Rupture): A Shakemap Rupture object.

        Returns:
            Distance object.
        """  # noqa
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
        return cls(gmpe, lon, lat, dep, rup)

    def getDistanceContext(self):
        """
        Returns:
            Openquake distance context
            `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html?highlight=distancescontext#openquake.hazardlib.gsim.base.DistancesContext>`__.
        """  # noqa
        return copy.deepcopy(self._distance_context)

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
        requires = set(['rhypo','rrup','rjb'])

        for ig in gmpe:
            if not isinstance(ig, GMPE):
                raise TypeError(
                    'getDistanceContext() cannot work with objects of '
                    'type "%s"' % type(ig))
            requires = requires | ig.REQUIRES_DISTANCES

        context = base.DistancesContext()

        if isinstance(self._rupture, EdgeRupture):
            ddict = get_distance(list(requires), lat, lon, dep, self._rupture,
                                 dx=self._rupture._mesh_dx)
        else:
            ddict = get_distance(list(requires), lat, lon, dep, self._rupture)
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


def get_distance(methods, lat, lon, dep, rupture, dx=0.5):
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
        methods (list): List of strings (or just a string) of distances to
            compute.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longidues.
        dep (array): A numpy array of depths (km).
        rupture (Rupture): A ShakeMap Rupture instance.
        dx (float): Mesh spacing for rupture; only used if rupture is an
            EdgeRupture subclass.

    Returns:
       dict: dictionary of numpy arrays of distances, size of lon.shape.
    """
    rupture._mesh_dx = dx

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
        raise ShakeLibException('lat, lon, and dep must have the same shape.')

    # -------------------------------------------------------------------------
    # Point distances
    # -------------------------------------------------------------------------
    if 'rhypo' in methods:
        distdict['rhypo'] = rupture.computeRhyp(lon, lat, dep)

    if 'repi' in methods:
        distdict['repi'] = rupture.computeRepi(lon, lat, dep)

    # -------------------------------------------------------------------------
    # Rupture distances
    # -------------------------------------------------------------------------
    gc2_distances = set(['rx', 'ry', 'ry0', 'U', 'T'])
    if 'rrup' in methods:
        distdict['rrup'] = rupture.computeRrup(lon, lat, dep)

    if 'rjb' in methods:
        distdict['rjb'] = rupture.computeRjb(lon, lat, dep)

    # If any of the GC2-related distances are requested, may as well do all
    if len(set(methods).intersection(gc2_distances)) > 0:
        distdict.update(rupture.computeGC2(lon, lat, dep))

    return distdict
