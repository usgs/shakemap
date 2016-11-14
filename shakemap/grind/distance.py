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
from shakemap.grind.rupture import get_quad_length
from shakemap.grind.rupture import reverse_quad
from shakemap.grind.source import rake_to_mech

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

    def __init__(self, gmpe, source, lat, lon, dep, use_median_distance=True):
        """
        :param gmpe:
            Concrete subclass of GMPE
            `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
            can be individual instance or list of instances.
        :param source:
            Shakemap Source instance. For point-source distances (Repi, Rhyp)
            the hypocenter is taken from the source instance. The finite-rupture
            distances (Rrup, Rjb, ...), are based on the rupture representation
            in the source instance if available. If a Rupture instance is
            unavailalbe, the finite-rupture distances are computed from 
            point-source distances, in which case the earthquake manitude from
            the Source instance is required for median distance corrections if
            use_median_distance is True.
        :param lat:
            A numpy array of site latitudes.
        :param lon:
            A numpy array of site longitudes.
        :param dep:
            A numpy array of site depths (km); down is positive. 
        :param use_median_distance:
            Boolean; only used if GMPE requests rupture distances and not
            rupture is availalbe. Default is True, meaning that point-source
            distances are adjusted based on magnitude to get the median rupture
            distance.
        :returns:
            Distance object.
        """
        self._source = source

        self._distance_context = self._calcDistanceContext(
            gmpe, lat, lon, dep, use_median_distance)

    @classmethod
    def fromSites(cls, gmpe, source, sites, use_median_distance=True):
        """
        Convenience class method to construct a Distance object from a sites
        object.

        :param gmpe:
            Concrete subclass of GMPE
            `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__;
            can be individual instance or list of instances.
        :param source:
            Shakeamp Source object.
        :param sites:
            Shakemap Sites object.
        :param use_median_distance:
            Boolean; only used if GMPE requests rupture distances and not rupture
            is availalbe. Default is True, meaning that point-source distances
            are adjusted based on magnitude to get the median rupture distance.
        :returns:
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
        return cls(gmpe, source, lat, lon, dep, use_median_distance)

    def getDistanceContext(self):
        """
        :returns:
            Openquake distance context
            `[link] <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html?highlight=distancescontext#openquake.hazardlib.gsim.base.DistancesContext>`__.
        """
        return copy.deepcopy(self._distance_context)

    def getSource(self):
        """
        :returns:
            Shakemap Source object. 
        """
        return copy.deepcopy(self._source)

    def _calcDistanceContext(self, gmpe, lat, lon, dep,
                             use_median_distance=True):
        """
        Create a DistancesContext object.

        :param gmpe:
            Concrete subclass of GMPE
            (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
            can be individual instance or list of instances.
        :param lat:
            Numpy array of latitudes.
        :param lon:
            Numpy array of longitudes.
        :param dep:
            Numpy array of depths (km).
        :param use_median_distance:
            Boolean; only used if GMPE requests rupture distances and not rupture
            is availalbe. Default is True, meaning that point-source distances
            are adjusted based on magnitude to get the median rupture distance.
        :returns:
            DistancesContext object with distance grids required by input 
            gmpe(s).
        :raises TypeError:
            If gmpe is not a subclass of GMPE. 
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

        ddict = get_distance(list(requires), lat, lon, dep, self._source,
                             use_median_distance=use_median_distance)

        for method in requires:
            (context.__dict__)[method] = ddict[method]

        return context

def get_distance_measures():
    """
    Returns a list of strings specifying the distance measure types 
    (e.g. "repi", "rhypo") available for the "methods" argument to 
    get_diistance().

    :returns:
        A list of strings.
    """

    return ['repi', 'rhypo', 'rjb', 'rrup', 'rx', 'ry', 'ry0', 'U', 'T']

def get_distance(methods, lat, lon, dep, source,
                 use_median_distance=True):
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

    :param methods:
        List of strings (or just a string) of distances to compute.
    :param lat:
       A numpy array of latitudes.
    :param lon:
       A numpy array of longidues.
    :param dep:
       A numpy array of depths (km).
    :param source:
       source instance.
    :param use_median_distance:
        Boolean; only used if GMPE requests rupture distances and not rupture is
        availalbe. Default is True, meaning that point-source distances are
        adjusted based on magnitude to get the median rupture distance.
    :returns:
       dictionary of numpy arrays of distances, size of lon.shape
       IMPORTANT: If a finite rupture is not supplied, and the distance
       measures requested include rx, ry, ry0, U, or T, then zeros will
       be returned; if rjb is requested, repi will be returned; if rrup
       is requested, rhypo will be returned.
    """
    rupture = source.getRupture()
    hypo = source.getHypo()
    if rupture is not None:
        quadlist = rupture.getQuadrilaterals()
        # Need a copy for GC2 since order of verticies/quads needs to be modivied.
        quadgc2 = copy.deepcopy(quadlist)
    else:
        quadlist = None

    # Dictionary for holding the distances
    distdict = dict()

    if not isinstance(methods, list):
        methods = [methods]

    methods_available = set(get_distance_measures())
    if not set(methods).issubset(methods_available):
        raise NotImplementedError(
            'One or more requested distance method is not '
            'valid or is not implemented yet')

    if (lat.shape == lon.shape) and (lat.shape == dep.shape):
        pass
    else:
        raise ShakeMapException('lat, lon, and dep must have the same shape.')

    oldshape = lon.shape

    if len(oldshape) == 2:
        newshape = (oldshape[0] * oldshape[1], 1)
    else:
        newshape = (oldshape[0], 1)

    if ('rrup' in methods) or ('rjb' in methods):
        x, y, z = latlon2ecef(lat, lon, dep)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

    # Define a projection that spands sites and rupture
    if rupture is None:
        all_lat = lat
        all_lon = lon
    else:
        all_lat = np.append(lat, rupture.getLats())
        all_lon = np.append(lon, rupture.getLons())

    west = np.nanmin(all_lon)
    east = np.nanmax(all_lon)
    south = np.nanmin(all_lat)
    north = np.nanmax(all_lat)
    proj = get_orthographic_projection(west, east, north, south)

    # ---------------------------------------------
    # Distances that do not require loop over quads
    # ---------------------------------------------

    if ('repi' in methods) or \
       (('rjb' in methods) and (quadlist is None)) or \
       (('rrup' in methods) and (quadlist is None)) or \
       (('ry0' in methods) and (quadlist is None)) or \
       (('rx' in methods) and (quadlist is None)) or \
       (('T' in methods) and (quadlist is None)) or \
       (('U' in methods) and (quadlist is None)):
        # I don't think this error check makes sense any more because hypo
        # is assigned above with source.getHypo() that constructs it from
        # source._event_dict entries. 
        if hypo is None:
            raise ShakeMapException('Cannot calculate epicentral distance '
                                    'without a point object')
        repidist = geodetic.distance(hypo.longitude, hypo.latitude, 0.0,
                                     lon, lat, dep)
        repidist = repidist.reshape(oldshape)
        distdict['repi'] = repidist

    if ('rhypo' in methods) or \
       (('rrup' in methods) and (quadlist is None)):
        if hypo is None:
            raise ShakeMapException('Cannot calculate epicentral distance '
                                    'without a point object')
        rhypodist = geodetic.distance(
            hypo.longitude, hypo.latitude, hypo.depth, lon, lat, dep)
        rhypodist = rhypodist.reshape(oldshape)
        distdict['rhypo'] = rhypodist

    # --------------------------------------------------------
    # Loop over quadlist for those distances that require loop
    # --------------------------------------------------------
    if 'rrup' in methods:
        minrrup = np.ones(newshape, dtype=lon.dtype) * 1e16
    if 'rjb' in methods:
        minrjb = np.ones(newshape, dtype=lon.dtype) * 1e16
    if ('rx' in methods) or ('ry' in methods) or \
       ('ry0' in methods) or ('U' in methods) or ('T' in methods):
        totweight = np.zeros(newshape, dtype=lon.dtype)
        GC2T = np.zeros(newshape, dtype=lon.dtype)
        GC2U = np.zeros(newshape, dtype=lon.dtype)

        if quadlist is not None:
            #-----------------------------------------------------------------
            # For these distances, we need to sort out strike discordance and
            # nominal strike prior to starting the loop if there is more than
            # one segment.
            #-----------------------------------------------------------------

            segind = rupture._getSegmentIndex()
            segindnp = np.array(segind)
            uind = np.unique(segind)
            nseg = len(uind)

            #-------------------------------------------------------------------
            # The first thing we need to worry about is finding the coordinate
            # shift. U's origin is " selected from the two endpoints most
            # distant from each other." 
            #-------------------------------------------------------------------

            if nseg > 1:
                # Need to get index of first and last quad
                # for each segment
                iq0 = np.zeros(nseg, dtype='int16')
                iq1 = np.zeros(nseg, dtype='int16')
                for k in uind:
                    ii = [i for i, j in enumerate(segind) if j == uind[k]]
                    iq0[k] = int(np.min(ii))
                    iq1[k] = int(np.max(ii))

                #---------------------------------------------------------------
                # This is an iterator for each possible combination of segments
                # including segment orientations (i.e., flipped). 
                #---------------------------------------------------------------

                it_seg = it.product(it.combinations(uind, 2),
                                    it.product([0, 1], [0, 1]))

                # Placeholder for the segment pair/orientation that gives the
                # largest distance. 
                dist_save = 0

                for k in it_seg:
                    s0ind = k[0][0]
                    s1ind = k[0][1]
                    p0ind = k[1][0]
                    p1ind = k[1][1]
                    if p0ind == 0:
                        P0 = quadlist[iq0[s0ind]][0]
                    else:
                        P0 = quadlist[iq1[s0ind]][1]
                    if p1ind == 0:
                        P1 = quadlist[iq1[s1ind]][0]
                    else:
                        P1 = quadlist[iq0[s1ind]][1]

                    dist = geodetic.distance(P0.longitude, P0.latitude, 0.0,
                                             P1.longitude, P1.latitude, 0.0)
                    if dist > dist_save:
                        dist_save = dist
                        A0 = P0
                        A1 = P1

                #---------------------------------------------------------------
                # A0 and A1 are the furthest two segment endpoints, but we still
                # need to sort out which one is the "origin".
                #---------------------------------------------------------------

                # Goofy while-loop is to adjust the side of the rupture where the
                # origin is located
                dummy = -1
                while dummy < 0:
                    A0.depth = 0
                    A1.depth = 0
                    p_origin = Vector.fromPoint(A0)
                    a0 = Vector.fromPoint(A0)
                    a1 = Vector.fromPoint(A1)
                    ahat = (a1 - a0).norm()

                    # Loop over traces
                    e_j = np.zeros(nseg)
                    b_prime = [None] * nseg
                    for j in range(nseg):
                        P0 = quadlist[iq0[j]][0]
                        P1 = quadlist[iq1[j]][1]
                        P0.depth = 0
                        P1.depth = 0
                        p0 = Vector.fromPoint(P0)
                        p1 = Vector.fromPoint(P1)
                        b_prime[j] = p1 - p0
                        e_j[j] = ahat.dot(b_prime[j])
                    E = np.sum(e_j)

                    # List of discordancy
                    dc = [np.sign(a) * np.sign(E) for a in e_j]
                    b = Vector(0, 0, 0)
                    for j in range(nseg):
                        b.x = b.x + b_prime[j].x * dc[j]
                        b.y = b.y + b_prime[j].y * dc[j]
                        b.z = b.z + b_prime[j].z * dc[j]
                    bhat = b.norm()
                    dummy = bhat.dot(ahat)
                    if dummy < 0:
                        tmpA0 = copy.copy(A0)
                        tmpA1 = copy.copy(A1)
                        A0 = tmpA1
                        A1 = tmpA0

                # To fix discordancy, need to flip quads and rearrange
                # the order of quadgc2

                # 1) flip quads
                for i in range(len(quadgc2)):
                    if dc[segind[i]] < 0: #***************U*UUUUUUUSDFUSDfkjjhakjsdhfljkhn
                        quadgc2[i] = reverse_quad(quadgc2[i])

                # 2) rearrange quadlist to remove discordancy
                qind = np.arange(len(quadgc2))
                segnp = np.array(segind)
                for i in range(nseg):
                    qsel = qind[segnp == uind[i]]
                    if dc[i] < 0:
                        qrev = qsel[::-1]
                        qind[segnp == uind[i]] = qrev

                quadgc2old = copy.deepcopy(quadgc2)
                for i in range(len(qind)):
                    quadgc2[i] = quadgc2old[qind[i]]


    if quadlist is not None:
        # Length of prior segments
        s_i = 0.0
        l_i = np.zeros(len(quadlist))
        for i in range(len(quadlist)):
            P0, P1, P2, P3 = quadlist[i]
            G0, G1, G2, G3 = quadgc2[i]

            if 'rrup' in methods:
                rrupdist = _calc_rupture_distance(P0, P1, P2, P3, sites_ecef)
                minrrup = np.minimum(minrrup, rrupdist)

            if 'rjb' in methods:
                S0 = copy.deepcopy(P0)
                S1 = copy.deepcopy(P1)
                S2 = copy.deepcopy(P2)
                S3 = copy.deepcopy(P3)
                S0.depth = 0.0
                S1.depth = 0.0
                S2.depth = 0.0
                S3.depth = 0.0
                rjbdist = _calc_rupture_distance(S0, S1, S2, S3, sites_ecef)
                minrjb = np.minimum(minrjb, rjbdist)

            if ('rx' in methods) or ('ry' in methods) or \
               ('ry0' in methods) or ('U' in methods) or ('T' in methods):
                # Rx, Ry, and Ry0 are all computed if one is requested since
                # they all require similar information for the weights. This
                # isn't necessary for a single segment rupture though.
                # Note, we are basing these calculations on GC2 coordinates U
                # and T as described in:
                # Spudich and Chiou (2015)
                # http://dx.doi.org/10.3133/ofr20151028.

                # Compute u_i and t_i for this quad
                t_i = __calc_t_i(G0, G1, lat, lon, proj)
                u_i = __calc_u_i(G0, G1, lat, lon, proj)

                # Quad length
                l_i[i] = get_quad_length(quadlist[i])

                # Weight of segment, three cases
                # Case 3: t_i == 0 and 0 <= u_i <= l_i
                w_i = np.zeros_like(t_i)
                # Case 1:
                ix = t_i != 0
                w_i[ix] = (1.0 / t_i[ix]) * (np.arctan((l_i[i] -
                    u_i[ix]) / t_i[ix]) - np.arctan(-u_i[ix] / t_i[ix]))
                # Case 2:
                ix = (t_i == 0) & ((u_i < 0) | (u_i > l_i[i]))
                w_i[ix] = 1 / (u_i[ix] - l_i[i]) - 1 / u_i[ix]

                totweight = totweight + w_i
                GC2T = GC2T + w_i * t_i
                if nseg == 1:
                    GC2U = GC2U + w_i * (u_i + s_i)
                else:
                    if i == 0:
                        qind = np.array(range(len(quadlist)))
                        l_kj = 0
                        s_ij_1 = 0
                    else:
                        l_kj = l_i[(segindnp == segindnp[i]) & (qind < i)]
                        s_ij_1 = np.sum(l_kj)

                    p1 = Vector.fromPoint(quadgc2[iq0[segind[i]]][0])
                    s_ij_2 = (p1 - p_origin).dot(np.sign(E)*ahat) / 1000.0
                    # Above is GC2N, for GC2T use:
#                    s_ij_2 = (p1 - p_origin).dot(bhat) / 1000.0


                    s_ij = s_ij_1 + s_ij_2
                    GC2U = GC2U + w_i * (u_i + s_ij)
                s_i = s_i + l_i[i]

        # Collect distances from loop into the distance dict
        if 'rjb' in methods:
            minrjb = minrjb.reshape(oldshape)
            distdict['rjb'] = minrjb

        if ('rx' in methods) or ('ry' in methods) or \
           ('ry0' in methods) or ('U' in methods) or ('T' in methods):
            # Normalize by sum of quad weights
            GC2T = GC2T / totweight
            GC2U = GC2U / totweight
            distdict['T'] = copy.deepcopy(GC2T).reshape(oldshape)
            distdict['U'] = copy.deepcopy(GC2U).reshape(oldshape)

            # Take care of Rx
            Rx = copy.deepcopy(GC2T)  # preserve sign (no absolute value)
            Rx = Rx.reshape(oldshape)
            distdict['rx'] = Rx

            # Ry
            Ry = GC2U - s_i / 2.0
            Ry = Ry.reshape(oldshape)
            distdict['ry'] = Ry

            # Ry0
            Ry0 = np.zeros_like(GC2U)
            ix = GC2U < 0
            Ry0[ix] = np.abs(GC2U[ix])
            if nseg > 1:
                s_i = s_ij + l_i[-1]
            ix = GC2U > s_i
            Ry0[ix] = GC2U[ix] - s_i
            Ry0 = Ry0.reshape(oldshape)
            distdict['ry0'] = Ry0

        if 'rrup' in methods:
            minrrup = minrrup.reshape(oldshape)
            distdict['rrup'] = minrrup

    else:
        if 'rjb' in methods:
            if use_median_distance:
                warnings.warn(
                    'No rupture; Replacing rjb with median rjb given M and repi.')
                cdir, tmp = os.path.split(__file__)

                # -------------------
                # Sort out file names
                # -------------------
                mech = source.getEventParam('mech')
                if not hasattr(source, '_tectonic_region'):
                    rf = os.path.join(
                        cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
                    vf = os.path.join(
                        cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")
                elif source._tectonic_region == 'Active Shallow Crust':
                    if mech == 'ALL':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p7_seis0_20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p7_seis0_20_Var.csv")
                    elif mech == 'RS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechR_ar1p7_seis0_20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechR_ar1p7_seis0_20_Var.csv")
                    elif mech == 'NM':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechN_ar1p7_seis0_20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechN_ar1p7_seis0_20_Var.csv")
                    elif mech == 'SS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechSS_ar1p7_seis0_20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_WC94_mechSS_ar1p7_seis0_20_Var.csv")
                elif source._tectonic_region == 'Stable Shallow Crust':
                    if mech == 'ALL':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechA_ar1p0_seis0_15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechA_ar1p0_seis0_15_Var.csv")
                    elif mech == 'RS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechR_ar1p0_seis0_15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechR_ar1p0_seis0_15_Var.csv")
                    elif mech == 'NM':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechN_ar1p0_seis0_15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechN_ar1p0_seis0_15_Var.csv")
                    elif mech == 'SS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechSS_ar1p0_seis0_15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rjb_S14_mechSS_ar1p0_seis0_15_Var.csv")
                else:
                    warnings.warn(
                        'Unsupported tectonic region; using coefficients for unknown'
                        'tectonic region.')
                    rf = os.path.join(
                        cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p0_seis0_20_Ratios.csv")
                    vf = os.path.join(
                        cdir, "data", "ps2ff", "Rjb_WC94_mechA_ar1p0_seis0_20_Var.csv")

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
                repis = distdict['repi']
                mags = np.ones_like(repis) * source.getEventParam('mag')
                rjb_hat = repi2rjb_tbl(repis, mags)
                distdict['rjb'] = rjb_hat
                # -------------------
                # Additional Variance
                # -------------------
                repi2rjbvar_ratios_tbl = pd.read_csv(vf, comment='#')
                repi2rjbvar_grid = repi2rjbvar_ratios_tbl.values[:, 1:]
                repi2rjbvar_obj = spint.RectBivariateSpline(
                    dist_list, mag_list, repi2rjbvar_grid, kx=1, ky=1)
                rjbvar = repi2rjbvar_obj.ev(np.log(repis), mags)
                distdict['rjbvar'] = rjbvar
            else:
                warnings.warn('No rupture; Replacing rjb with repi')
                distdict['rjb'] = distdict['repi'].copy()
        if 'rrup' in methods:
            if use_median_distance:
                warnings.warn(
                    'No rupture; Replacing rrup with median rrup given M and repi.')
                cdir, tmp = os.path.split(__file__)

                # -------------------
                # Sort out file names
                # -------------------
                rake = source._event_dict.get('rake')
                mech = rake_to_mech(rake)
                if not hasattr(source, '_tectonic_region'):
                    rf = os.path.join(
                        cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
                    vf = os.path.join(
                        cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")
                elif source._tectonic_region == 'Active Shallow Crust':
                    if mech == 'ALL':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p7_seis0-20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p7_seis0-20_Var.csv")
                    elif mech == 'RS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechR_ar1p7_seis0-20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechR_ar1p7_seis0-20_Var.csv")
                    elif mech == 'NM':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechN_ar1p7_seis0-20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechN_ar1p7_seis0-20_Var.csv")
                    elif mech == 'SS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechSS_ar1p7_seis0-20_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_WC94_mechSS_ar1p7_seis0-20_Var.csv")
                elif source._tectonic_region == 'Stable Shallow Crust':
                    if mech == 'ALL':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechA_ar1p0_seis0-15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechA_ar1p0_seis0-15_Var.csv")
                    elif mech == 'RS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechR_ar1p0_seis0-15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechR_ar1p0_seis0-15_Var.csv")
                    elif mech == 'NM':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechN_ar1p0_seis0-15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechN_ar1p0_seis0-15_Var.csv")
                    elif mech == 'SS':
                        rf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechSS_ar1p0_seis0-15_Ratios.csv")
                        vf = os.path.join(
                            cdir, "data", "ps2ff", "Rrup_S14_mechSS_ar1p0_seis0-15_Var.csv")
                else:
                    warnings.warn(
                        'Unsupported tectonic region; using coefficients for unknown'
                        'tectonic region.')
                    rf = os.path.join(
                        cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p0_seis0-20_Ratios.csv")
                    vf = os.path.join(
                        cdir, "data", "ps2ff", "Rrup_WC94_mechA_ar1p0_seis0-20_Var.csv")

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
                repis = distdict['repi']
                mags = np.ones_like(repis) * source.getEventParam('mag')
                rrup_hat = repi2rrup_tbl(repis, mags)
                distdict['rrup'] = rrup_hat

                # -------------------
                # Additional Variance
                # -------------------
                repi2rrupvar_ratios_tbl = pd.read_csv(vf, comment='#')
                repi2rrupvar_grid = repi2rrupvar_ratios_tbl.values[:, 1:]
                repi2rrupvar_obj = spint.RectBivariateSpline(
                    dist_list, mag_list, repi2rrupvar_grid, kx=1, ky=1)
                rrupvar = repi2rrupvar_obj.ev(np.log(repis), mags)
                distdict['rrupvar'] = rrupvar
            else:
                warnings.warn('No rupture; Replacing rrup with rhypo')
                distdict['rrup'] = distdict['rhypo'].copy()
        if 'rx' in methods:
            warnings.warn('No rupture; Setting Rx to zero.')
            distdict['rx'] = np.zeros_like(distdict['repi'])
        if 'ry0' in methods:
            warnings.warn('No rupture; Setting ry0 to zero')
            distdict['ry0'] = np.zeros_like(distdict['repi'])
        if 'ry' in methods:
            warnings.warn('No rupture; Setting ry to zero')
            distdict['ry'] = np.zeros_like(distdict['repi'])
        if 'U' in methods:
            warnings.warn('No rupture; Setting U to zero')
            distdict['U'] = np.zeros_like(distdict['repi'])
        if 'T' in methods:
            warnings.warn('No rupture; Setting T to zero')
            distdict['T'] = np.zeros_like(distdict['repi'])

    return distdict


def _distance_sq_to_segment(p0, p1):
    """
    Calculate the distance^2 from the origin to a segment defined by two vectors

    :param p0: 
        Numpy array Nx3 (ECEF).
    :param p1: 
        Numpy array Nx3 (ECEF).
    :returns:
        The squared distance from the origin to a segment.
    """
    # /*
    #  * This algorithm is from (Vince's) CS1 class.
    #  * It returns the distance^2 from the origin to a segment defined
    #  * by two vectors
    #  */

    dist = np.zeros_like(p1[:, 0])
    # /* Are the two points equal? */
    idx_equal = (p0[:, 0] == p1[:, 0]) & (
        p0[:, 1] == p1[:, 1]) & (p0[:, 2] == p1[:, 2])
    dist[idx_equal] = np.sqrt(p0[idx_equal, 0]**2 +
                              p0[idx_equal, 1]**2 + p0[idx_equal, 2]**2)

    v = p1 - p0

    # /*
    #  * C1 = c1/|v| is the projection of the origin O on line (P0,P1).
    #  * If C1 is negative, then O is outside the segment and
    #  * closer to the P0 side.
    #  * If C1 is positive and >V then O is on the other side.
    #  * If C1 is positive and <V then O is inside.
    #  */

    c1 = -1 * np.sum(p0 * v, axis=1)
    idx_neg = c1 <= 0
    dist[idx_neg] = p0[idx_neg, 0]**2 + p0[idx_neg, 1]**2 + p0[idx_neg, 2]**2

    c2 = np.sum(v * v, axis=1)
    idx_less_c1 = c2 <= c1
    dist[idx_less_c1] = p1[idx_less_c1, 0]**2 + \
        p1[idx_less_c1, 1]**2 + p1[idx_less_c1, 2]**2

    idx_other = np.logical_not(idx_neg | idx_equal | idx_less_c1)

    nr, nc = p0.shape
    t1 = c1 / c2
    t1.shape = (nr, 1)
    tmp = p0 + (v * t1)
    dist[idx_other] = tmp[idx_other, 0]**2 + \
        tmp[idx_other, 1]**2 + tmp[idx_other, 2]**2

    return dist

# call this once per quad


def _calc_rupture_distance(P0, P1, P2, P3, points):
    """
    Calculate the shortest distance from a set of points to a rupture surface.

    :param P0:
        Point object, representing the first top-edge vertex of a rupture
        quadrilateral.
    :param P1:
        Point object, representing the second top-edge vertex of a rupture
        quadrilateral.
    :param P2:
        Point object, representing the first bottom-edge vertex of a rupture
        quadrilateral.
    :param P3:
        Point object, representing the second bottom-edge vertex of a rupture
        quadrilateral.
    :param points:
        Numpy array Nx3 of points (ECEF) to calculate distance from.
    :returns:
        Array of size N of distances (in km) from input points to rupture
        surface.
    """
    # Convert to ecef
    p0 = Vector.fromPoint(P0)
    p1 = Vector.fromPoint(P1)
    p2 = Vector.fromPoint(P2)
    p3 = Vector.fromPoint(P3)

    # Make a unit vector normal to the plane
    normalVector = (p1 - p0).cross(p2 - p0).norm()

    dist = np.ones_like(points[:, 0]) * np.nan

    p0d = p0.getArray() - points
    p1d = p1.getArray() - points
    p2d = p2.getArray() - points
    p3d = p3.getArray() - points

    # Create 4 planes with normals pointing outside rectangle
    n0 = (p1 - p0).cross(normalVector).getArray()
    n1 = (p2 - p1).cross(normalVector).getArray()
    n2 = (p3 - p2).cross(normalVector).getArray()
    n3 = (p0 - p3).cross(normalVector).getArray()

    sgn0 = np.signbit(np.sum(n0 * p0d, axis=1))
    sgn1 = np.signbit(np.sum(n1 * p1d, axis=1))
    sgn2 = np.signbit(np.sum(n2 * p2d, axis=1))
    sgn3 = np.signbit(np.sum(n3 * p3d, axis=1))

    inside_idx = (sgn0 == sgn1) & (sgn1 == sgn2) & (sgn2 == sgn3)
    dist[inside_idx] = np.power(np.abs(
        np.sum(p0d[inside_idx, :] * normalVector.getArray(), axis=1)), 2)

    outside_idx = np.logical_not(inside_idx)
    s0 = _distance_sq_to_segment(p0d, p1d)
    s1 = _distance_sq_to_segment(p1d, p2d)
    s2 = _distance_sq_to_segment(p2d, p3d)
    s3 = _distance_sq_to_segment(p3d, p0d)

    smin = np.minimum(np.minimum(s0, s1), np.minimum(s2, s3))
    dist[outside_idx] = smin[outside_idx]
    dist = np.sqrt(dist) / 1000.0
    shp = dist.shape
    if len(shp) == 1:
        dist.shape = (shp[0], 1)
    if np.any(np.isnan(dist)):
        raise ShakeMapException("Could not calculate some distances!")
    dist = np.fliplr(dist)
    return dist


def __calc_u_i(P0, P1, lat, lon, proj):
    """
    Calculate u_i distance. See Spudich and Chiou OFR 2015-1028. This is the 
    distance along strike from the first vertex (P0) of the i-th segment.

    :param P0:
        Point object, representing the first top-edge vertex of a rupture 
        quadrilateral.
    :param P1:
        Point object, representing the second top-edge vertex of a rupture 
        quadrilateral.
    :param lat:
        A numpy array of latitude.
    :param lon:
        A numpy array of longitude.
    :param proj:
        An orthographic projection from 
        openquake.hazardlib.geo.utils.get_orthographic_projection. 
    :returns:
        Array of size lat.shape of distances (in km).
    """

    # projected coordinates are in km
    p0x, p0y = proj(P0.x, P0.y)
    p1x, p1y = proj(P1.x, P1.y)

    # Unit vector pointing along strike
    u_i_hat = Vector(p1x - p0x, p1y - p0y, 0).norm()

    # Convert sites to Cartesian
    sx, sy = proj(lon, lat)
    sx1d = np.reshape(sx, (-1,))
    sy1d = np.reshape(sy, (-1,))

    # Vectors from P0 to sites
    r = np.zeros([len(sx1d), 2])
    r[:, 0] = sx1d - p0x
    r[:, 1] = sy1d - p0y

    # Dot product gives t_i
    u_i = np.sum(u_i_hat.getArray()[0:2] * r, axis=1)
    shp = u_i.shape
    if len(shp) == 1:
        u_i.shape = (shp[0], 1)
    u_i = np.fliplr(u_i)

    return u_i


def __calc_t_i(P0, P1, lat, lon, proj):
    """
    Calculate t_i distance. See Spudich and Chiou OFR 2015-1028. This is the
    distance measured normal to strike from the i-th segment. Values on the
    hanging-wall are positive and those on the foot-wall are negative.

    :param P0:
        Point object, representing the first top-edge vertex of a rupture
        quadrilateral.
    :param P1:
        Point object, representing the second top-edge vertex of a rupture
        quadrilateral.
    :param lat:
        A numpy array of latitudes.
    :param lon:
        A numpy array of longitudes.
    :param proj:
        An orthographic projection from 
        openquake.hazardlib.geo.utils.get_orthographic_projection. 
    :returns:
        Array of size N of distances (in km) from input points to rupture
        surface.
    """

    # projected coordinates are in km
    p0x, p0y = proj(P0.x, P0.y)
    p1x, p1y = proj(P1.x, P1.y)

    # Unit vector pointing normal to strike
    t_i_hat = Vector(p1y - p0y, -(p1x - p0x), 0).norm()

    # Convert sites to Cartesian
    sx, sy = proj(lon, lat)
    sx1d = np.reshape(sx, (-1,))
    sy1d = np.reshape(sy, (-1,))

    # Vectors from P0 to sites
    r = np.zeros([len(sx1d), 2])
    r[:, 0] = sx1d - p0x
    r[:, 1] = sy1d - p0y

    # Dot product gives t_i
    t_i = np.sum(t_i_hat.getArray()[0:2] * r, axis=1)
    shp = t_i.shape
    if len(shp) == 1:
        t_i.shape = (shp[0], 1)
    t_i = np.fliplr(t_i)
    return t_i
