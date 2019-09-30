#!/usr/bin/env python3

# stdlib modules
import itertools as it
import copy

# third party imports
import numpy as np
from openquake.hazardlib.geo.utils import OrthographicProjection
from openquake.hazardlib.geo import geodetic

from impactutils.vectorutils.vector import Vector
from shakelib.rupture.utils import reverse_quad, get_quad_length


def _computeGC2(rupture, lon, lat, depth):
    """
    Method for computing GC2 from a ShakeMap Rupture instance.

    Args:
        rupture (Rupture): ShakeMap rupture object.
        lon (array): Numpy array of site longitudes.
        lat (array): Numpy array of site latitudes.
        depth (array): Numpy array of site depths.

    Returns:
        dict: Dictionary of GC2 distances. Keys include "T", "U", "rx"
            "ry", "ry0".
    """

    quadlist = rupture.getQuadrilaterals()
    quadgc2 = copy.deepcopy(quadlist)

    oldshape = lon.shape

    if len(oldshape) == 2:
        newshape = (oldshape[0] * oldshape[1], 1)
    else:
        newshape = (oldshape[0], 1)

    # -------------------------------------------------------------------------
    # Define a projection that spans sites and rupture
    # -------------------------------------------------------------------------

    all_lat = np.append(lat, rupture.lats)
    all_lon = np.append(lon, rupture.lons)

    west = np.nanmin(all_lon)
    east = np.nanmax(all_lon)
    south = np.nanmin(all_lat)
    north = np.nanmax(all_lat)
    proj = OrthographicProjection(west, east, north, south)

    totweight = np.zeros(newshape, dtype=lon.dtype)
    GC2T = np.zeros(newshape, dtype=lon.dtype)
    GC2U = np.zeros(newshape, dtype=lon.dtype)

    # -------------------------------------------------------------------------
    # First sort out strike discordance and nominal strike prior to
    # starting the loop if there is more than one group/trace.
    # -------------------------------------------------------------------------
    group_ind = rupture._getGroupIndex()

    # Need group_ind as numpy array for sensible indexing...
    group_ind_np = np.array(group_ind)
    uind = np.unique(group_ind_np)
    n_groups = len(uind)

    if n_groups > 1:
        # ---------------------------------------------------------------------
        # The first thing we need to worry about is finding the coordinate
        # shift. U's origin is "selected from the two endpoints most
        # distant from each other."
        # ---------------------------------------------------------------------

        # Need to get index of first and last quad
        # for each segment
        iq0 = np.zeros(n_groups, dtype='int16')
        iq1 = np.zeros(n_groups, dtype='int16')
        for k in uind:
            ii = [i for i, j in enumerate(group_ind) if j == uind[k]]
            iq0[k] = int(np.min(ii))
            iq1[k] = int(np.max(ii))

        # ---------------------------------------------------------------------
        # This is an iterator for each possible combination of traces
        # including trace orientations (i.e., flipped).
        # ---------------------------------------------------------------------

        it_seg = it.product(it.combinations(uind, 2),
                            it.product([0, 1], [0, 1]))

        # Placeholder for the trace pair/orientation that gives the
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

        # ---------------------------------------------------------------------
        # A0 and A1 are the furthest two segment endpoints, but we still
        # need to sort out which one is the "origin".
        # ---------------------------------------------------------------------

        # This goofy while-loop is to adjust the side of the rupture where the
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
            e_j = np.zeros(n_groups)
            b_prime = [None] * n_groups
            for j in range(n_groups):
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
            for j in range(n_groups):
                b.x = b.x + b_prime[j].x * dc[j]
                b.y = b.y + b_prime[j].y * dc[j]
                b.z = b.z + b_prime[j].z * dc[j]
            bhat = b.norm()
            dummy = bhat.dot(ahat)
            if dummy < 0:
                tmpA0 = copy.deepcopy(A0)
                A0 = copy.deepcopy(A1)
                A1 = tmpA0

        # ---------------------------------------------------------------------
        # To fix discordancy, need to flip quads and rearrange
        # the order of quadgc2
        # ---------------------------------------------------------------------

        # 1) flip quads
        for i in range(len(quadgc2)):
            if dc[group_ind[i]] < 0:
                quadgc2[i] = reverse_quad(quadgc2[i])

        # 2) rearrange quadlist order
        qind = np.arange(len(quadgc2))
        for i in range(n_groups):
            qsel = qind[group_ind_np == uind[i]]
            if dc[i] < 0:
                qrev = qsel[::-1]
                qind[group_ind_np == uind[i]] = qrev

        quadgc2old = copy.deepcopy(quadgc2)
        for i in range(len(qind)):
            quadgc2[i] = quadgc2old[qind[i]]

        # End of if-statement for adjusting group discordancy

    s_i = 0.0
    l_i = np.zeros(len(quadgc2))

    for i in range(len(quadgc2)):
        G0, G1, G2, G3 = quadgc2[i]

        # Compute u_i and t_i for this quad
        t_i = __calc_t_i(G0, G1, lat, lon, proj)
        u_i = __calc_u_i(G0, G1, lat, lon, proj)

        # Quad length (top edge)
        l_i[i] = get_quad_length(quadgc2[i])

        # ---------------------------------------------------------------------
        # Weight of segment, three cases
        # ---------------------------------------------------------------------

        # Case 3: t_i == 0 and 0 <= u_i <= l_i
        w_i = np.zeros_like(t_i)

        # To avoid division by zero in totweight later on:
        ix = (t_i == 0) & (0 <= u_i) & (u_i <= l_i[i])
        totweight[ix] = 1.0

        # Case 1:
        ix = t_i != 0
        w_i[ix] = (1.0 / t_i[ix]) * (np.arctan(
            (l_i[i] - u_i[ix]) / t_i[ix]) - np.arctan(-u_i[ix] / t_i[ix]))

        # Case 2:
        ix = (t_i == 0) & ((u_i < 0) | (u_i > l_i[i]))
        w_i[ix] = 1 / (u_i[ix] - l_i[i]) - 1 / u_i[ix]

        totweight = totweight + w_i
        GC2T = GC2T + w_i * t_i

        if n_groups == 1:
            GC2U = GC2U + w_i * (u_i + s_i)
        else:
            if i == 0:
                qind = np.array(range(len(quadgc2)))
                l_kj = 0
                s_ij_1 = 0
            else:
                l_kj = l_i[(group_ind_np == group_ind_np[i]) & (qind < i)]
                s_ij_1 = np.sum(l_kj)

            # First endpoint in the current 'group' (or 'trace' in GC2 terms)
            p1 = Vector.fromPoint(quadgc2[iq0[group_ind[i]]][0])
            s_ij_2 = (p1 - p_origin).dot(np.sign(E) * ahat) / 1000.0

            # Above is GC2N, for GC2T use:
            # s_ij_2 = (p1 - p_origin).dot(bhat) / 1000.0

            s_ij = s_ij_1 + s_ij_2
            GC2U = GC2U + w_i * (u_i + s_ij)

        s_i = s_i + l_i[i]

    GC2T = GC2T / totweight
    GC2U = GC2U / totweight

    # Dictionary for holding the distances
    distdict = dict()

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
    if n_groups > 1:
        s_i = s_ij + l_i[-1]
    ix = GC2U > s_i
    Ry0[ix] = GC2U[ix] - s_i
    Ry0 = Ry0.reshape(oldshape)
    distdict['ry0'] = Ry0

    return distdict


def __calc_u_i(P0, P1, lat, lon, proj):
    """
    Calculate u_i distance. See Spudich and Chiou OFR 2015-1028. This is the
    distance along strike from the first vertex (P0) of the i-th segment.

    Args:
        P0 (point): OQ Point object, representing the first top-edge vertex of
            a rupture quadrilateral.
        P1 (point): OQ Point object, representing the second top-edge vertex of
            a rupture quadrilateral.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longitudes.
        proj (object): An orthographic projection from
            openquake.hazardlib.geo.utils.OrthographicProjection.

    Returns:
        array: Array of size N of distances (in km) from input points to
            rupture surface.

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

    # Dot product gives u_i
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

    Args:
        P0 (point): OQ Point object, representing the first top-edge vertex of
            a rupture quadrilateral.
        P1 (point): OQ Point object, representing the second top-edge vertex of
            a rupture quadrilateral.
        lat (array): A numpy array of latitudes.
        lon (array): A numpy array of longitudes.
        proj (object): An orthographic projection from
            openquake.hazardlib.geo.utils.OrthographicProjection.

    Returns:
        array: Array of size N of distances (in km) from input points to
            rupture surface.
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
