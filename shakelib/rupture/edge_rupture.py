#!/usr/bin/env python

# stdlib modules
import copy

# third party imports
import numpy as np
from openquake.hazardlib.geo.point import Point

from impactutils.vectorutils.ecef import latlon2ecef
from impactutils.vectorutils.vector import Vector
from impactutils.time.ancient_time import HistoricTime
from shakelib.rupture.base import Rupture
import shakelib.rupture.utils as utils
import shakelib.rupture.gc2 as gc2


class EdgeRupture(Rupture):
    """
    Rupture class that representst the rupture surface by specifying the top
    edge and the bottom edges. These edges do not need to be horizontal. The
    freedom to allow for non-horizontal edges (as compared to QuadRupture)
    comes at the cost of slower distance calculations. This is because the
    rupture must be discretized and then the distances are compued in a brute
    force fashion based on this mesh, which can be quite large.
    """

    def __init__(self, d, origin, mesh_dx=0.5):
        """
        Initialization of an EdgeRupture from a GeoJSON dictionary and an
        Origin.

        Args:
            d (dict): Rupture GeoJSON dictionary.
            origin (Origin): Reference to a ShakeMap Origin object.
            mesh_dx (float): Target spacing (in km) for rupture discretization;
                default is 0.5 km and it is only used if the rupture file is an
                EdgeRupture.

        Returns:
            EdgeRupture instance.

        """
        polys = d['features'][0]['geometry']['coordinates'][0]
        n_polygons = len(polys)
        toplons = np.empty(shape=(0, 0))
        toplats = np.empty(shape=(0, 0))
        topdeps = np.empty(shape=(0, 0))
        botlons = np.empty(shape=(0, 0))
        botlats = np.empty(shape=(0, 0))
        botdeps = np.empty(shape=(0, 0))
        g_ind = 0
        group_index = []
        for i in range(n_polygons):
            p = polys[i]
            n_points = len(p)
            n_pairs = int((n_points - 1) / 2)

            p_lons = [pt[0] for pt in p][0:-1]
            p_lats = [pt[1] for pt in p][0:-1]
            p_depths = [pt[2] for pt in p][0:-1]

            tlon = np.array(p_lons[0:n_pairs])
            blon = np.array(p_lons[(n_pairs):])[::-1]
            tlat = np.array(p_lats[0:n_pairs])
            blat = np.array(p_lats[(n_pairs):])[::-1]
            tdep = np.array(p_depths[0:n_pairs])
            bdep = np.array(p_depths[(n_pairs):])[::-1]

            toplons = np.append(toplons, tlon)
            toplats = np.append(toplats, tlat)
            topdeps = np.append(topdeps, tdep)
            botlons = np.append(botlons, blon)
            botlats = np.append(botlats, blat)
            botdeps = np.append(botdeps, bdep)

            group_index.extend([g_ind] * n_pairs)
            g_ind = g_ind + 1

        reference = d['metadata']['reference']

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v
        d['metadata']['mesh_dx'] = mesh_dx

        self._geojson = d

        self._toplons = np.array(toplons)
        self._toplats = np.array(toplats)
        self._topdeps = np.array(topdeps)
        self._botlons = np.array(botlons)
        self._botlats = np.array(botlats)
        self._botdeps = np.array(botdeps)
        self._origin = origin
        self._group_index = np.array(group_index)
        self._mesh_dx = mesh_dx
        self._reference = reference
        self._computeStikeDip()

    @classmethod
    def fromArrays(cls, toplons, toplats, topdeps, botlons, botlats, botdeps,
                   origin, group_index=None, mesh_dx=0.5, reference=''):
        """
        Class method to initialize an EdgeRupture class from arrays of
        longitude, latitude, and depth along the top and bottom edges.

        Args:
            toplons (ndarray): Array of top edge longitudes.
            toplats (ndarray): Array of top edge latitudes.
            topdeps (ndarray): Array of top edge depths (km).
            botlons (ndarray): Array of bot edge longitudes.
            botlats (ndarray): Array of bot edge latitudes.
            botdeps (ndarray): Array of bot edge depths (km).
            origin (Origin): Reference to a ShakeMap Origin object.
            group_index (ndarray): Optional array of group index.
                If None, then assume only single group.
            mesh_dx (float): Target spacing (in km) for rupture discretization.
            reference (str): Citable reference for rupture.

        Returns:
            EdgeRupture instance.

        """
        toplons = np.array(toplons)
        toplats = np.array(toplats)
        topdeps = np.array(topdeps)
        botlons = np.array(botlons)
        botlats = np.array(botlats)
        botdeps = np.array(botdeps)
        if group_index is not None:
            group_index = np.array(group_index)
        else:
            group_index = np.zeros_like(toplons)

        coords = []
        u_groups = np.unique(group_index)
        n_groups = len(u_groups)
        for i in range(n_groups):
            ind = np.where(u_groups[i] == group_index)[0]
            lons = np.concatenate([toplons[ind], botlons[ind][::-1],
                                   toplons[ind][0].reshape((1,))])
            lats = np.concatenate([toplats[ind], botlats[ind][::-1],
                                   toplats[ind][0].reshape((1,))])
            deps = np.concatenate([topdeps[ind], botdeps[ind][::-1],
                                   topdeps[ind][0].reshape((1,))])
            poly = []
            for lon, lat, dep in zip(lons, lats, deps):
                poly.append([lon, lat, dep])
            coords.append(poly)

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
                     "type": "MultiPolygon",
                     "coordinates": [coords]
                 }
             }]}

        # Add origin information to metadata
        odict = origin.__dict__
        for k, v in odict.items():
            if isinstance(v, HistoricTime):
                d['metadata'][k] = v.strftime('%Y-%m-%dT%H:%M:%SZ')
            else:
                d['metadata'][k] = v
        d['metadata']['mesh_dx'] = mesh_dx

        return cls(d, origin, mesh_dx)

    def getLength(self):
        """
        Compute length of rupture (km). For EdgeRupture, we compute the length
        as the length of the top edge projected to the surface.

        Returns:
            float: Rupture length in km.
        """
        lons = self._toplons
        lats = self._toplats
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        rlength = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(lons[ind], lats[ind])
                P1 = Point(lons[ind + 1], lats[ind + 1])
                dist = P0.distance(P1)
                rlength = rlength + dist
        return rlength

    def getWidth(self):
        """
        Compute average rupture width in km. For EdgeRupture, we compute this
        as the rupture area divided by teh rupture length.

        Returns:
            float: Rupture width in km.
        """
        area = self.getArea()
        length = self.getLength()
        return area / length

    def getArea(self):
        """
        Compute the rupture area. For EdgeRupture, we compute this by grouping
        the traces into "quadrilaterals" for which the verticies may not be
        co-planar. We then break up the quadrilaterals into triangles for which
        we can compute area.

        Returns:
            float: Rupture area in square km.
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        area = 0
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                p0 = latlon2ecef(self._toplats[ind],
                                 self._toplons[ind],
                                 self._topdeps[ind])
                p1 = latlon2ecef(self._toplats[ind + 1],
                                 self._toplons[ind + 1],
                                 self._topdeps[ind + 1])
                p2 = latlon2ecef(self._botlats[ind + 1],
                                 self._botlons[ind + 1],
                                 self._botdeps[ind + 1])
                p3 = latlon2ecef(self._botlats[ind],
                                 self._botlons[ind],
                                 self._botdeps[ind])
                a = np.sqrt((p1[0] - p0[0])**2 +
                            (p1[1] - p0[1])**2 +
                            (p1[2] - p0[2])**2)
                b = np.sqrt((p2[0] - p0[0])**2 +
                            (p2[1] - p0[1])**2 +
                            (p2[2] - p0[2])**2)
                c = np.sqrt((p2[0] - p1[0])**2 +
                            (p2[1] - p1[1])**2 +
                            (p2[2] - p1[2])**2)
                s = (a + b + c) / 2
                A1 = np.sqrt(s * (s - a) * (s - b) * (s - c))
                a = np.sqrt((p0[0] - p3[0])**2 +
                            (p0[1] - p3[1])**2 +
                            (p0[2] - p3[2])**2)
                b = np.sqrt((p2[0] - p3[0])**2 +
                            (p2[1] - p3[1])**2 +
                            (p2[2] - p3[2])**2)
                c = np.sqrt((p0[0] - p2[0])**2 +
                            (p0[1] - p2[1])**2 +
                            (p0[2] - p2[2])**2)
                s = (a + b + c) / 2
                A2 = np.sqrt(s * (s - a) * (s - b) * (s - c))
                area = area + (A1 + A2) / 1000 / 1000
        return area

    def getStrike(self):
        """
        Return representative strike for this rupture. Note that strike
        can vary along the rupture.

        Returns:
            float: Strike angle in degrees.
        """
        return self._strike

    def getDip(self):
        """
        Representative dip for this rupture. Note that dip
        can vary along the rupture.

        Returns:
            float: dip angle in degrees.
        """
        return self._dip

    @property
    def lats(self):
        """
        Return an array of latitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans.
        """
        lats = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_lats = self._toplats[groups == u_groups[i]]
            top0 = top_lats[0].reshape((1,))
            bot_lats = self._botlats[groups == u_groups[i]]
            lats = np.concatenate((lats, top_lats, bot_lats[::-1], top0, nan))
        return np.array(lats)

    @property
    def lons(self):
        """
        Return an array of longitudes for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop longitude values; disconnected
                segments are separated by nans.
        """
        lons = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_lons = self._toplons[groups == u_groups[i]]
            top0 = top_lons[0].reshape((1,))
            bot_lons = self._botlons[groups == u_groups[i]]
            lons = np.concatenate((lons, top_lons, bot_lons[::-1], top0, nan))
        return np.array(lons)

    @property
    def depths(self):
        """
        Return an array of depths for the rupture verticies arranged for
        plotting purposes; will give an outline of each group connected
        segments.

        Returns:
            array: Numpy array of closed-loop latitude values; disconnected
                segments are separated by nans.
        """
        deps = np.empty(shape=(0,))
        groups = self._group_index
        u_groups = np.unique(groups)
        ng = len(u_groups)
        nan = np.array(np.nan).reshape(1,)
        for i in range(ng):
            top_deps = self._topdeps[groups == u_groups[i]]
            top0 = top_deps[0].reshape((1,))
            bot_deps = self._botdeps[groups == u_groups[i]]
            deps = np.concatenate((deps, top_deps, bot_deps[::-1], top0, nan))
        return np.array(deps)

    def _getGroupIndex(self):
        """
        Return a list of segment group indexes.

        Returns:
            list: Segment group indexes; length equals the number of
                quadrilaterals.
        """
        return copy.deepcopy(self._group_index)

    def _computeStikeDip(self):
        """
        Loop over all triangles and get the average normal, north, and up
        vectors in ECEF. Use these to compute a representative strike and dip.
        """
        seg = self._group_index
        groups = np.unique(seg)
        ng = len(groups)
        norm_vec = Vector(0, 0, 0)
        north_vec = Vector(0, 0, 0)
        up_vec = Vector(0, 0, 0)
        for i in range(ng):
            group_segments = np.where(groups[i] == seg)[0]
            nseg = len(group_segments) - 1
            for j in range(nseg):
                ind = group_segments[j]
                P0 = Point(self._toplons[ind],
                           self._toplats[ind],
                           self._topdeps[ind])
                P1 = Point(self._toplons[ind + 1],
                           self._toplats[ind + 1],
                           self._topdeps[ind + 1])
                P2 = Point(self._botlons[ind + 1],
                           self._botlats[ind + 1],
                           self._botdeps[ind + 1])
                P3 = Point(self._botlons[ind],
                           self._botlats[ind],
                           self._botdeps[ind])
                P1up = Point(self._toplons[ind + 1],
                             self._toplats[ind + 1],
                             self._topdeps[ind + 1] - 1.0)
                P1N = Point(self._toplons[ind + 1],
                            self._toplats[ind + 1] + 0.001,
                            self._topdeps[ind + 1])
                P3up = Point(self._botlons[ind],
                             self._botlats[ind],
                             self._botdeps[ind] - 1.0)
                P3N = Point(self._botlons[ind],
                            self._botlats[ind] + 0.001,
                            self._botdeps[ind])
                p0 = Vector.fromPoint(P0)
                p1 = Vector.fromPoint(P1)
                p2 = Vector.fromPoint(P2)
                p3 = Vector.fromPoint(P3)
                p1up = Vector.fromPoint(P1up)
                p1N = Vector.fromPoint(P1N)
                p3up = Vector.fromPoint(P3up)
                p3N = Vector.fromPoint(P3N)

                # Sides
                s01 = p1 - p0
                s02 = p2 - p0
                s03 = p3 - p0
                s21 = p1 - p2
                s23 = p3 - p2

                # First triangle
                t1norm = (s02.cross(s01)).norm()
                a = s01.mag()
                b = s02.mag()
                c = s21.mag()
                s = (a + b + c) / 2
                A1 = np.sqrt(s * (s - a) * (s - b) * (s - c)) / 1000

                # Second triangle
                t2norm = (s03.cross(s02)).norm()
                a = s03.mag()
                b = s23.mag()
                c = s02.mag()
                s = (a + b + c) / 2
                A2 = np.sqrt(s * (s - a) * (s - b) * (s - c)) / 1000

                # Up and North
                p1up = (p1up - p1).norm()
                p3up = (p3up - p3).norm()
                p1N = (p1N - p1).norm()
                p3N = (p3N - p3).norm()

                # Combine
                norm_vec = norm_vec + A1 * t1norm + A2 * t2norm
                north_vec = north_vec + A1 * p1N + A2 * p3N
                up_vec = up_vec + A1 * p1up + A2 * p3up

        norm_vec = norm_vec.norm()
        north_vec = north_vec.norm()
        up_vec = up_vec.norm()

        # Do I need to flip the vector because it is pointing down (i.e.,
        # right-hand rule is violated)?
        flip = np.sign(up_vec.dot(norm_vec))
        norm_vec = flip * norm_vec

        # Angle between up_vec and norm_vec is dip
        self._dip = np.arcsin(up_vec.cross(norm_vec).mag()) * 180 / np.pi

        # Normal vector projected to horizontal plane
        nvph = (norm_vec - up_vec.dot(norm_vec) * up_vec).norm()

        # Dip direction is angle between nvph and north; strike is orthogonal.
        cp = nvph.cross(north_vec)
        sign = np.sign(cp.dot(up_vec))
        dp = nvph.dot(north_vec)
        strike = np.arctan2(sign * cp.mag(), dp) * 180 / np.pi - 90
        if strike < -180:
            strike = strike + 360
        self._strike = strike

    def getDepthToTop(self):
        """
        Returns:
            float: Depth to top of rupture in km.
        """
        return np.min(self._topdeps)

    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals. Unlike QuadRupture, these
        quadrilaterals are not restricted to be coplanar or have
        horizontal top/bottom edges.

        Return:
            list: List of quadrilaterals, where each quadrilateral is
                a list of OQ Points.
        """
        ugroup = np.unique(self._group_index)
        ngroup = len(ugroup)
        qlist = []
        for i in range(ngroup):
            ind = np.where(self._group_index == ugroup[i])[0]
            nq = len(ind) - 1
            for j in range(nq):
                P0 = Point(self._toplons[j],
                           self._toplats[j],
                           self._topdeps[j])
                P1 = Point(self._toplons[j + 1],
                           self._toplats[j + 1],
                           self._topdeps[j + 1])
                P2 = Point(self._botlons[j + 1],
                           self._botlats[j + 1],
                           self._botdeps[j + 1])
                P3 = Point(self._botlons[j],
                           self._botlats[j],
                           self._botdeps[j])
                qlist.append([P0, P1, P2, P3])

        return qlist

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

        mesh_dx = self._mesh_dx

        # ---------------------------------------------------------------------
        # Sort out sites
        # ---------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        # ---------------------------------------------------------------------
        # Get mesh
        # ---------------------------------------------------------------------
        mx = []
        my = []
        mz = []
        u_groups = np.unique(self._group_index)
        n_groups = len(u_groups)
        for j in range(n_groups):
            g_ind = np.where(u_groups[j] == self._group_index)[0]
            nq = len(self._toplats[g_ind]) - 1
            for i in range(nq):
                q = [Point(self._toplons[g_ind[i]],
                           self._toplats[g_ind[i]],
                           self._topdeps[g_ind[i]]),
                     Point(self._toplons[g_ind[i + 1]],
                           self._toplats[g_ind[i + 1]],
                           self._topdeps[g_ind[i + 1]]),
                     Point(self._botlons[g_ind[i + 1]],
                           self._botlats[g_ind[i + 1]],
                           self._botdeps[g_ind[i + 1]]),
                     Point(self._botlons[g_ind[i]],
                           self._botlats[g_ind[i]],
                           self._botdeps[g_ind[i]])
                     ]
                mesh = utils.get_quad_mesh(q, dx=mesh_dx)
                mx.extend(list(np.reshape(mesh['x'], (-1,))))
                my.extend(list(np.reshape(mesh['y'], (-1,))))
                mz.extend(list(np.reshape(mesh['z'], (-1,))))
        mesh_mat = np.array([np.array(mx), np.array(my), np.array(mz)])

        # ---------------------------------------------------------------------
        # Compute distance
        # ---------------------------------------------------------------------
        dist = np.zeros_like(x)
        for i in range(len(x)):
            sitecol = sites_ecef[i, :].reshape([3, 1])
            dif = sitecol - mesh_mat
            distarray = np.sqrt(np.sum(dif * dif, axis=0))
            dist[i] = np.min(distarray) / 1000.0  # convert to km

        dist = np.reshape(dist, oldshape)

        return dist

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

        mesh_dx = self._mesh_dx

        # ---------------------------------------------------------------------
        # Sort out sites
        # ---------------------------------------------------------------------
        oldshape = lon.shape

        if len(oldshape) == 2:
            newshape = (oldshape[0] * oldshape[1], 1)
        else:
            newshape = (oldshape[0], 1)

        x, y, z = latlon2ecef(lat, lon, depth)
        x.shape = newshape
        y.shape = newshape
        z.shape = newshape
        sites_ecef = np.hstack((x, y, z))

        # ---------------------------------------------------------------------
        # Get mesh
        # ---------------------------------------------------------------------
        mx = []
        my = []
        mz = []
        u_groups = np.unique(self._group_index)
        n_groups = len(u_groups)
        for j in range(n_groups):
            g_ind = np.where(u_groups[j] == self._group_index)[0]
            nq = len(self._toplats[g_ind]) - 1
            for i in range(nq):
                q = [Point(self._toplons[g_ind[i]],
                           self._toplats[g_ind[i]],
                           0),
                     Point(self._toplons[g_ind[i + 1]],
                           self._toplats[g_ind[i + 1]],
                           0),
                     Point(self._botlons[g_ind[i + 1]],
                           self._botlats[g_ind[i + 1]],
                           0),
                     Point(self._botlons[g_ind[i]],
                           self._botlats[g_ind[i]],
                           0)
                     ]
                mesh = utils.get_quad_mesh(q, dx=mesh_dx)
                mx.extend(list(np.reshape(mesh['x'], (-1,))))
                my.extend(list(np.reshape(mesh['y'], (-1,))))
                mz.extend(list(np.reshape(mesh['z'], (-1,))))
        mesh_mat = np.array([np.array(mx), np.array(my), np.array(mz)])

        # ---------------------------------------------------------------------
        # Compute distance
        # ---------------------------------------------------------------------
        dist = np.zeros_like(x)
        for i in range(len(x)):
            sitecol = sites_ecef[i, :].reshape([3, 1])
            dif = sitecol - mesh_mat
            distarray = np.sqrt(np.sum(dif * dif, axis=0))
            dist[i] = np.min(distarray) / 1000.0  # convert to km

        dist = np.reshape(dist, oldshape)

        return dist

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
        # This just hands off to the module-level method
        # NOTE: Not sure if the non-horizontal top edges of EdgeRupture will
        #       case problems. Should do some checking. It might be okay to
        #       bring quad vertices up to the surface in this case.
        dict = gc2._computeGC2(self, lon, lat, depth)
        return dict
