#!/usr/bin/env python

import numpy as np
import openquake.hazardlib.geo as geo
import copy

from shakelib.distance import get_distance
from shakelib.rupture.utils import _distance_sq_to_segment
import impactutils.vectorutils.ecef as ecef
from impactutils.vectorutils.vector import Vector


class Bayless2013(object):
    """
    Implements the Bayless and Somerville (2013) directivity model. This is an
    update to Somerville et al. (1997), which defines some of the input
    parameters.

    Fd is the directivity effect parameter, and it is evaluated as

    Fd = (c0 + c1 * Fgeom) * Tcd * Tmw * Taz

    Note that Fd is intended to be used in GMPEs as:

    ln(IM_dir) = ln(IM) + Fd

    where

    - Fd is the directivity effect
    - IM is the intensity measure predicted by the GMPE.
    - IM_dir is the directivity-adjusted IM.

    The model is seprated into strike-slip and dip-slip categories:
    SS: (abs(rake) >  0 & abs(rake) < 30) | (abs(rake) > 150 & abs(rake) < 180)
    DS:  abs(rake) > 60 & abs(rake) < 120

    Notes from Somerville et al. (1997):

    * d = length of dipping fault rupturing towards site
    * Y = d/W
    * s = length of striking fault rupturing towards site
    * X = s/L

    Bayless and Somerville state that each quadrilateral should have a pseudo-
    hypocenter:

        "Define the pseudo-hypocenter for rupture of successive segments as the
        point on the side edge of the fault segment that is closest to the side
        edge of the previous segment, and that lies half way between the top
        and bottom of the fault. We assume that the fault is segmented along
        strike, not updip. All geometric parameters are computed relative to
        the pseudo-hypocenter."

    To do:

        - Interpolate for arbitrary periods
        - Inherit from Directivity base class.

    References:
        Bayless, J., and Somerville, P. (2013). Bayless-Somerville Directivity
        Model, Chapter 3 of PEER Report No. 2013/09, P. Spudich (Editor),
        Pacific Earthquake Engineering Research Center, Berkeley, CA.
        `[link] <http://peer.berkeley.edu/publications/peer_reports/reports_2013/webPEER-2013-09-Spudich.pdf>`__

        Somerville, P. G., Smith, N. F., Graves, R. W., and Abrahamson, N. A.
        (1997). Modification of empirical strong ground motion attenuation
        relations to include the amplitude and duration effects of rupture
        directivity, Seismol. Res. Lett. 68, 199–222.
        `[link] <http://srl.geoscienceworld.org/content/68/1/199.abstract>`__

    """  # noqa
    # ------------------------------------------------------------------------
    # C0 adn C1 are for RotD50. One set for each mechanism (SS vs DS).
    # FN and FP are also available
    # ------------------------------------------------------------------------
    __periods = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10])
    __c0ss = np.array([0.0, 0.0, -0.12, -0.175, -0.21, -0.235, -0.255,
                       -0.275, -0.29, -0.3])
    __c1ss = np.array([0.0, 0.0, 0.075, 0.090, 0.095, 0.099, 0.103, 0.108,
                       0.112, 0.115])
    __c0ds = np.array([0.0, 0.0, 0.0, 0.0, 0.0, -0.033, -0.089, -0.133,
                       -0.16, -0.176])
    __c1ds = np.array([0.0, 0.0, 0.0, 0.0, 0.034, 0.093, 0.128, 0.15,
                       0.165, 0.179])

    def __init__(self, origin, rup, lat, lon, depth, T):
        """
        Constructor for bayless2013.

        Args:
            origin: Shakemap Origin object.
            rup: Shakemap Rupture object.
            lat: Numpy array of latitudes.
            lon: Numpy array of longitudes.
            depth: Numpy array of depths (km); down is positive.
            T: Period; Currently, only acceptable values are
                0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 7.5, 10.
        """
        self._origin = origin
        self._rup = rup
        self._rake = origin.rake
        self._M = origin.mag
        self._hyp = origin.getHypo()
        self._lon = lon
        self._lat = lat
        self._dep = depth
        self._T = T

        # Lists of widths and lengths for each quad in the rupture
        self._W = self._rup.getIndividualWidths()
        self._L = self._rup.getIndividualTopLengths()

        # Number of quads
        self._nq = len(self._W)

        # Currently assuming that the rake is the same on all subruptures .
        self.__getSlipCategory()

        # Rupture weights are supposed to be based on seismic moment.
        # Since moment is proportional to area, lets just use area
        # for now
        area = self._W * self._L
        self.weights = area / np.sum(area)

        # Put in pseudo-hypocenters for each quad
        self.__setPseudoHypocenters()

        self._fd = 0
        for i in range(self._nq):
            self.i = i

            # Compute some genral stuff that is required for all mechanisms
            dtypes = ['rrup', 'rx', 'ry0']
            dists = get_distance(dtypes, self._lat, self._lon, self._dep,
                                 self._rup)
            self.__Rrup = np.reshape(dists['rrup'], self._lat.shape)
            self.__Rx = np.reshape(dists['rx'], self._lat.shape)
            self.__Ry = np.reshape(dists['ry0'], self._lat.shape)
            # NOTE: use Rx and Ry to compute Az in 'computeAz'. It is probably
            #       possible to make this a lot faster by avoiding the
            #       calculation of these distances each time.

            # Az is the NGA definition of source-to-site azimuth for a finite
            # rupture. See Kaklamanos et al. (2011) Figure 2 for illustration.

            self.__computeAz()  # uses Rx and Ry, which are for the i-th quad.

            # Magnitude taper (does not depend on mechanism)
            if self._M <= 5.0:
                self._T_Mw = 0.0
            elif (self._M > 5.0) and (self._M < 6.5):
                self._T_Mw = 1.0 - (6.5 - self._M) / 1.5
            else:
                self._T_Mw = 1.0

            if self.SlipCategory == 'SS':
                self.__computeSS()
                self._fd = self._fd + self.weights[i] * self._fd_SS

            elif self.SlipCategory == 'DS':
                self.__computeDS()
                self._fd = self._fd + self.weights[i] * self._fd_DS

            else:
                # Compute both SS and DS
                self.__computeSS()
                self.__computeDS()

                # Normalize rake to reference angle
                sintheta = np.abs(np.sin(np.radians(self._rake)))
                costheta = np.abs(np.cos(np.radians(self._rake)))
                refrake = np.arctan2(sintheta, costheta)

                # Compute weights:
                DipWeight = refrake / (np.pi / 2.0)
                StrikeWeight = 1.0 - DipWeight
                fdcombined = StrikeWeight * self._fd_SS + DipWeight * \
                    self._fd_DS
                self._fd = self._fd + self.weights[i] * fdcombined

    def __setPseudoHypocenters(self):
        """ Set a pseudo-hypocenter.

        Adapted from ShakeMap 3.5 src/contour/directivity.c

        From Bayless and Somerville:

        "Define the pseudo-hypocenter for rupture of successive segments as
        the point on the side edge of the rupture segment that is closest to
        the side edge of the previous segment, and that lies half way
        between the top and bottom of the rupture. We assume that the rupture
        is segmented along strike, not updip. All geometric parameters are
        computed relative to the pseudo-hypocenter."
        """
        hyp_ecef = Vector.fromPoint(geo.point.Point(
            self._hyp.longitude, self._hyp.latitude, self._hyp.depth))
        # Loop over each quad
        self.phyp = [None] * self._nq
        for i in range(self._nq):
            P0, P1, P2, P3 = self._rup.getQuadrilaterals()[i]
            p0 = Vector.fromPoint(P0)  # convert to ECEF
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)

            # Create 4 planes with normals pointing outside rectangle
            hpnp = Vector.cross(p1 - p0, p2 - p0).norm()
            hpp = -hpnp.x * p0.x - hpnp.y * p0.y - hpnp.z * p0.z
            n0 = Vector.cross(p1 - p0, hpnp)
            n1 = Vector.cross(p2 - p1, hpnp)
            n2 = Vector.cross(p3 - p2, hpnp)
            n3 = Vector.cross(p0 - p3, hpnp)

            # -----------------------------------------------------------------
            # Is the hypocenter inside the projected rectangle?
            # Dot products show which side the origin is on.
            # If origin is on same side of all the planes, then it is 'inside'
            # -----------------------------------------------------------------

            sgn0 = np.signbit(Vector.dot(n0, p0 - hyp_ecef))
            sgn1 = np.signbit(Vector.dot(n1, p1 - hyp_ecef))
            sgn2 = np.signbit(Vector.dot(n2, p2 - hyp_ecef))
            sgn3 = np.signbit(Vector.dot(n3, p3 - hyp_ecef))

            if (sgn0 == sgn1) and (sgn1 == sgn2) and (sgn2 == sgn3):
                # Origin is inside. Use distance-to-plane formula.
                d = Vector.dot(hpnp, hyp_ecef) + hpp
                d = d * d

                # Put the pseudo hypocenter on the plane
                D = Vector.dot(hpnp, hyp_ecef) + hpp
                self.phyp[i] = hyp_ecef - hpnp * D

            else:
                # Origin is outside. Find distance to edges
                p0p = np.reshape(p0.getArray() - hyp_ecef.getArray(), [1, 3])
                p1p = np.reshape(p1.getArray() - hyp_ecef.getArray(), [1, 3])
                p2p = np.reshape(p2.getArray() - hyp_ecef.getArray(), [1, 3])
                p3p = np.reshape(p3.getArray() - hyp_ecef.getArray(), [1, 3])
                s1 = _distance_sq_to_segment(p1p, p2p)
                s3 = _distance_sq_to_segment(p3p, p0p)

                # Assuming that the rupture is segmented along strike and not
                # updip (as described by Bayless and somerville), we only
                # need to consider s1 and s3:
                if s1 > s3:
                    e30 = p0 - p3
                    e30norm = e30.norm()
                    mag = e30.mag()
                    self.phyp[i] = p3 + e30norm * (0.5 * mag)
                else:
                    e21 = p1 - p2
                    e21norm = e21.norm()
                    mag = e21.mag()
                    self.phyp[i] = p2 + e21norm * (0.5 * mag)

    def __computeDS(self):
        # d is the length of dipping rupture rupturing toward site;
        # Note: max[(Y*W),exp(0)] -- just apply a min of 1?
        self.__computeD(self.i)

        # Geometric directivity predictor:
        RxoverW = (self.__Rx /
                   self._W[self.i]).clip(
            min=-np.pi / 2.0, max=2.0 * np.pi / 3.0)
        f_geom = np.log(self.d) * np.cos(RxoverW)

        # Distance taper
        T_CD = np.ones_like(self._lat)
        ix = [(self.__Rrup / self._W[self.i] > 1.5) &
              (self.__Rrup / self._W[self.i] < 2.0)]
        T_CD[ix] = 1.0 - (self.__Rrup[ix] / self._W[self.i] - 1.5) / 0.5
        T_CD[self.__Rrup / self._W[self.i] >= 2.0] = 0.0

        # Azimuth taper
        T_Az = np.sin(np.abs(self.Az))**2

        # Select Coefficients
        ix = [self._T == self.__periods]
        C0 = self.__c0ds[ix]
        C1 = self.__c1ds[ix]

        self._fd_DS = (C0 + C1 * f_geom) * T_CD * self._T_Mw * T_Az

    def __computeSS(self):
        # s is the length of striking fault rupturing toward site;
        # max[(X*L),exp(1)]
        # theta (see Figure 5 in SSGA97)
        self.__computeThetaAndS(self.i)

        # Geometric directivity predictor:
        f_geom = np.log(self.s) * (0.5 * np.cos(2 * self.theta) + 0.5)

        # Distance taper
        T_CD = np.ones_like(self._lat)
        ix = [
            (self.__Rrup /
             self._L[
                 self.i] > 0.5) & (
                self.__Rrup /
                self._L[
                    self.i] < 1.0)]
        T_CD[ix] = 1 - (self.__Rrup[ix] / self._L[self.i] - 0.5) / 0.5
        T_CD[self.__Rrup / self._L[self.i] >= 1.0] = 0.0

        # Azimuth taper
        T_Az = 1.0

        # Select Coefficients
        ix = [self._T == self.__periods]
        C0 = self.__c0ss[ix]
        C1 = self.__c1ss[ix]
        self._fd_SS = (C0 + C1 * f_geom) * T_CD * self._T_Mw * T_Az

    def __computeAz(self):
        Az = np.ones_like(self.__Rx) * np.pi / 2.0
        Az = Az * np.sign(self.__Rx)
        ix = [self.__Ry > 0.0]
        Az[ix] = np.arctan(self.__Rx[ix] / self.__Ry[ix])
        self.Az = Az

    def __computeD(self, i):
        """Compute d for the i-th quad/segment.

        Y = d/W, where d is the portion (in km) of the width of the fault which
        ruptures up-dip from the hypocenter to the top of the fault.

        Args:
            i (int): index of segment for which d is to be computed.
        """
        hyp_ecef = self.phyp[i]  # already in ECEF
        hyp_col = np.array([[hyp_ecef.x], [hyp_ecef.y], [hyp_ecef.z]])

        # First compute "updip" vector
        P0, P1, P2 = self._rup.getQuadrilaterals()[i][0:3]
        p1 = Vector.fromPoint(P1)  # convert to ECEF
        p2 = Vector.fromPoint(P2)
        e21 = p1 - p2
        e21norm = e21.norm()
        hp1 = p1 - hyp_ecef
        # convert to km (used as max later)
        udip_len = Vector.dot(hp1, e21norm) / 1000.0
        udip_col = np.array(
            [[e21norm.x], [e21norm.y], [e21norm.z]])  # ECEF coords

        # Sites
        slat = self._lat
        slon = self._lon

        # Convert sites to ECEF:
        site_ecef_x = np.ones_like(slat)
        site_ecef_y = np.ones_like(slat)
        site_ecef_z = np.ones_like(slat)

        # Make a 3x(#number of sites) matrix of site locations
        # (rows are x, y, z) in ECEF
        site_ecef_x, site_ecef_y, site_ecef_z = ecef.latlon2ecef(
            slat, slon, np.zeros(slon.shape))
        site_mat = np.array([np.reshape(site_ecef_x, (-1,)),
                             np.reshape(site_ecef_y, (-1,)),
                             np.reshape(site_ecef_z, (-1,))])

        # Hypocenter-to-site matrix
        h2s_mat = site_mat - hyp_col  # in ECEF

        # Dot hypocenter-to-site with updip vector
        d_raw = np.abs(np.sum(h2s_mat * udip_col, axis=0)) / \
            1000.0  # convert to km
        d_raw = np.reshape(d_raw, self._lat.shape)
        self.d = d_raw.clip(min=1.0, max=udip_len)

    def __computeThetaAndS(self, i):
        """
        Args:
            i (int): Compute d for the i-th quad/segment.
        """
        # self.phyp is in ECEF
        tmp = ecef.ecef2latlon(self.phyp[i].x, self.phyp[i].y, self.phyp[i].z)
        epi_ecef = Vector.fromPoint(geo.point.Point(tmp[1], tmp[0], 0.0))
        epi_col = np.array([[epi_ecef.x], [epi_ecef.y], [epi_ecef.z]])

        # First compute along strike vector
        P0, P1, P2, P3 = self._rup.getQuadrilaterals()[i]
        p0 = Vector.fromPoint(P0)  # convert to ECEF
        p1 = Vector.fromPoint(P1)
        e01 = p1 - p0
        e01norm = e01.norm()
        hp0 = p0 - epi_ecef
        hp1 = p1 - epi_ecef
        strike_min = Vector.dot(hp0, e01norm) / 1000.0  # convert to km
        strike_max = Vector.dot(hp1, e01norm) / 1000.0  # convert to km
        strike_col = np.array(
            [[e01norm.x], [e01norm.y], [e01norm.z]])  # ECEF coords

        # Sites
        slat = self._lat
        slon = self._lon

        # Convert sites to ECEF:
        site_ecef_x = np.ones_like(slat)
        site_ecef_y = np.ones_like(slat)
        site_ecef_z = np.ones_like(slat)

        # Make a 3x(#number of sites) matrix of site locations
        # (rows are x, y, z) in ECEF
        site_ecef_x, site_ecef_y, site_ecef_z = ecef.latlon2ecef(
            slat, slon, np.zeros(slon.shape))
        site_mat = np.array([np.reshape(site_ecef_x, (-1,)),
                             np.reshape(site_ecef_y, (-1,)),
                             np.reshape(site_ecef_z, (-1,))])

        # Epicenter-to-site matrix
        e2s_mat = site_mat - epi_col  # in ECEF
        mag = np.sqrt(np.sum(e2s_mat * e2s_mat, axis=0))

        # Avoid division by zero
        mag[mag == 0] = 1e-12
        e2s_norm = e2s_mat / mag

        # Dot epicenter-to-site with along-strike vector
        s_raw = np.sum(e2s_mat * strike_col, axis=0) / 1000.0  # conver to km

        # Put back into a 2d array
        s_raw = np.reshape(s_raw, self._lat.shape)
        self.s = np.abs(s_raw.clip(min=strike_min,
                                   max=strike_max)).clip(min=np.exp(1))

        # Compute theta
        sdots = np.sum(e2s_norm * strike_col, axis=0)
        theta_raw = np.arccos(sdots)

        # But theta is defined to be the reference angle
        # (i.e., the equivalent angle between 0 and 90 deg)
        sintheta = np.abs(np.sin(theta_raw))
        costheta = np.abs(np.cos(theta_raw))
        theta = np.arctan2(sintheta, costheta)
        self.theta = np.reshape(theta, self._lat.shape)

    def getFd(self):
        """
        Returns:
            Numpy array of Fd; Fd is the directivity amplification factor.
        """
        return copy.deepcopy(self._fd)

    def __getSlipCategory(self):
        """
        Sets self.SlipCategory based on rake angle. Can be SS for
        strike slip, DS for dip slip, or Unspecified.
        """
        arake = np.abs(self._rake)
        self.SlipCategory = 'Unspecified'
        if ((arake >= 0) and (arake <= 30)) or (
                (arake >= 150) and (arake <= 180)):
            self.SlipCategory = 'SS'
        if (arake >= 60) and (arake <= 120):
            self.SlipCategory = 'DS'
