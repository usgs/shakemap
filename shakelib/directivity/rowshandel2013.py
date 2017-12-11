"""
Implements the Rowshandel (2013) directivity model. Note that the report
provides many options for computing the directivity factor and we have
also used some coefficients that are unpublished updates by Rowshandel.
Some of the implementation options are controlled by arguments (e.g.,
mtype).

Fd is the directivity function. One of the options is whether or not to
use the 'centering' term (the argument is 'centered' and defaults to
True). If it is True then

Fd = C1 * (XiPrime - Xic) * LD * DT * WP

otherwise

Fd = (C1 * Xiprime + C2) * LD * DT * WP

where

- C1 and C2 are regression coefficients.
- XiPrime is a factor that accounts for rupture propagation and slip
  direction.
- Xic is the centering term.
- LD is the rupture length de-normalization factor.
- DT is the distance-taper.
- WP is the narrow-band multiplier.

Note that Fd is intended to be used in GMPEs as:

ln(IM_dir) = ln(IM) + Fd

where

- Fd is the directivity effect.
- IM is the intensity measure predicted by the GMPE.
- IM_dir is the directivity-adjusted IM.

To do
    - Add checks on function arguments (e.g., mtype) for valid values.
    - Interpolate periods.

References:
    Rowshandel, B. (2013). Rowshandel’s NGA-West2 directivity model,
    Chapter 3 of PEER Report No. 2013/09, P. Spudich (Editor), Pacific
    Earthquake Engineering Research Center, Berkeley, CA.
    `[link] <http://peer.berkeley.edu/publications/peer_reports/reports_2013/webPEER-2013-09-Spudich.pdf>`__
"""  # noqa

import numpy as np
import copy

import openquake.hazardlib.geo as geo
from openquake.hazardlib.geo.utils import get_orthographic_projection

from shakelib.rupture import utils
from shakelib.distance import get_distance
from impactutils.vectorutils.ecef import latlon2ecef
from impactutils.vectorutils.ecef import ecef2latlon
from impactutils.vectorutils.vector import Vector


class Rowshandel2013(object):

    __c1 = [0.35, 0.75, 0.95, 1.28, 1.60]
    __c2 = [-0.035, -0.10, -0.15, -0.26, -0.30]
    __periods = [1.0, 3.0, 5.0, 7.5, 10.0]

    def __init__(self, origin, rup, lat, lon, dep, dx, T, a_weight=0.5,
                 mtype=1, simpleDT=False, centered=True):
        """
        Constructor for rowshandel2013.

        Args:
            origin: Origin instance.
            rup: Rupture instance.
            lat (ndarray): Numpy array of site latitudes.
            lon (ndarray): Numpy array of site longitudes.
            dep (ndarray): Numpy array of site depths (km); positive down.
            dx (float): Target mesh spacing for subruptures in km. The mesh
                snaps to the edges of the quadrilaterals so the actual mesh
                spacing will not equal this exactly; spacing in x and y will
                not be equal.
            T (list): List floats (or a float) for periods to compute fd;
                Currently, only acceptable values are 1, 2, 3, 4, 5, 7.5.
            a_weight (float): Weighting factor for how p-dot-q and s-dot-q are
                averaged; 0 for only p-dot-q (propagation factor) and 1 for
                only s-dot-q (slip factor).
            mtype (int): Either 1 or 2; 1 for adding only positive components
                of dot products, 2 for adding all components of dot products.
            simpleDT (bool): Should the simpler DT equation be used? Usually
                False.
            centered (bool): Should the centered directivity parameter be used?
        """
        self._origin = origin
        self._rup = rup
        self._hyp = origin.getHypo()
        self._lat = lat
        self._lon = lon
        self._dep = dep
        self._rake = origin.rake
        self._dx = dx
        self._M = origin.mag
        if isinstance(T, float):
            self._T = [T]
        else:
            self._T = T
        if (a_weight > 1.0) or (a_weight < 0):
            raise ValueError('a_weight must be between 0 and 1.')
        self._a_weight = a_weight
        if (mtype != 1) and (mtype != 2):
            raise ValueError('mtype can onlty be 1 or 2.')
        self._mtype = mtype
        self._simpleDT = simpleDT
        self._centered = centered

        # Period independent parameters
        self.__computeWrup()
        self.__computeLD()
        self.__computeXiPrime()
        self.__getCenteringTerm()

        self.__computeFd()

    @classmethod
    def fromSites(cls, origin, rup, sites, dx, T, a_weight=0.5,
                  mtype=1, simpleDT=False, centered=True):
        """Construct a rowshandel2013 instance from a sites instance.

        Class method for constructing a rowshandel2013 instance from
        a sites instance.

        Args:
            origin: Origin instance.
            rup: Rupture instance.
            sites: Sites object.
            dx: Float for target mesh spacing for subruptures in km. The mesh
                snaps to the edges of the quadrilaterals so the actual mesh
                spacing will not equal this exactly; spacing in x and y will
                not be equal.
            T: Period; Currently, only acceptable values are 1, 2, 3, 4, 5,
                7.5.
            a_weight: Weighting factor for how p-dot-q and s-dot-q are
                averaged; 0 for only p-dot-q (propagation factor) and 1 for
                only s-dot-q (slip factor).
            mtype: Integer, either 1 or 2; 1 for adding only positive
                components of dot products, 2 for adding all components of dot
                products.
            simpleDT: Boolean; should the simpler DT equation be used? Usually
                False.
            centered: Boolean; should the centered directivity parameter be
                used

        Returns:
            Rowshandel2013 directivity class.
        """
        sm_dict = sites.getVs30Grid().getGeoDict()
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
        return cls(origin, rup, lat, lon, dep, dx, T,
                   a_weight, mtype, simpleDT, centered)

    def getFd(self):
        """
        :returns:
            Numpy array of Fd; this is the directivity amplification factor
            (log units).
        """
        return copy.deepcopy(self._fd)

    def getXiPrime(self):
        """
        Returns:
            Numpy array of Xi'; this is the variable that accounts for
            geometric factors in Fd.
        """
        return copy.deepcopy(self._xi_prime)

    def getDT(self):
        """
        Returns:
           List of numpy arrays of the distance taper factor. Length is the
           number of periods.
        """
        return copy.deepcopy(self._DTlist)

    def getWP(self):
        """
        Returns:
            List of numpy array of the narrow-band multiplier. Length is the
            number of periods.
        """
        return copy.deepcopy(self._WPlist)

    def getLD(self):
        """
        Returns:
            Numpy array of the rupture length de-normalization factor.
        """
        return copy.deepcopy(self._LD)

    def __computeFd(self):
        """
        Fd is the term that modifies the GMPE output as an additive term
        (in log space).
        ln(Y) = f(M, R, ...) + Fd
        """

        xcipc = self._xi_prime - self._xi_c

        # fd is a list with same length as T
        self._fd = [None] * len(self._T)
        self._DTlist = [None] * len(self._T)
        self._WPlist = [None] * len(self._T)
        for i in range(len(self._T)):
            period = self._T[i]
            c1sel = self.__c1[self.__periods == period]
            c2sel = self.__c2[self.__periods == period]

            # Period dependent parameters
            self.__computeWP(period)
            self.__computeDT(period)

            self._DTlist[i] = copy.deepcopy(self._DT)
            self._WPlist[i] = copy.deepcopy(self._WP)

            if self._centered:
                self._fd[i] = c1sel * self._WP * self._DT * xcipc
            else:
                self._fd[i] = (c1sel * self._xi_prime + c2sel) * \
                    self._LD * self._DT * self._WP

    def __computeWrup(self):
        """
        Compute the the portion (in km) of the width of the rupture which
        ruptures up-dip from the hypocenter to the top of the rupture.

        Wrup is the portion (in km) of the width of the rupture which
        ruptures up-dip from the hypocenter to the top of the rupture.

        * This is ambiguous for ruptures with varible top of rupture (not
          allowed in NGA). For now, lets just compute this for the
          quad where the hypocenter is located.
        * Alternative is to compute max Wrup for the different quads.

        """
        nquad = len(self._rup.getQuadrilaterals())

        # ---------------------------------------------------------------------
        # First find which quad the hypocenter is on
        # ---------------------------------------------------------------------

        x, y, z = latlon2ecef(
            self._hyp.latitude, self._hyp.longitude, self._hyp.depth)
        hyp_ecef = np.array([[x, y, z]])
        qdist = np.zeros(nquad)
        for i in range(0, nquad):
            qdist[i] = utils._quad_distance(
                self._rup.getQuadrilaterals()[i], hyp_ecef)
        ind = int(np.where(qdist == np.min(qdist))[0][0])
        # *** check that this doesn't break with more than one quad
        q = self._rup.getQuadrilaterals()[ind]

        # ---------------------------------------------------------------------
        # Compute Wrup on that quad
        # ---------------------------------------------------------------------

        pp0 = Vector.fromPoint(geo.point.Point(
            q[0].longitude, q[0].latitude, q[0].depth))
        hyp_ecef = Vector.fromPoint(geo.point.Point(
            self._hyp.longitude, self._hyp.latitude, self._hyp.depth))
        hp0 = hyp_ecef - pp0
        ddv = utils.get_quad_down_dip_vector(q)
        self._Wrup = Vector.dot(ddv, hp0) / 1000

    def __computeXiPrime(self):
        """
        Computes the xi' value.
        """
        hypo_ecef = Vector.fromPoint(geo.point.Point(
            self._hyp.longitude, self._hyp.latitude, self._hyp.depth))

        slat = self._lat
        slon = self._lon

        # Convert site to ECEF:
        site_ecef_x = np.ones_like(slat)
        site_ecef_y = np.ones_like(slat)
        site_ecef_z = np.ones_like(slat)

        # Make a 3x(#number of sites) matrix of site locations
        # (rows are x, y, z) in ECEF
        site_ecef_x, site_ecef_y, site_ecef_z = latlon2ecef(
            slat, slon, np.zeros(slon.shape))
        site_mat = np.array([np.reshape(site_ecef_x, (-1,)),
                             np.reshape(site_ecef_y, (-1,)),
                             np.reshape(site_ecef_z, (-1,))])

        xi_prime_unscaled = np.zeros_like(slat)

        # Normalize by total number of subruptures. For mtype == 1, the number
        # of subruptures will vary with site and be different for xi_s and
        # xi_p, so keep two variables and sum them for each quad.
        nsubs = np.zeros(np.product(slat.shape))
        nsubp = np.zeros(np.product(slat.shape))

        xi_prime_s = np.zeros(np.product(slat.shape))
        xi_prime_p = np.zeros(np.product(slat.shape))

        for k in range(len(self._rup.getQuadrilaterals())):
            # Select a quad
            q = self._rup.getQuadrilaterals()[k]

            # Quad mesh (ECEF coords)
            mesh = utils.get_quad_mesh(q, self._dx)

            # Rupture plane normal vector (ECEF coords)
            rpnv = utils.get_quad_normal(q)
            rpnvcol = np.array([[rpnv.x],
                                [rpnv.y],
                                [rpnv.z]])

            cp_mat = np.array([np.reshape(mesh['cpx'], (-1,)),
                               np.reshape(mesh['cpy'], (-1,)),
                               np.reshape(mesh['cpz'], (-1,))])

            # Compute matrix of p vectors
            hypcol = np.array([[hypo_ecef.x],
                               [hypo_ecef.y],
                               [hypo_ecef.z]])
            pmat = cp_mat - hypcol

            # Project pmat onto quad
            ndotp = np.sum(pmat * rpnvcol, axis=0)
            pmat = pmat - ndotp * rpnvcol

            mag = np.sqrt(np.sum(pmat * pmat, axis=0))
            pmatnorm = pmat / mag  # like r1

            # According to Rowshandel:
            #   "The choice of the +/- sign in the above equations
            #    depends on the (along-the-strike and across-the-dip)
            #    location of the rupturing sub-fault relative to the
            #    location of the hypocenter."
            # and:
            #   "for the along the strike component of the slip unit
            #    vector, the choice of the sign should result in the
            #    slip unit vector (s) being exactly the same as  the
            #    rupture unit vector (p) for a pure strike-slip case"

            # Strike slip and dip slip components of unit slip vector
            # (ECEF coords)
            ds_mat, ss_mat = _get_quad_slip_ds_ss(
                q, self._rake, cp_mat, pmatnorm)

            slpmat = (ds_mat + ss_mat)
            mag = np.sqrt(np.sum(slpmat * slpmat, axis=0))
            slpmatnorm = slpmat / mag

            # Loop over sites
            for i in range(site_mat.shape[1]):
                sitecol = np.array([[site_mat[0, i]],
                                    [site_mat[1, i]],
                                    [site_mat[2, i]]])

                qmat = sitecol - cp_mat  # 3x(ni*nj), like r2
                mag = np.sqrt(np.sum(qmat * qmat, axis=0))
                qmatnorm = qmat / mag

                # Propagation dot product
                pdotqraw = np.sum(pmatnorm * qmatnorm, axis=0)

                # Slip vector dot product
                sdotqraw = np.sum(slpmatnorm * qmatnorm, axis=0)

                if self._mtype == 1:
                    # Only sum over (+) directivity effect subruptures

                    # xi_p_prime
                    pdotq = pdotqraw.clip(min=0)
                    nsubp[i] = nsubp[i] + np.sum(pdotq > 0)

                    # xi_s_prime
                    sdotq = sdotqraw.clip(min=0)
                    nsubs[i] = nsubs[i] + np.sum(sdotq > 0)

                elif self._mtype == 2:
                    # Sum over contributing subruptures

                    # xi_p_prime
                    pdotq = pdotqraw
                    nsubp[i] = nsubp[i] + cp_mat.shape[1]

                    # xi_s_prime
                    sdotq = sdotqraw
                    nsubs[i] = nsubs[i] + cp_mat.shape[1]

                # Normalize by n sub ruptures later
                xi_prime_s[i] = xi_prime_s[i] + np.sum(sdotq)
                xi_prime_p[i] = xi_prime_p[i] + np.sum(pdotq)

        # Apply a water level to nsubp and nsubs to avoid division by
        # zero. This should only occur when the numerator is also zero
        # and so the resulting value should be zero.
        nsubs = np.maximum(nsubs, 1)
        nsubp = np.maximum(nsubp, 1)

        # We are outside the 'k' loop over nquads.
        # o Normalize xi_prime_s and xi_prime_p
        # o Reshape them
        # o Add them together with the 'a' weights
        xi_prime_tmp = (self._a_weight) * (xi_prime_s / nsubs) + \
                       (1 - self._a_weight) * (xi_prime_p / nsubp)
        xi_prime_unscaled = xi_prime_unscaled + \
            np.reshape(xi_prime_tmp, slat.shape)

        # Scale so that xi_prime has range (0, 1)
        if self._mtype == 1:
            xi_prime = xi_prime_unscaled
        elif self._mtype == 2:
            xi_prime = 0.5 * (xi_prime_unscaled + 1)

        self._xi_prime = xi_prime

    def __computeLD(self):
        """
        Computes LD -- the rupture length de-normalization term.
        """
        # Use an orthographic projection
        latmin = self._lat.min()
        latmax = self._lat.max()
        lonmin = self._lon.min()
        lonmax = self._lon.max()
        proj = geo.utils.get_orthographic_projection(
            lonmin, lonmax, latmax, latmin)

        # Get epi projection
        epi_x, epi_y = proj(self._hyp.longitude, self._hyp.latitude)

        # Get the lines for the top edge of the rupture
        qds = self._rup.getQuadrilaterals()
        nq = len(qds)
        top_lat = np.array([[]])
        top_lon = np.array([[]])
        top_dep = np.array([[]])
        for j in range(0, nq):
            top_lat = np.append(top_lat, qds[j][0].latitude)
            top_lat = np.append(top_lat, qds[j][1].latitude)
            top_lon = np.append(top_lon, qds[j][0].longitude)
            top_lon = np.append(top_lon, qds[j][1].longitude)
            top_dep = np.append(top_dep, qds[j][0].depth)
            top_dep = np.append(top_dep, qds[j][1].depth)
        top_x, top_y = proj(top_lon, top_lat)

        Lrup_max = 400
        slat = self._lat
        slon = self._lon
        LD = np.zeros_like(slat)
        Ls = np.zeros_like(slat)

        ni, nj = LD.shape
        for i in range(ni):
            for j in range(nj):
                # -------------------------------------------------------------
                # Compute Ls
                # -------------------------------------------------------------
                # Convert to local orthographic
                site_x, site_y = proj(slon[i, j], slat[i, j])

                # Shift so center is at epicenter
                site_x2 = site_x - epi_x
                site_y2 = site_y - epi_y
                top_x2 = top_x - epi_x
                top_y2 = top_y - epi_y

                # Angle to rotate to put site on x-axis
                alpha = np.arctan2(site_y2, site_x2)
                axis = [0, 0, 1]
                rmat = _rotation_matrix(axis, -alpha)
                llr = [None] * len(top_lat)

                # Apply rotation to each point on trace
                for k in range(len(top_lat)):
                    llr[k] = np.dot(rmat, [top_x2[k], top_y2[k], 0])
                    site3 = np.dot(rmat, [site_x2, site_y2, 0])
                Li = np.min([np.max([a[0] for a in llr]), site3[0]])

                # -------------------------------------------------------------
                # Compute LD and save results into matrices
                # -------------------------------------------------------------

                Lrup = np.sqrt(Li * Li + self._Wrup * self._Wrup)
                LD[i, j] = np.log(Lrup) / np.log(Lrup_max)
                Ls[i, j] = Li
        self._LD = LD
        self._Ls = Ls

    def __computeDT(self, period):
        """
        Computes DT -- the distance taper term.
        """
        slat = self._lat
        slon = self._lon
        site_z = np.zeros_like(slat)
        ddict = get_distance('rrup', slat, slon, site_z, self._rup)
        Rrup = np.reshape(ddict['rrup'], (-1, ))
        nsite = len(Rrup)

        if self._simpleDT:   # eqn 3.10
            R1 = 35
            R2 = 70
            DT = np.ones(nsite)
            ix = [(Rrup > R1) & (Rrup < R2)]
            DT[ix] = 2 - Rrup[ix] / R1
            DT[Rrup >= R2] = 0
        else:                  # eqn 3.9
            if period >= 1:
                R1 = 20 + 10 * np.log(period)
                R2 = 2 * (20 + 10 * np.log(period))
            else:
                R1 = 20
                R2 = 40
            DT = np.ones(nsite)
            # As written in report:
            # DT[(Rrup > R1) & (Rrup < R2)] = \
            #    2 - Rrup[(Rrup > R1) & (Rrup < R2)]/(20 + 10 * np.log(period))
            # Modification:
            DT[(Rrup > R1) & (Rrup < R2)] = 2 - \
                Rrup[(Rrup > R1) & (Rrup < R2)] / R1
            # Note: it is not clear if the above modification is 'correct'
            #       but it gives results that make more sense
            DT[Rrup >= R2] = 0
        DT = np.reshape(DT, slat.shape)
        self._DT = DT

    def __getCenteringTerm(self):
        L = self._rup.getLength()
        CSSLP = -0.32 + 0.15 * np.log(L)
        CTRST = -0.20 + 0.09 * np.log(L)
        CNRML = -0.20 + 0.08 * np.log(L)
        if self._rake > 0:
            act = (CSSLP - CTRST) / 2
            cct = (CSSLP + CTRST) / 2
            self._xi_c = act * np.cos(2 * np.radians(self._rake)) + cct
        else:
            act = (CSSLP - CNRML) / 2
            cct = (CSSLP + CNRML) / 2
            self._xi_c = act * np.cos(2 * np.radians(self._rake)) + cct

    def __computeWP(self, period):
        """
        Computes WP -- the narrow-band multiplier.
        """
        # As written in report:
        # sig = 0.6
        # Tp = np.exp(1.27*self._M - 7.28)
        # Update (Rowshandel, personal communication)
        sig = 0.55
        Tp = np.exp(1.04 * self._M - 6.0)
        self._WP = np.exp(-(np.log10(period / Tp) *
                            np.log10(period / Tp)) / (2 * sig * sig))


def _get_quad_slip_ds_ss(q, rake, cp, p):
    """
    Compute the DIP SLIP and STRIKE SLIP components of the unit slip vector in
    ECEF coords for a quad and rake angle.
    :param q:
        A quadrilateral.
    :param rake:
        Direction of motion of the hanging wall relative to the
        foot wall, as measured by the angle (deg) from the strike vector.
    :param cp:
        A 3x(n sub rupture) array giving center points of each sub rupture
        in ECEF coords.
    :param p:
        A 3x(n sub rupture) array giving the unit vectors of the propagation
        vector on each sub rupture in ECEF coords.
    :returns:
        List of dip slip and strike slip components (each is a matrix)
        of the unit slip vector in ECEF space.
    """
    # Get quad vertices, strike, dip
    P0, P1 = q[0:2]
    strike = P0.azimuth(P1)
    dip = utils.get_quad_dip(q)

    # Slip unit vectors in 'local' (i.e., z-up, x-east) coordinates
    d1_local = utils.get_local_unit_slip_vector_DS(strike, dip, rake)
    s1_local = utils.get_local_unit_slip_vector_SS(strike, dip, rake)

    # Convert to a column array
    d1_col = np.array([[d1_local.x],
                       [d1_local.y],
                       [d1_local.z]])
    s1_col = np.array([[s1_local.x],
                       [s1_local.y],
                       [s1_local.z]])

    # Define 'local' coordinate system
    qlats = [a.latitude for a in q]
    qlons = [a.longitude for a in q]
    proj = get_orthographic_projection(
        np.min(qlons), np.max(qlons), np.min(qlats), np.max(qlats))

    # Convert p and cp to geographic coords
    p0lat, p0lon, p0z = ecef2latlon(cp[0, ], cp[1, ], cp[2, ])
    p1lat, p1lon, p1z = ecef2latlon(cp[0, ] + p[0, ],
                                    cp[1, ] + p[1, ],
                                    cp[2, ] + p[2, ])

    # Convert p to 'local' coords
    p0x, p0y = proj(p0lon, p0lat)
    p1x, p1y = proj(p1lon, p1lat)
    px = p1x - p0x
    py = p1y - p0y
    pz = p1z - p0z

    # Apply sign changes in 'local' coords
    s1mat = np.array([[np.abs(s1_col[0]) * np.sign(px)],
                      [np.abs(s1_col[1]) * np.sign(py)],
                      [np.abs(s1_col[2]) * np.sign(pz)]])
#                      [np.abs(s1_col[2])*np.ones_like(pz)]])

    dipsign = -np.sign(np.sin(np.radians(rake)))
    d1mat = np.array([[d1_col[0] * np.ones_like(px) * dipsign],
                      [d1_col[1] * np.ones_like(py) * dipsign],
                      [d1_col[2] * np.ones_like(pz) * dipsign]])

    # Need to track 'origin'
    s0 = np.array([[0], [0], [0]])

    # Convert from 'local' to geographic coords
    s1_ll = proj(s1mat[0, ], s1mat[1, ], reverse=True)
    d1_ll = proj(d1mat[0, ], d1mat[1, ], reverse=True)
    s0_ll = proj(s0[0], s0[1], reverse=True)

    # And then back to ECEF:
    s1_ecef = latlon2ecef(s1_ll[1], s1_ll[0], s1mat[2, ])
    d1_ecef = latlon2ecef(d1_ll[1], d1_ll[0], d1mat[2, ])
    s0_ecef = latlon2ecef(s0_ll[1], s0_ll[0], s0[2])
    s00 = s0_ecef[0].reshape(-1)
    s01 = s0_ecef[1].reshape(-1)
    s02 = s0_ecef[2].reshape(-1)
    d_mat = np.array([d1_ecef[0].reshape(-1) - s00,
                      d1_ecef[1].reshape(-1) - s01,
                      d1_ecef[2].reshape(-1) - s02])
    s_mat = np.array([s1_ecef[0].reshape(-1) - s00,
                      s1_ecef[1].reshape(-1) - s01,
                      s1_ecef[2].reshape(-1) - s02])
    return d_mat, s_mat


def _rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    Source: Response by 'unutbu' in this thread:
    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2)
    b, c, d = -axis * np.sin(theta / 2)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
