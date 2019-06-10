import numpy as np
from scipy.interpolate import RegularGridInterpolator
import itertools as it

from shakelib.correlation.ccf_base import CrossCorrelationBase
from shakemap.c.clib import eval_lb_correlation

# Periods that apply to the axes in the cross-correlation tables
Tlist = np.array([0.01, 0.1, 0.2, 0.5, 1, 2, 5, 7.5, 10.0001])

# Table II. Short range coregionalization matrix, B1
B1 = np.array([
    [0.305, 0.24, 0.23,  0.22, 0.16, 0.07,  0.03, 0,    0],
    [0.24,  0.275, 0.19, 0.13, 0.08, 0,     0,    0,    0],
    [0.23,  0.19, 0.26,  0.19, 0.12, 0.04,  0,    0,    0],
    [0.22,  0.13, 0.19,  0.32, 0.23, 0.14,  0.09, 0.06, 0.04],
    [0.16,  0.08, 0.12,  0.23, 0.32, 0.22,  0.13, 0.09, 0.07],
    [0.07,  0,    0.04,  0.14, 0.22, 0.335, 0.23, 0.19, 0.16],
    [0.03,  0,    0,     0.09, 0.13, 0.23,  0.34, 0.29, 0.24],
    [0,     0,    0,     0.06, 0.09, 0.19,  0.29, 0.30, 0.25],
    [0,     0,    0,     0.04, 0.07, 0.16,  0.24, 0.25, 0.24]
])

# Table III. Long range coregionalization matrix, B2
B2 = np.array([
    [0.315, 0.26,  0.27, 0.24, 0.17, 0.11,  0.08, 0.06, 0.05],
    [0.26,  0.295, 0.22, 0.15, 0.07, 0,     0,    0,   -0.03],
    [0.27,  0.22,  0.29, 0.24, 0.15, 0.09,  0.03, 0.02, 0],
    [0.24,  0.15,  0.24, 0.33, 0.27, 0.23,  0.17, 0.14, 0.14],
    [0.17,  0.07,  0.15, 0.27, 0.38, 0.34,  0.23, 0.19, 0.21],
    [0.11,  0,     0.09, 0.23, 0.34, 0.445, 0.33, 0.29, 0.32],
    [0.08,  0,     0.03, 0.17, 0.23, 0.33,  0.45, 0.42, 0.42],
    [0.06,  0,     0.02, 0.14, 0.19, 0.29,  0.42, 0.47, 0.47],
    [0.05, -0.03,  0,    0.14, 0.21, 0.32,  0.42, 0.47, 0.54]
])

# Table IV. Nugget effect coregionalization matrix, B3
B3 = np.array([
    [0.38, 0.36,  0.35, 0.17,  0.04,  0.04,  0,     0.03, 0.08],
    [0.36, 0.43,  0.35, 0.13,  0,     0.02,  0,     0.02, 0.08],
    [0.35, 0.35,  0.45, 0.11, -0.04, -0.02, -0.04, -0.02, 0.03],
    [0.17, 0.13,  0.11, 0.35,  0.2,   0.06,  0.02,  0.04, 0.02],
    [0.04, 0,    -0.04, 0.20,  0.30,  0.14,  0.09,  0.12, 0.04],
    [0.04, 0.02, -0.02, 0.06,  0.14,  0.22,  0.12,  0.13, 0.09],
    [0,    0,    -0.04, 0.02,  0.09,  0.12,  0.21,  0.17, 0.13],
    [0.03, 0.02, -0.02, 0.04,  0.12,  0.13,  0.17,  0.23, 0.10],
    [0.08, 0.08,  0.03, 0.02,  0.04,  0.09,  0.13,  0.10, 0.22]
])


class LothBaker2013(CrossCorrelationBase):
    """
    Created by Christophe Loth, 12/18/2012
    Pythonized and vectorized by C. Bruce Worden, 3/15/2017
    Compute the spatial correlation of epsilons for the NGA ground motion
    models

    The function is strictly empirical, fitted over the range the range
    0.01s <= t1, t2 <= 10s

    Documentation is provided in the following paper:
    Loth, C., and Baker, J. W. (2013). “A spatial cross-correlation model of
    ground motion spectral accelerations at multiple periods.”
    Earthquake Engineering & Structural Dynamics, 42, 397-417.
    """

    def __init__(self, periods):
        """
        Create an instance of LB13.

        Args:
            periods (numpy.array): An array of periods that will be requested
                from the function. Values must be [0.01 -> 10.0], and must be
                sorted from smallest to largest.

        Returns:
            An instance of :class:`LothBaker2013`.
        """

        if np.any(periods < 0.01):
            raise ValueError('The periods must be greater or equal to 0.01s')
        if np.any(periods > 10):
            raise ValueError('The periods must be less or equal to 10s')

        nper = np.size(periods)

        rgi1 = RegularGridInterpolator((Tlist, Tlist), B1, method='linear')
        rgi2 = RegularGridInterpolator((Tlist, Tlist), B2, method='linear')
        rgi3 = RegularGridInterpolator((Tlist, Tlist), B3, method='linear')

        #
        # Build new tables with entries at the periods we will use
        #
        tlist = list(it.product(periods, periods))
        self.b1 = rgi1(tlist).reshape((nper, nper))
        self.b2 = rgi2(tlist).reshape((nper, nper))
        self.b3 = rgi3(tlist).reshape((nper, nper))

        db1 = np.diagonal(B1)
        db2 = np.diagonal(B2)
        db3 = np.diagonal(B3)

        rgid1 = RegularGridInterpolator((Tlist,), db1, method='linear')
        rgid2 = RegularGridInterpolator((Tlist,), db2, method='linear')
        rgid3 = RegularGridInterpolator((Tlist,), db3, method='linear')

        id1 = rgid1(periods)
        id2 = rgid2(periods)
        id3 = rgid3(periods)

        np.fill_diagonal(self.b1, id1)
        np.fill_diagonal(self.b2, id2)
        np.fill_diagonal(self.b3, id3)

    def getCorrelation(self, ix1, ix2, h):
        """
        Compute the correlation between two periods and a separation distance
        of h.

        The indices (ix1 and ix2) and h must have the same dimensions. The
        indices may be equal, and there is no restriction on which one is
        larger. The indices refer to periods in the 'period' argument to the
        class constructor. The result is stored in h.

        Args:
            ix1 (2D, C-contiguous numpy array)):
                The indices of the first period of interest.
            ix2 (2D, C-contiguous numpy array)):
                The indices of the second period of interest.
            h (2D, C-contiguous numpy array)):
                The separation distance between two sites (units of km).
                h will be returned with the result, so it must be
                copied if the values in h are to be preserved.

        Returns:
            h (numpy array):
                The predicted correlation coefficient. The output array
                will have the same shape as the inputs.

        """
        # Verify the validity of input arguments
        if np.any(h < 0):
            raise ValueError('The separation distance must be positive')
        if np.shape(ix1) != np.shape(ix2) != np.shape(h):
            raise ValueError(
                'The input arguments must all have the same dimensions')

        #
        # Index into the arrays to get the coefficients corresponding to the
        # periods of interest.
        #
        # These variables are used in ne.evaluate but unseen by linter
        # b1 = self.b1[ix1, ix2]  # noqa
        # b2 = self.b2[ix1, ix2]  # noqa
        # b3 = self.b3[ix1, ix2]  # noqa
        # afact = -3.0 / 20.0  # noqa
        # bfact = -3.0 / 70.0  # noqa
        #
        # Compute the correlation coefficient (Equation 42)
        #
        # rho = ne.evaluate(
        #     "b1 * exp(h * afact) + b2 * exp(h * bfact) + (h == 0) * b3")
        rho = eval_lb_correlation(self.b1, self.b2, self.b3, ix1, ix2, h)

        return rho
