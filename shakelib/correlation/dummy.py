import numpy as np

from shakelib.correlation.ccf_base import CrossCorrelationBase


class DummyCorrelation(CrossCorrelationBase):
    """
    Simplified correlation module for testing purposes. Should not be used
    in productions runs as it does not produce valid correlations.
    """

    def __init__(self, periods):
        """
        Initialize the cross-correlation object.

        Args:
            periods (ndarray): An array of periods that will be requested
                from the function. Values must be in the range [0.01, 10.0],
                and must me sorted from smallest to largest.

        Returns:
            An instance of :class:`DummyCorrelation`.
        """

        if np.any(periods < 0.01):
            raise ValueError("The periods must be greater or equal to 0.01s")
        if np.any(periods > 10):
            raise ValueError("The periods must be less or equal to 10s")

        self.periods = periods.copy()

    def getCorrelation(self, ix1, ix2, h):
        """
        Compute the correlation between two periods and a separation distance
        of h km. The result returned is::

          rho = T1/T2 * exp(-h/10)

        where rho is the correlation, T1 is the smaller period, T2 is the
        larger period, and h is the distance between the points of interest.

        The index arrays (ix1 and ix2) and h array must have the same
        dimensions. The
        indices may be equal, and there is no restriction on which one is
        larger. The indices refer to periods in the 'period' argument to the
        class constructor.

        Args:
            ix1 (ndarray):
                The indices of the first period of interest.
            ix2 (ndarrays):
                The indices of the second period of interest.
            h (ndarray):
                The separation distance between two sites (units of km).
                The result is placed in h.

        Returns:
            h (ndarray): The predicted correlation coefficient. The output
            array will have the same shape as the inputs.

        """
        # Verify the validity of input arguments
        if np.any(h < 0):
            raise ValueError("The separation distance must be positive")
        if np.shape(ix1) != np.shape(ix2) or np.shape(ix1) != np.shape(h):
            raise ValueError("The input arguments must all have the same dimensions")

        p1 = self.periods[ix1]
        p2 = self.periods[ix2]

        rho = p1 / p2
        invix = rho > 1.0
        rho[invix] = 1.0 / rho[invix]
        ad = np.abs(h)
        rho = rho * np.exp(-ad / 10)
        h[:] = rho[:]
        return h
