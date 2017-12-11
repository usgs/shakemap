
import copy
import numpy as np


class NewmarkHall1982(object):

    """
    Class for conversion between PGA and PSA10 by Newmark and Hall (1982).

    - PSA10 stands for spectral acceleration with oscillator period of 1.0 sec
    - PGV is peak ground velocity.

    The conversion factor is 37.27*2.54.
    Note that 2.54 is the conversion factor to convert from cm/s to in/s and
    37.27 is derived from SA(f) = 1.65*(2*pi*V*f)/386.09, where:

    - SA(f) is spectral acceleration at frequency f (in g)
    - f is the frequency of interest
    - V is the velocity
    - 1.65 is the N&H amplification factor for velocity at 5% damping
    - 386.09 is the acceleration of gravity in inches per second per g

    The sigma factor was computed from an average sigma value determined
    by plotting PGV/PSA10 versus Distance for earthquakes with magnitudes
    greater than or equal to 5.0.

    To do
        - Inherit from ConvertIMT class.

    References:
        Newmark, N. M., & Hall, W. J. (1982). Earthquake spectra and design.
        Earthquake Engineering Research Institute, El Cerrito, California.
    """

    __lnfact = np.log(37.27 * 2.54)
    __lnsigma = 0.5146578

    @staticmethod
    def psa102pgv(psa10, sigma):
        """
        Convert PSA10 (spectral acceleration with oscillator period of 1.0 sec)
        to PGV.

        **Important:** PSA10 and sigma must be logarithmic units.

        :param psa10:
            Numpy array or float of PSA10 values; units must be natural log
            of g.
        :param sigma:
            Numpy array or float of standard deviation of PGV from a GMPE;
            units must be natural log.
        :returns:
            Two arrays
                - Array of PGV iwth natural log of cm/s units.
                - Array of its standard deviation with natural log units.
        """
        pgv = psa10 + NewmarkHall1982.__lnfact

        sigmaTot = np.sqrt((sigma ** 2) +
                           (NewmarkHall1982.__lnsigma ** 2))

        return pgv, sigmaTot

    @staticmethod
    def getConversionFactor():
        """
        :returns:
            The Newmark and Hall (1982) multiplicative conversion factor for
            convering Sa(1.0) to PGV (cm/s).
        """
        return np.exp(NewmarkHall1982.__lnfact)

    @staticmethod
    def getLnSigma():
        """
        :returns:
            The the estimate of the logarithmic standard deviation of the PGV
            predicted by Newmark and Hall (1982).
        """
        return copy.copy(NewmarkHall1982.__lnsigma)
