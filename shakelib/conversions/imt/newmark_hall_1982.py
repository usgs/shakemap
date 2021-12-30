# Standard library imports
import copy

# Third party imports
import numpy as np

# Local imports
from shakelib.conversions.convert_imt import IMTConverter


class NewmarkHall1982(IMTConverter):
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

    References:
        Newmark, N. M., & Hall, W. J. (1982). Earthquake spectra and design.
        Earthquake Engineering Research Institute, El Cerrito, California.
    """

    def __init__(self):
        super().__init__()
        # output_input dictionary where the key is the output
        # and the value is a list of the possible inputs
        self.output_input = {"PGV": ["PSA10"]}
        self._lnfact = np.log(37.27 * 2.54)
        self.conversion_factor = np.exp(self._lnfact)
        self._lnsigma = 0.5146578

    def convertAmps(self, imt_in, imt_out, imt):
        """
        Returns an array of converted IMT amplitude values.

        Args:
            imt_in (str): OQ intensity measure type. Same as type as the input
                values defined by the imt variable.
            imt_out (str): OQ intensity measure type that the values will
                be converted to.
            imt (OpenQuake IMT): The intensity measurements of the input
                ground motions. Valid IMTs are PGV, and SA.

        Returns:
            array: Numpy array of amps converted from imt_in to imt_out.

        Raises:
            ValueError: If not a valid conversion.
        """
        self._verifyConversion(imt_in, imt_out)
        imt_in = imt_in.upper().strip()
        imt_out = imt_out.upper().strip()
        conversion_factor = self._lnfact
        if imt_in == "PSA10" and imt_out == "PGV":
            new_imt = self._convertToPGV(imt, conversion_factor)
        else:
            raise ValueError(f"No conversion available from {imt_in!r} to {imt_out!r}")
        return new_imt

    def convertSigmas(self, imt_in, imt_out, sigma):
        """
        Convert standard deviation to natural log units.

        Args:
            imt_in (str): OQ intensity measure type. Same as type as the input
                values defined by the imt variable.
            imt_out (str): OQ intensity measure type that the values will
                be converted to.
            sigma (array): Numpy array or float of standard deviation of PGV from a GMPE;
                units must be natural log.

        Returns:
            array: Converted standard deviations with natural log units.
        """
        self._verifyConversion(imt_in, imt_out)
        lnsigma = self._lnsigma
        sigmaTot = np.sqrt((sigma ** 2) + (lnsigma ** 2))
        return sigmaTot

    def getLnSigma(self):
        """
        Returns:
            float: The the estimate of the logarithmic standard deviation of the PGV
                predicted by Newmark and Hall (1982).
        """
        return copy.copy(self._lnsigma)

    @staticmethod
    def _convertToPGV(psa10, lnfact):
        """
        Convert PSA10 (spectral acceleration with oscillator period of 1.0 sec)
        to PGV.

        **Important:** PSA10 and sigma must be logarithmic units.

        Args:
            psa10 (array): OQ intensity measure type. Same as type as the input
                values defined by the imt variable.
            lnfact (float): Conversion factor.

        Returns:
            array: Converted PGV values with natural log of cm/s units.
        """
        pgv = psa10 + lnfact
        return pgv
