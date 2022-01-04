# Local imports
from shakelib.conversions.convert_imt import IMTConverter


class BommerAlarcon2006(IMTConverter):
    """
    Class for conversion between PGV (units of cm/s) and PSA05 (units of g)
    by Bommer and Alarcon (2006).

    - PSA05 stands for spectral acceleration with oscillator period of 0.5 sec
    - PGV is peak ground velocity.

    References:
        Bommer, J. J., & Alarcon, J. E. (2006). The prediction and use of peak
        ground velocity. Journal of Earthquake Engineering, 10(01), 1-31.
        `[link] <http://www.worldscientific.com/doi/abs/10.1142/S1363246906002463>`__
    """

    def __init__(self):
        super().__init__()
        # output_input dictionary where the key is the output
        # and the value is a list of the possible inputs
        self.output_input = {"PGV": ["PSA05"], "PSA05": ["PGV"]}
        self.conversion_factor = 1.0 / (20.0) * 100.0 * 9.81

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
        # Verify that this is a valid conversion
        self._verifyConversion(imt_in, imt_out)
        imt_in = imt_in.upper().strip()
        imt_out = imt_out.upper().strip()
        conversion_factor = self.conversion_factor
        # Check which method to use
        if imt_in == "PSA05" and imt_out == "PGV":
            new_imt = self._convertToPGV(imt, conversion_factor)
        elif imt_in == "PGV" and imt_out == "PSA05":
            new_imt = self._convertToPSA05(imt, conversion_factor)
        else:
            raise ValueError(f"No conversion available from {imt_in!r} to {imt_out!r}")
        return new_imt

    @staticmethod
    def _convertToPGV(psa05, conversion_factor):
        """
        Convert PSA05 (spectral acceleration with oscillator period of 0.5 sec)
        in g to PGV cm/s.
        **Important:** PSA10 must be linear units.

        Args:
            psa05 (array): Numpy array or float of PSA05 values; linear units.

        Returns
            array: Numpy array or float of PGV converted from psa05.
        """
        return psa05 * conversion_factor

    @staticmethod
    def _convertToPSA05(pgv, conversion_factor):
        """
        Convert PGV in cm/s to PSA05 in g.
        **Important:** PGV must be linear units.

        Args:
            pgv (array): Numpy array or float of PGV values; linear units.

        Returns:
            array: Numpy array or float of PSA05 (spectral acceleration with
            oscillator period of 0.5 sec) converted from PGV.
        """
        return pgv / conversion_factor
