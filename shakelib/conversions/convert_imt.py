# Standard library imports
from abc import ABC, abstractmethod
import logging


class IMTConverter(ABC):
    """Base class for implementing conversions between intensity
    measurement types (IMT)."""

    @abstractmethod
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
        pass

    def getConversionFactor(self):
        """
        Helper method that returns the conversion factor.

        Returns:
            float: Conversion factor.
        """
        factor = self.conversion_factor
        return factor

    def getInputIMT(self, imt_out):
        """
        Get valid input IMT types that can be converted to the specified
        imt_out.

        Args:
            imt_out (str): OQ intensity measure type.

        Returns:
            list: List of valid input IMT types. If not available types
                None is returned.
        """
        imt_out = imt_out.upper().strip()
        try:
            imt_in = self.output_input[imt_out]
            return imt_in
        except KeyError:
            logging.info(f"No available conversion to {imt_out!r}")
            return None

    def _verifyConversion(self, imt_in, imt_out):
        """
        Helper method used to verify that a conversion is valid.

        Args:
            imt_in (str): OQ intensity measure type. Same as type as the input
                values defined by the imt variable.
            imt_out (str): OQ intensity measure type that the values will
                be converted to.

        Raises:
            ValueError: If the conversion is not valid.
        """
        valid_inputs = self.getInputIMT(imt_out)
        imt_in = imt_in.upper().strip()
        imt_out = imt_out.upper().strip()
        if imt_in not in valid_inputs and imt_in != imt_out:
            raise ValueError(f"No conversion available from {imt_in!r} to {imt_out!r}")
