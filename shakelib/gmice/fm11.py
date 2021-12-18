# third party imports
import numpy as np

# stdlib imports
from openquake.hazardlib.imt import PGA, PGV, SA
from shakelib.gmice.gmice import GMICE


class FM11(GMICE):
    """
    Implements the ground motion intensity conversion equations (GMICE) of
    Faenza and Michelini (2010, 2011).

    References:
        Faenza and Michelini,(2010). Regression analysis of MCS intensity and
        ground motion parameters in Italy and its application in ShakeMap.
        GJI, 180, 1138-1152, doi: 10.1111/j.1365-246X.2009.04467.xself
        and Faenza and Michelini (2011). Regression analysis of MCS intensity
        and ground motion spectral accelerations (SAs) in Italy. GJI, 186,
        1415-1430, doi: 10.1111/j.1365-246X.2011.05125.x
    """

    # -----------------------------------------------------------------------
    # MMI = C1 + C2 * log10 (Y)
    #
    # Limit the distance residuals to between 0 and 160 km.
    # Limit the magnitude residuals to between M3.0 and M6.9.
    # These are the default in the other gmice
    #      Limit the distance residuals to between 10 and 300 km.
    #      Limit the magnitude residuals to between M3.0 and M7.3.
    #
    # These are calcualted on the basis of the maximum horizontal component.
    # For psa 03, 10 and 30 the regression for the geometrical mean are
    # available but not implemented in this modeles since the one for PGA
    # and PGV are not available.
    # -----------------------------------------------------------------------
    def __init__(self):
        super().__init__()
        self.min_max = (1.0, 10.0)
        self.name = "Faenza and Michelini (2010, 2011)"
        self.scale = "scale_fm11.ps"
        self._constants = {
            self._pga: {"C1": 1.68, "C2": 2.58, "SMMI": 0.18, "SPGM": 0.31},
            self._pgv: {"C1": 5.11, "C2": 2.35, "SMMI": 0.14, "SPGM": 0.22},
            self._sa03: {"C1": 1.24, "C2": 2.47, "SMMI": 0.30, "SPGM": 0.42},
            self._sa10: {"C1": 3.12, "C2": 2.05, "SMMI": 0.21, "SPGM": 0.31},
            self._sa30: {"C1": 4.31, "C2": 2.00, "SMMI": 0.14, "SPGM": 0.26},
        }

        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGA, PGV, SA])

        self.DEFINED_FOR_SA_PERIODS = set([0.3, 1.0, 3.0])

    def getMIfromGM(self, amps, imt, dists=None, mag=None):
        """
        Function to compute macroseismic intensity from ground-motion
        intensity. Supported ground-motion IMTs are PGA, PGV and PSA
        at 0.3, 1.0, and 3.0 sec periods.

        Args:
            amps (ndarray):
                Ground motion amplitude; natural log units; g for PGA and
                PSA, cm/s for PGV.
            imt (OpenQuake IMT):
                Type the input amps (must be one of PGA, PGV, or SA).
                Supported SA periods are 0.3, 1.0, and 3.0 sec.
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`
            dists (ndarray):
                Not used
            mag (float):
                Not used

        Returns:
            ndarray of Modified Mercalli Intensity and ndarray of
            dMMI / dln(amp) (i.e., the slope of the relationship at the
            point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c = self._getConsts(imt)

        #
        # Convert (for accelerations) from ln(g) to cm/s^2
        # then take the log10
        #
        if imt != self._pgv:
            units = 981.0
        else:
            units = 1.0
        #
        # Math: log10(981 * exp(amps)) = log10(981) + log10(exp(amps))
        # = log10(981) + amps * log10(e)
        # For PGV, just convert ln(amp) to log10(amp) by multiplying
        # by log10(e)
        #
        lamps = np.log10(units) + amps * lfact

        mmi = c["C1"] + c["C2"] * lamps
        dmmi_damp = np.full_like(lamps, c["C2"] * lfact)

        mmi = np.clip(mmi, 1.0, 10.0)
        mmi[np.isnan(amps)] = np.nan
        return mmi, dmmi_damp

    def getGMfromMI(self, mmi, imt, dists=None, mag=None):
        """
        Function to tcompute ground-motion intensity from macroseismic
        intensity. Supported IMTs are PGA, PGV and PSA for 0.3, 1.0, and
        3.0 sec periods.

        Args:
            mmi (ndarray):
                Macroseismic intensity.
            imt (OpenQuake IMT):
                IMT of the requested ground-motions intensities (must be
                one of PGA, PGV, or SA).
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`
            dists (ndarray):
                Not used
            mag (float):
                Not used

        Returns:
            Ndarray of ground motion intensity in natural log of g for PGA
            and PSA, and natural log cm/s for PGV; ndarray of dln(amp) / dMMI
            (i.e., the slope of the relationship at the point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c = self._getConsts(imt)
        mmi = mmi.copy()
        # Set nan values to 1
        ix_nan = np.isnan(mmi)
        mmi[ix_nan] = 1.0

        pgm = np.zeros_like(mmi)
        dpgm_dmmi = np.zeros_like(mmi)
        dummy_variable = np.ones(len(mmi))

        #
        # MMI to PGM
        #
        pgm = np.power(10, (mmi - c["C1"]) / c["C2"])
        dpgm_dmmi = 1.0 / (c["C2"] * lfact) * dummy_variable

        if imt != self._pgv:
            units = 981.0
        else:
            units = 1.0

        # Return a ln(amp) value. Convert PGA to from cm/s^2 to g
        pgm /= units
        pgm = np.log(pgm)

        # Set nan values back from 1 to nan
        pgm[ix_nan] = np.nan
        dpgm_dmmi[ix_nan] = np.nan

        return pgm, dpgm_dmmi

    def getGM2MIsd(self):
        """
        Return a dictionary of standard deviations for the ground-motion
        to MMI conversion. The keys are the ground motion types.

        Returns:
            Dictionary of GM to MI sigmas (in MMI units).
        """
        return {
            self._pga: self._constants[self._pga]["SMMI"],
            self._pgv: self._constants[self._pgv]["SMMI"],
            self._sa03: self._constants[self._sa03]["SMMI"],
            self._sa10: self._constants[self._sa10]["SMMI"],
            self._sa30: self._constants[self._sa30]["SMMI"],
        }

    def getMI2GMsd(self):
        """
        Return a dictionary of standard deviations for the MMI
        to ground-motion conversion. The keys are the ground motion
        types.

        Returns:
            Dictionary of MI to GM sigmas (ln(PGM) units).
        """
        #
        # Need to convert log10 to ln units
        #
        lfact = np.log(10.0)
        return {
            self._pga: lfact * self._constants[self._pga]["SPGM"],
            self._pgv: lfact * self._constants[self._pgv]["SPGM"],
            self._sa03: lfact * self._constants[self._sa03]["SPGM"],
            self._sa10: lfact * self._constants[self._sa10]["SPGM"],
            self._sa30: lfact * self._constants[self._sa30]["SPGM"],
        }

    def _getConsts(self, imt):
        """
        Helper function to get the constants.
        """

        if (
            imt != self._pga
            and imt != self._pgv
            and imt != self._sa03
            and imt != self._sa10
            and imt != self._sa30
        ):
            raise ValueError("Invalid IMT " + str(imt))
        c = self._constants[imt]
        return c
