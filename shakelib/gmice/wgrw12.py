# third party imports
import numpy as np

# stdlib imports
from openquake.hazardlib.imt import PGA, PGV, SA
from shakelib.gmice.gmice import GMICE


class WGRW12(GMICE):
    """
    Implements the ground motion intensity conversion equations (GMICE) of
    Worden et al. (2012).

    References:
        Worden, C. B., Gerstenberger, M. C., Rhoades, D. A., & Wald, D. J.
        (2012). Probabilistic relationships between ground‐motion parameters
        and modified Mercalli intensity in California. Bulletin of the
        Seismological Society of America, 102(1), 204-221.
    """

    # -----------------------------------------------------------------------
    #
    # MMI = c2->C1 + c2->C2 * log(Y)  for log(Y) <= c2->T1
    # MMI = C1 + C2 * log(Y)          for c2->T1 < log(Y) <= T1
    # MMI = C3 + C4 * log(Y)          for log(Y) > T1
    #
    # or
    #
    # MMI = c2->C1 + c2->C2 * log(Y) + C5 + C6 * log(D) + C7 * M
    #                            for log(Y) <= c2->T1
    # MMI = C1 + C2 * log(Y) + C5 + C6 * log(D) + C7 * M
    #                            for c2->T1 < log(Y) <= T1
    # MMI = C3 + C4 * log(Y) + C5 + C6 * log(D) + C7 * M for log(Y) > T1
    #
    # Limit the distance residuals to between 10 and 300 km.
    # Limit the magnitude residuals to between M3.0 and M7.3.
    #
    # -----------------------------------------------------------------------
    def __init__(self):
        super().__init__()
        self.min_max = (1.0, 10.0)
        self.name = "Worden et al. (2012)"
        self.scale = "scale_wgrw12.ps"
        self._constants = {
            self._pga: {
                "C1": 1.78,
                "C2": 1.55,
                "C3": -1.60,
                "C4": 3.70,
                "C5": -0.91,
                "C6": 1.02,
                "C7": -0.17,
                "T1": 1.57,
                "T2": 4.22,
                "SMMI": 0.66,
                "SPGM": 0.35,
            },
            self._pgv: {
                "C1": 3.78,
                "C2": 1.47,
                "C3": 2.89,
                "C4": 3.16,
                "C5": 0.90,
                "C6": 0.00,
                "C7": -0.18,
                "T1": 0.53,
                "T2": 4.56,
                "SMMI": 0.63,
                "SPGM": 0.38,
            },
            self._sa03: {
                "C1": 1.26,
                "C2": 1.69,
                "C3": -4.15,
                "C4": 4.14,
                "C5": -1.05,
                "C6": 0.60,
                "C7": 0.00,
                "T1": 2.21,
                "T2": 4.99,
                "SMMI": 0.82,
                "SPGM": 0.44,
            },
            self._sa10: {
                "C1": 2.50,
                "C2": 1.51,
                "C3": 0.20,
                "C4": 2.90,
                "C5": 2.27,
                "C6": -0.49,
                "C7": -0.29,
                "T1": 1.65,
                "T2": 4.98,
                "SMMI": 0.75,
                "SPGM": 0.47,
            },
            self._sa30: {
                "C1": 3.81,
                "C2": 1.17,
                "C3": 1.99,
                "C4": 3.01,
                "C5": 1.91,
                "C6": -0.57,
                "C7": -0.21,
                "T1": 0.99,
                "T2": 4.96,
                "SMMI": 0.89,
                "SPGM": 0.64,
            },
        }

        self._constants2 = {
            self._pga: {"C1": 1.71, "C2": 2.08, "T1": 0.14, "T2": 2.0},
            self._pgv: {"C1": 4.62, "C2": 2.17, "T1": -1.21, "T2": 2.0},
            self._sa03: {"C1": 1.15, "C2": 1.92, "T1": 0.44, "T2": 2.0},
            self._sa10: {"C1": 2.71, "C2": 2.17, "T1": -0.33, "T2": 2.0},
            self._sa30: {"C1": 7.35, "C2": 3.45, "T1": -1.55, "T2": 2.0},
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
                Numpy array of distances from rupture (km).
            mag (float):
                Earthquake magnitude.

        Returns:
            ndarray of Modified Mercalli Intensity and ndarray of
            dMMI / dln(amp) (i.e., the slope of the relationship at the
            point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c, c2 = self._getConsts(imt)

        if dists is not None and mag is not None:
            doresid = True
            ldd = np.log10(np.clip(dists, 10, 300))
            lmm = np.clip(mag, 3.0, 7.3)
        else:
            doresid = False

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
        mmi = np.zeros_like(amps)
        dmmi_damp = np.zeros_like(amps)
        #
        # This is the MMI 1 to 2 range that is discussed in the paper but not
        # specifically implemented there
        #
        idx = lamps < c2["T1"]
        mmi[idx] = c2["C1"] + c2["C2"] * lamps[idx]
        dmmi_damp[idx] = c2["C2"] * lfact
        #
        # This is the lower segment of the bi-linear fit
        #
        idx = (lamps >= c2["T1"]) & (lamps < c["T1"])
        mmi[idx] = c["C1"] + c["C2"] * lamps[idx]
        dmmi_damp[idx] = c["C2"] * lfact
        #
        # This is the upper segment of the bi-linear fit
        #
        idx = lamps >= c["T1"]
        mmi[idx] = c["C3"] + c["C4"] * lamps[idx]
        dmmi_damp[idx] = c["C4"] * lfact

        if doresid:
            mmi += c["C5"] + c["C6"] * ldd + c["C7"] * lmm

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
                Rupture distances (km) to the corresponding MMIs.
            mag (float):
                Earthquake magnitude.

        Returns:
            Ndarray of ground motion intensity in natural log of g for PGA
            and PSA, and natural log cm/s for PGV; ndarray of dln(amp) / dMMI
            (i.e., the slope of the relationship at the point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c, c2 = self._getConsts(imt)
        mmi = mmi.copy()
        ix_nan = np.isnan(mmi)
        mmi[ix_nan] = 1.0

        if dists is not None and mag is not None:
            doresid = True
            ldd = np.log10(np.clip(dists, 10, 300))
            lmm = np.clip(mag, 3.0, 7.3)
        else:
            doresid = False

        if doresid:
            mmi -= c["C5"] + c["C6"] * ldd + c["C7"] * lmm

        pgm = np.zeros_like(mmi)
        dpgm_dmmi = np.zeros_like(mmi)

        #
        # MMI 1 to 2
        #
        idx = mmi < 2.0
        pgm[idx] = np.power(10, (mmi[idx] - c2["C1"]) / c2["C2"])
        dpgm_dmmi[idx] = 1.0 / (c2["C2"] * lfact)
        #
        # Lower segment of bi-linear relationship
        #
        idx = (mmi >= 2.0) & (mmi < c["T2"])
        pgm[idx] = np.power(10, (mmi[idx] - c["C1"]) / c["C2"])
        dpgm_dmmi[idx] = 1.0 / (c["C2"] * lfact)
        #
        # Upper segment of bi-linear relationship
        #
        idx = mmi >= c["T2"]
        pgm[idx] = np.power(10, (mmi[idx] - c["C3"]) / c["C4"])
        dpgm_dmmi[idx] = 1.0 / (c["C4"] * lfact)

        if imt != self._pgv:
            units = 981.0
        else:
            units = 1.0
        pgm /= units
        pgm = np.log(pgm)
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
        c2 = self._constants2[imt]
        return (c, c2)
