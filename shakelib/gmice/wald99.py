# third party imports
import numpy as np

# stdlib imports
from openquake.hazardlib import imt
from openquake.hazardlib.imt import PGA, PGV
from shakelib.gmice.gmice import GMICE


class Wald99(GMICE):
    """
    Implements the ground motion intensity conversion equations (GMICE) of
    Wald et al. (1999). This module implements a simplified version in that
    it only uses one of PGV or PGA, and not a combination of the two (PGV
    for higher intensities, PGA for lower) as is recommended in the reference.

    References:
        Wald, D.J., V. Quitoriano, T.H. Heaton, and H. Kanamori (1999).
        Relationships between peak gorund acceleration, peak ground
        velocity, and Modified Mercalli Intensity in California,
        Earthquake Spectra, Volume15, No. 3, August 1999.
    """

    # -----------------------------------------------------------------------
    #
    # Imm = 3.66 * log10(PGA) - 1.66 for Imm >= V (e.g., log10(PGA) >= 1.82)
    # Imm = 2.20 * log10(PGA) + 1.0  for Imm < 5
    #
    # Imm = 3.47 * log10(PGV) + 2.35 for Imm >= V
    # Imm = 2.10 * log10(PGV) + 3.40  for Imm < 5
    #
    # -----------------------------------------------------------------------
    def __init__(self):
        super().__init__()
        self.min_max = (1.0, 10.0)
        self.name = "Wald et al. (1999)"
        self.scale = "scale_wald99.ps"
        self.__constants = {
            self._pga: {
                "C1": 3.66,
                "C2": -1.66,
                "C3": 2.20,
                "C4": 1.00,
                "T1": 1.82,
                "T2": 5.00,
                "SMMI": 1.08,
                "SPGM": 0.295,
            },
            self._pgv: {
                "C1": 3.47,
                "C2": 2.35,
                "C3": 2.10,
                "C4": 3.40,
                "T1": 0.76,
                "T2": 5.00,
                "SMMI": 0.98,
                "SPGM": 0.282,
            },
        }

        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGA, PGV])

        self.DEFINED_FOR_SA_PERIODS = set([])

    def getPreferredMI(self, df, dists=None, mag=None):
        if "PGA" not in df or "PGV" not in df:
            if "PGV" in df:
                oqimt = imt.from_string("PGV")
                return self.getMIfromGM(df["PGV"], oqimt, dists, mag)[0]
            elif "PGA" in df:
                oqimt = imt.from_string("PGA")
                return self.getMIfromGM(df["PGA"], oqimt, dists, mag)[0]
            else:
                return None
        oqimt = imt.from_string("PGA")
        mmi_pga = self.getMIfromGM(df["PGA"], oqimt, dists, mag)[0]
        ix_nan_pga = np.isnan(mmi_pga)
        oqimt = imt.from_string("PGV")
        mmi_pgv = self.getMIfromGM(df["PGV"], oqimt, dists, mag)[0]
        ix_nan_pgv = np.isnan(mmi_pgv)
        ix_nan = ix_nan_pga | ix_nan_pgv
        vscale = np.zeros_like(mmi_pga)
        vscale[~ix_nan] = (mmi_pga[~ix_nan] - 5) / 2
        vscale[vscale < 0] = 0
        vscale[vscale > 1] = 1
        ascale = 1 - vscale
        mmi = np.full_like(mmi_pga, np.nan)
        mmi[~ix_nan] = np.clip(
            mmi_pga[~ix_nan] * ascale[~ix_nan] + mmi_pgv[~ix_nan] * vscale[~ix_nan],
            1,
            10,
        )
        mmi[ix_nan_pgv] = mmi_pga[ix_nan_pgv]
        mmi[ix_nan_pga] = mmi_pgv[ix_nan_pga]
        ix_nan = np.isnan(mmi)
        mmi95 = np.full_like(mmi, False, dtype=bool)
        mmi95[~ix_nan] = mmi[~ix_nan] > 9.5
        mmi[mmi95] = 10.0
        return mmi

    def getMIfromGM(self, amps, imt, dists=None, mag=None):
        """
        Function to compute macroseismic intensity from ground-motion
        intensity. Supported ground-motion IMTs are PGA and PGV

        Args:
            amps (ndarray):
                Ground motion amplitude; natural log units; g for PGA and
                PSA, cm/s for PGV.
            imt (OpenQuake IMT):
                Type the input amps (must be one of PGA or PGV).
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`
            dists (ndarray):
                Numpy array of distances from rupture (km). This parameter
                is ignored by this GMICE.
            mag (float):
                Earthquake magnitude. This parameter is ignored by this
                GMICE.

        Returns:
            ndarray of Modified Mercalli Intensity and ndarray of
            dMMI / dln(amp) (i.e., the slope of the relationship at the
            point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c = self._getConsts(imt)
        ix_nan = np.isnan(amps)

        #
        # Convert (for accelerations) from ln(g) to cm/s^2
        # then take the log10
        #
        if imt == self._pga:
            units = 981.0
        else:
            units = 1.0
        #
        # Math: log10(981 * exp(amps)) = log10(981) + log10(exp(amps))
        # = log10(981) + amps * log10(e)
        # For PGV, just convert ln(amp) to log10(amp) by multiplying
        # by log10(e)
        #
        lamps = np.zeros_like(amps)
        lamps[~ix_nan] = np.log10(units) + amps[~ix_nan] * lfact
        mmi = np.full_like(amps, np.nan)
        dmmi_damp = np.zeros_like(amps)
        #
        # This is the upper segment of the bi-linear fit
        #
        idx = lamps >= c["T1"]
        mmi[idx] = c["C2"] + c["C1"] * lamps[idx]
        dmmi_damp[idx] = c["C1"] * lfact
        #
        # This is the lower segment of the bi-linear fit
        #
        idx = lamps < c["T1"]
        mmi[idx] = c["C4"] + c["C3"] * lamps[idx]
        dmmi_damp[idx] = c["C3"] * lfact

        mmi = np.clip(mmi, 1.0, 10.0)
        mmi[ix_nan] = np.nan
        return mmi, dmmi_damp

    def getGMfromMI(self, mmi, imt, dists=None, mag=None):
        """
        Function to tcompute ground-motion intensity from macroseismic
        intensity. Supported IMTs are PGA and PGV.

        Args:
            mmi (ndarray):
                Macroseismic intensity.
            imt (OpenQuake IMT):
                IMT of the requested ground-motions intensities (must be
                one of PGA or PGV).
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`
            dists (ndarray):
                Rupture distances (km) to the corresponding MMIs. This
                parameter is ignored by this GMICE.
            mag (float):
                Earthquake magnitude. This parameter is ignored by this
                GMICE.

        Returns:
            Ndarray of ground motion intensity in natural log of g for PGA
            and PSA, and natural log cm/s for PGV; ndarray of dln(amp) / dMMI
            (i.e., the slope of the relationship at the point in question).
        """  # noqa
        lfact = np.log10(np.e)
        c = self._getConsts(imt)
        mmi = mmi.copy()
        ix_nan = np.isnan(mmi)
        mmi[ix_nan] = 1.0

        pgm = np.zeros_like(mmi)
        dpgm_dmmi = np.zeros_like(mmi)

        #
        # Upper segment of bi-linear relationship
        #
        idx = mmi >= c["T2"]
        pgm[idx] = np.power(10, (mmi[idx] - c["C2"]) / c["C1"])
        dpgm_dmmi[idx] = 1.0 / (c["C1"] * lfact)
        #
        # Lower segment of bi-linear relationship
        #
        idx = mmi < c["T2"]
        pgm[idx] = np.power(10, (mmi[idx] - c["C4"]) / c["C3"])
        dpgm_dmmi[idx] = 1.0 / (c["C3"] * lfact)

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
            self._pga: self.__constants[self._pga]["SMMI"],
            self._pgv: self.__constants[self._pgv]["SMMI"],
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
            self._pga: lfact * self.__constants[self._pga]["SPGM"],
            self._pgv: lfact * self.__constants[self._pgv]["SPGM"],
        }

    def _getConsts(self, imt):
        """
        Helper function to get the constants.
        """

        if imt != self._pga and imt != self._pgv:
            raise ValueError("Invalid IMT " + str(imt))
        return self.__constants[imt]
