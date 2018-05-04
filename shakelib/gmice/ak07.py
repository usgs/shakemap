import numpy as np

from openquake.hazardlib.imt import PGA, PGV, SA


class AK07(object):
    """
    Implements the ground motion intensity conversion equations (GMICE) of
    Atkinson and Kaka (2007). This module implements a simplified version.

    References:
        Atkinson, Gail & Kaka, SanLinn. (2007). Relationships between
        Felt Intensity and Instrumental Ground Motion in the Central
        United States and California. Bulletin of The Seismological
        Society of America - BULL SEISMOL SOC AMER.
        97. 497-510. 10.1785/0120060154.
    """

    # -----------------------------------------------------------------------
    #
    # MMI = C1 + C2 * log(Y)          for c2->T1 < log(Y) <= T1
    # MMI = C3 + C4 * log(Y)          for log(Y) > T1
    #
    # or
    #
    # MMI = c2->C1 + c2->C2 * log(Y) + C5 + C6 * M + C7 * log(D)
    #                            for log(Y) <= c2->T1
    # MMI = C1 + C2 * log(Y) + C5 + C6 * log(D) + C7 * M
    #                            for c2->T1 < log(Y) <= T1
    # MMI = C3 + C4 * log(Y) + C5 + C6 * M + C7 * log(D)
    #                           for log(Y) > T1
    #
    # Limit the distance residuals to between 10 and 300 km.
    # Limit the magnitude residuals to between M3.0 and M7.3.
    #
    # Constants taken from Table 4 and 5 in the reference material
    #
    # -----------------------------------------------------------------------

    __pga = PGA()
    __pgv = PGV()
    __sa03 = SA(0.3)
    __sa10 = SA(1.0)
    __sa30 = SA(3.0)
    __constants = {
        __pga: {'C1':  2.65, 'C2':  1.39, 'C3':  -1.91, 'C4': 4.09,
                'C5':  -1.96, 'C6':  0.02, 'C7': 0.98, 'T1':  1.69,
                'T2': 5, 'SMMI': 0.89},
        __pgv: {'C1':  4.37, 'C2':  1.32, 'C3': 3.54, 'C4': 3.03,
                'C5': 0.47, 'C6':  -0.19, 'C7': 0.26, 'T1': 0.48,
                'T2': 5, 'SMMI': 0.76},
        __sa03: {'C1':  2.4, 'C2':  1.36, 'C3': -1.83, 'C4': 3.56,
                 'C5': -0.11, 'C6':  -0.2, 'C7':  0.64, 'T1':  1.92,
                 'T2': 5, 'SMMI': 0.79},
        __sa10: {'C1':  3.23, 'C2':  1.18, 'C3':  0.57, 'C4': 2.95,
                 'C5':  1.92, 'C6': -0.39, 'C7': 0.04, 'T1':  1.5,
                 'T2': 5, 'SMMI': 0.73},
        __sa30: {'C1':  3.72, 'C2':  1.29, 'C3':  1.99, 'C4': 3.0,
                 'C5':  2.24, 'C6': -0.33, 'C7': -0.31, 'T1':  1,
                 'T2': 5, 'SMMI': 0.72}
    }

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        PGA,
        PGV,
        SA
    ])

    DEFINED_FOR_SA_PERIODS = set([0.3, 1.0, 3.0])

    def getMIfromGM(self, amps, imt, dists=None, mag=None):
        """
        Function to compute macroseismic intensity from ground-motion
        intensity. Supported ground-motion IMTs are PGA, PGV and PSA
        at 0.3, 1.0, and 2.0 sec periods.

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
        """
        lfact = np.log10(np.e)
        c = self.__getConsts(imt)
        if dists is not None and mag is not None:
            doresid = True
            # Limit distances and take the log10
            ldd = np.log10(np.clip(dists, 10, 300))
            # Llimit magnitudes
            lmm = np.clip(mag, 3.0, 7.3)
        else:
            doresid = False

        # Account for ln form of amplitudes
        # Convert (for accelerations) from ln(g) to cm/s^2
        # then take the log10
        #
        if imt != self.__pgv:
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
        # This is the lower segment of the bi-linear fit
        # where log(Y) is less than T1
        #
        idx = (lamps < c['T1'])
        mmi[idx] = c['C1'] + c['C2'] * lamps[idx]
        dmmi_damp[idx] = c['C2'] * lfact
        #
        # This is the upper segment of the bi-linear fit
        # where log(Y) is greater than or equal to T1
        #
        idx = lamps >= c['T1']
        mmi[idx] = c['C3'] + c['C4'] * lamps[idx]
        dmmi_damp[idx] = c['C4'] * lfact

        # Inclusion of residuals if magnitude and
        # distance information is available
        if doresid:
            mmi += c['C5'] + c['C6'] * lmm + c['C7'] * ldd

        # Limit mmi values
        mmi = np.clip(mmi, 1.0, 10.0)
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
        """
        lfact = np.log10(np.e)
        c = self.__getConsts(imt)
        mmi = mmi.copy()
        # Set nan values to 1
        ix_nan = np.isnan(mmi)
        mmi[ix_nan] = 1.0

        if dists is not None and mag is not None:
            doresid = True
            # Limit distances and take the log10
            ldd = np.log10(np.clip(dists, 10, 300))
            # Limit magnitudes
            lmm = np.clip(mag, 3.0, 7.3)
        else:
            doresid = False

        # Inclusion of residuals if magnitude and
        # distance information is available
        if doresid:
            mmi -= c['C5'] + c['C6'] * lmm + c['C7'] * ldd

        pgm = np.zeros_like(mmi)
        dpgm_dmmi = np.zeros_like(mmi)
        #
        # This is the lower segment of the bi-linear fit
        # where MMI is less than I5
        #
        idx = mmi < c['T2']
        pgm[idx] = np.power(10, (mmi[idx] - c['C1']) / c['C2'])
        dpgm_dmmi[idx] = 1.0 / (c['C2'] * lfact)
        #
        # This is the upper segment of the bi-linear fit
        # where MMI is greater than or equal to I5
        #
        idx = mmi >= c['T2']
        pgm[idx] = np.power(10, (mmi[idx] - c['C3']) / c['C4'])
        dpgm_dmmi[idx] = 1.0 / (c['C4'] * lfact)

        if imt != self.__pgv:
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
        return {self.__pga: self.__constants[self.__pga]['SMMI'],
                self.__pgv: self.__constants[self.__pgv]['SMMI'],
                self.__sa03: self.__constants[self.__sa03]['SMMI'],
                self.__sa10: self.__constants[self.__sa10]['SMMI'],
                self.__sa30: self.__constants[self.__sa30]['SMMI']}

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
        return {self.__pga: lfact * 0.57,
                self.__pgv: lfact * 0.52,
                self.__sa03: lfact * 0.63,
                self.__sa10: lfact * 57,
                self.__sa30: lfact * 81}

    @staticmethod
    def getName():
        """
        Get the name of this GMICE.

        Returns:
            String containing name of this GMICE.
        """
        return 'Atkinson and Kaka (2007)'

    @staticmethod
    def getScale():
        """
        Get the name of the PostScript file containing this GMICE's
        scale.

        Returns:
            Name of GMICE scale file.
        """
        return 'scale_ak07.ps'

    @staticmethod
    def getMinMax():
        """
        Get the minimum and maximum MMI values produced by this GMICE.

        Returns:
            Tuple of min and max values of GMICE.
        """
        return (1.0, 10.0)

    @staticmethod
    def getDistanceType():
        return 'rrup'

    def __getConsts(self, imt):
        """
        Helper function to get the constants.
        """

        if (imt != self.__pga and imt != self.__pgv and imt != self.__sa03 and
                imt != self.__sa10 and imt != self.__sa30):
            raise ValueError("Invalid IMT " + str(imt))
        c = self.__constants[imt]
        return c
