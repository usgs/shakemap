# stdlib imports
from abc import ABC, abstractmethod

# third party imports
from openquake.hazardlib.imt import PGA, PGV, SA, from_string

# local imports
from openquake.hazardlib import imt
from openquake.hazardlib.imt import PGA, PGV, SA


class GMICE(ABC):
    """
    Base class called to implement ground motion intensity conversion
    equations (GMICE).

    Inherited by AK07, Wald99, WGRW12.
    """

    def __init__(self):
        self._pga = PGA()
        self._pgv = PGV()
        self._sa03 = SA(0.3)
        self._sa10 = SA(1.0)
        self._sa30 = SA(3.0)
        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set()

    @staticmethod
    def getDistanceType():
        return "rrup"

    def supports(self, imt):
        """Determine whether input IMT is supported by GMICE instance.

        Args:
            imt (str): Valid IMT string - 'MMI', 'PGV', 'PGA', 'SA(0.3)', etc.
        Returns:
            bool: True if gmice is defined for input IMT (and period), False
                  if not.
        """
        for imtcomp in self.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
            thisimt = from_string(imt)
            if thisimt.string == imtcomp.__name__:
                return True
            if "SA" in thisimt.string:
                for period in self.DEFINED_FOR_SA_PERIODS:
                    if period == thisimt.period:
                        return True
                return False
        return False

    def getMinMax(self):
        """
        Get the minimum and maximum MMI values produced by this GMICE.

        Returns:
            Tuple of min and max values of GMICE.
        """
        return self.min_max

    def getName(self):
        """
        Get the name of this GMICE.

        Returns:
            String containing name of GMICE.
        """
        return self.name

    def getScale(self):
        """
        Get the name of the PostScript file containing this GMICE's
        scale.

        Returns:
            Name of GMICE scale file.
        """
        return self.scale

    def getPreferredMI(self, df, dists=None, mag=None):
        """
        Function to compute macroseismic intensity from the preferred
        ground-motion intensity. The function uses PGV by default, but
        this may be overridden by individual classes.

        Args:
            df (dict):
                Dictionaary containing all of the available ground
                motions.
            dists (ndarray):
                Numpy array of distances from rupture (km).
            mag (float):
                Earthquake magnitude.

        Returns:
            ndarray of Modified Mercalli Intensity and ndarray of
            dMMI / dln(amp) (i.e., the slope of the relationship at the
            point in question).
        """
        if "PGV" not in df:
            return None
        oqimt = imt.from_string("PGV")
        return self.getMIfromGM(df["PGV"], oqimt, dists, mag)[0]

    def getPreferredSD(self):
        """
        Return an array of standard deviations for the preferred
        ground-motion to MMI conversion. Return None is the preferred
        IMT is not in the list of available IMTs for the GMICE.

        Returns:
            (numpy array):  Array of GM to MI sigmas (in MMI units).
        """
        oqimt = imt.from_string("PGV")
        sd = self.getGM2MIsd()
        if oqimt in sd:
            return sd[oqimt]
        else:
            return None

    @abstractmethod
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
        """  # noqa
        pass

    @abstractmethod
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
        pass

    @abstractmethod
    def getGM2MIsd(self):
        """
        Return a dictionary of standard deviations for the ground-motion
        to MMI conversion. The keys are the ground motion types.

        Returns:
            Dictionary of GM to MI sigmas (in MMI units).
        """
        pass

    @abstractmethod
    def getMI2GMsd(self):
        """
        Return a dictionary of standard deviations for the MMI
        to ground-motion conversion. The keys are the ground motion
        types.

        Returns:
            Dictionary of MI to GM sigmas (ln(PGM) units).
        """
        pass

    @abstractmethod
    def _getConsts(self, imt):
        """
        Helper function to get the constants.
        """
        pass
