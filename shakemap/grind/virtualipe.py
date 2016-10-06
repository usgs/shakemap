import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.imt import PGA, PGV, SA, MMI

from shakemap.utils.exception import ShakeMapException

class VirtualIPE(GMPE):
    """
    Implements a virtual IPE that is the combination of a MultiGMPE
    and a GMICE. Will first attempt to use PGV for the intensities,
    and then will try PGA, and then will bail out.

    """

    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([MMI])
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = None
    REQUIRES_DISTANCES = None
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    REQUIRES_RUPTURE_PARAMETERS = None
    REQUIRES_SITES_PARAMETERS = None

    @classmethod
    def fromFuncs(cls, gmpe, gmice):
        """
        Args:
            gmpe: An instance of the MultiGMPE object.
            gmice: An instance of a GMICE object.

        Returns:
            A new instance of a VirtualIPE object.
        """
        self = cls()
        self.gmpe = gmpe
        self.gmice = gmice

        if (PGV in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES and
            PGV in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
            self.imt = PGV()
        elif (PGA in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES and
              PGA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
            self.imt = PGA()
        elif (SA in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES and
              SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
                 self.imt = SA(1.0)
        else:
            raise ShakeMapException(
                    'The supplied GMPE and GMICE do not use a common IMT'
                )
        
        return self

    def get_mean_and_stddevs(self, sx, rx, dx, imt, stddev_types):
        """
        See superclass `method <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.GroundShakingIntensityModel.get_mean_and_stddevs>`__. 
        """

        if imt != MMI():
            raise ValueError("imt must be MMI")
        #
        # Get the mean ground motions and stddev for the preferred IMT
        #
        mgm, sdev = self.gmpe.get_mean_and_stddevs(sx, rx, dx, self.imt, 
                stddev_types)

        #
        # Get the MMI and the dMMI/dPGM from the GMICE
        #
        if hasattr(dx, 'rrup'):
            dist4gmice = dx.rrup
        else:
            dist4gmice = dx.rhypo

        mmi, dmda = self.gmice.getMIfromGM(mgm, self.imt, dist4gmice, rx.mag)

        #
        # Compute the uncertainty of the combined prediction/conversion
        #
        gm2mi_sd = self.gmice.getGM2MIsd()[self.imt]
        gm_var_in_mmi = dmda**2 * sdev[0]**2
        mmi_sd = np.sqrt(gm2mi_sd**2 + gm_var_in_mmi)
        
        return mmi, mmi_sd
