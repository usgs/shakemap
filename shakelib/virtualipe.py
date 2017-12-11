import copy

import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.imt import PGA, PGV, SA, MMI
from openquake.hazardlib import const

from shakelib.utils.exception import ShakeLibException


class VirtualIPE(GMPE):
    """
    Implements a virtual IPE that is the combination of a MultiGMPE
    and a GMICE. Will first attempt to use PGV for the intensities,
    and then will try PGA, and then SA(1.0), and then will bail out.

    Uncertainty is computed by combining the uncertainty of the GMPE
    with the uncertainty of the GMICE. Standard error propagation
    techniques are used (see the ShakeMap manual for a detailed
    explanation). For the intra- and inter-event components of total
    uncertainty, we assign all of GMICE uncertaninty to the intra-event
    term, and none to the inter-event term. This choice is conservative,
    and seems appropriate until GMICE are produced with separate inter-
    and intra-event terms.

    Note that the combined inter- and intra-event uncertainties will
    only approximately equal the total uncertainty because the GMPEs
    will only produce combined uncertainties that are approximately
    equal to their total uncertainty.

    """

    #: The OpenQuake IMT this module can produce ('MMI' only).
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([MMI])
    #: The OpenQuake standard deviation types that may be produced (will
    #: depend on the GMPE provided).
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = None
    #: Distance measures required (will depend on the GMPE provided).
    REQUIRES_DISTANCES = None
    #: OpenQuake IMC used (will depend on the GMPE, but "Larger" is 
    #: typical).
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    #: Determined by the GMPE selected.
    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    #: Determined by the GMPE selected.
    REQUIRES_RUPTURE_PARAMETERS = None
    #: Determined by the GMPE selected.
    REQUIRES_SITES_PARAMETERS = None

    @classmethod
    def fromFuncs(cls, gmpe, gmice):
        """
        Creates a new VirtualIPE object with the specified MultiGMPE and
        GMICE. There is no default constructor, you must use this method.

        Args:
            gmpe: An instance of the MultiGMPE object.
            gmice: An instance of a GMICE object.

        Returns:
            :class:`VirtualIPE`: A new instance of a VirtualIPE object.

        """
        self = cls()
        self.gmpe = gmpe
        self.gmice = gmice

        if (gmpe.ALL_GMPES_HAVE_PGV is True and
                PGV in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
            self.imt = PGV()
        elif (PGA in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES and
              PGA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
            self.imt = PGA()
        elif (SA in gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES and
              SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES):
            self.imt = SA(1.0)
        else:
            raise ShakeLibException(
                'The supplied GMPE and GMICE do not have a common IMT'
            )

        self.DEFINED_FOR_STANDARD_DEVIATION_TYPES = \
            gmpe.DEFINED_FOR_STANDARD_DEVIATION_TYPES.copy()

        self.REQUIRES_DISTANCES = gmpe.REQUIRES_DISTANCES.copy()
        self.REQUIRES_RUPTURE_PARAMETERS = \
                gmpe.REQUIRES_RUPTURE_PARAMETERS.copy()
        self.REQUIRES_SITES_PARAMETERS = \
                gmpe.REQUIRES_SITES_PARAMETERS.copy()
        self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = \
                copy.copy(gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT)
        self.DEFINED_FOR_TECTONIC_REGION_TYPE = \
                copy.copy(gmpe.DEFINED_FOR_TECTONIC_REGION_TYPE)

        return self

    def get_mean_and_stddevs(self, sx, rx, dx, imt, stddev_types, fd=None):
        """
        See superclass
        `method <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.GroundShakingIntensityModel.get_mean_and_stddevs>`__
        for parameter definitions. The only acceptable IMT is MMI.

        Additional subclass argument is "fd", which is the directivity
        amplification factor in natural log units. This is optional,
        and must be a numpy array with the same dimentions as the
        sites and is added to the ground motions before conversion to
        MMI.

        Returns:
            ndarray, list of ndarray:

            mmi (ndarray): Ground motions predicted by the MultiGMPE using
            the supplied parameters are converted to MMI using the GMICE.

            mmi_sd (list of ndarrays): The uncertainty of the combined
            prediction/conversion process. The prediction uncertainty will
            typically be either OpenQuake's TOTAL or INTRA_EVENT.  But can
            be any set that the MultiGMPE supports. See the ShakeMap manual
            for a detailed discussion of the way the uncertainty is computed.

        """ # noqa

        if imt != MMI():
            raise ValueError("imt must be MMI")
        #
        # Get the mean ground motions and stddev for the preferred IMT
        #
        mgm, sdev = self.gmpe.get_mean_and_stddevs(sx, rx, dx, self.imt,
                                                   stddev_types)

        if fd is not None:
            mgm = mgm + fd

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
        # Total and intra-event uncertanty are inflated by the uncertainty
        # of the conversion; inter-event uncertainty is not.
        #
        nsd = len(sdev)
        mmi_sd = [None] * nsd
        gm2mi_var = (self.gmice.getGM2MIsd()[self.imt])**2
        dmda *= dmda
        for i in range(nsd):
            gm_var_in_mmi = dmda * sdev[i]**2
            if stddev_types[i] == const.StdDev.INTER_EVENT:
                mmi_sd[i] = np.sqrt(gm_var_in_mmi)
            else:
                mmi_sd[i] = np.sqrt(gm2mi_var + gm_var_in_mmi)

        return mmi, mmi_sd
