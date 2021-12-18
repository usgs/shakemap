"""
A GMPE that returns a constant everywhere. Useful for testing, but
nothing else.
"""
import numpy as np

from openquake.hazardlib.gsim.base import GMPE, CoeffsTable
from openquake.hazardlib import const
from openquake.hazardlib.imt import PGA, PGV, SA


class NullGMPE(GMPE):
    """
    This is a GMPE for testing. It returns the mean and stddevs
    specified in the constructor.
    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.ACTIVE_SHALLOW_CRUST
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([PGA, PGV, SA])
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.GREATER_OF_TWO_HORIZONTAL
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set(
        [const.StdDev.TOTAL, const.StdDev.INTER_EVENT, const.StdDev.INTRA_EVENT]
    )
    REQUIRES_SITES_PARAMETERS = set(("vs30",))
    REQUIRES_RUPTURE_PARAMETERS = set(())
    REQUIRES_DISTANCES = set(("rjb",))

    def __init__(self, mean=0, phi=0.8, tau=0.6):
        """
        The default constructor takes three named arguments:

        Args:
            mean (float): the mean value returned by the GMPE (default=0).
                This value is returned for all locations, regardles of
                the IMT or the contents of sites, rupture and distance
                contexts.
            phi (float): the within-event standard deviation returned by the
                 GMPE (default=0.8)
            tau (float): the between-event standard deviation returned by the
                 GMPE (default=0.6)

        The total standard deviation returned will be ``sqrt(phi^2 + tau^2)``.
        """
        self.mean = mean
        self.phi = phi
        self.tau = tau
        self.sigma = np.sqrt(phi ** 2 + tau ** 2)

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        Implements the OpenQuake GroundShakingIntensityModel
        get_mean_and_stddevs interface. See superclass
        `method <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.GroundShakingIntensityModel.get_mean_and_stddevs>`__.

        Returns a constant values for all locations specified in
        the dists.rbj array, regardless of the contents of that array or
        any of the other contexts. The imt is also ignored.
        """  # noqa

        mean = np.full_like(dists.rjb, self.mean)
        stddevs = []
        for stddev_type in stddev_types:
            if stddev_type == const.StdDev.TOTAL:
                stddevs.append(
                    np.full_like(dists.rjb, np.sqrt(self.phi ** 2 + self.tau ** 2))
                )
            elif stddev_type == const.StdDev.INTRA_EVENT:
                stddevs.append(np.full_like(dists.rjb, self.phi))
            elif stddev_type == const.StdDev.INTER_EVENT:
                stddevs.append(np.full_like(dists.rjb, self.tau))

        return mean, stddevs

    #
    # Dummy COEFFS table so MultiGMPE won't complain
    #
    COEFFS = CoeffsTable(
        sa_damping=5,
        table="""\
IMT     c1
pga    0.
pgv    0.
0.01   0.
10     0.
""",
    )
