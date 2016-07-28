"""
Routines to convert GMPE amps to/from ShakeMap values
GMPE amps are in ln(g), in an IMC of their choosing
ShakeMap amps are %g and max horizontal
"""

from . import BeyerBommer2006 as BB

from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib import const


def ampShakeMapToGMPE(sm_amps, gmpe, imt):

    if not isinstance(imt, PGV):
        gmpe_amps = sm_amps / 100
    else:
        gmpe_amps = sm_amps.copy()

    gmpe_amps = BB.ampIMCtoIMC(gmpe_amps, const.IMC.GREATER_OF_TWO_HORIZONTAL,
                               gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, imt)
    gmpe_amps = gmpe.to_distribution_values(gmpe_amps)

    return gmpe_amps


def ampGmpeToShakeMap(gmpe_amps, gmpe, imt):

    sm_amps = gmpe.to_imt_unit_values(gmpe_amps)
    sm_amps = BB.ampIMCtoIMC(sm_amps, gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT,
                             const.IMC.GREATER_OF_TWO_HORIZONTAL, imt)

    if not isinstance(imt, PGV):
        sm_amps *= 100

    return sm_amps


def sigmaShakeMapToGMPE(sm_sigmas, gmpe, imt):

    gmpe_sigmas = gmpe.to_distribution_values(sm_sigmas)
    gmpe_sigmas = BB.sigmaIMCtoIMC(gmpe_sigmas, const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                   gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, imt)

    return gmpe_sigmas


def sigmaGmpeToShakeMap(gmpe_sigmas, gmpe, imt):

    sm_sigmas = BB.sigmaIMCtoIMC(gmpe_sigmas, gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT,
                                 const.IMC.GREATER_OF_TWO_HORIZONTAL, imt)
    sm_sigmas = gmpe.to_imt_unit_values(sm_sigmas)

    return sm_sigmas
