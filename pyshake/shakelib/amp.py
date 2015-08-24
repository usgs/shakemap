"""
Routines to convert GMPE amps to/from ShakeMap values
GMPE amps are in ln(g), in an IMC of their choosing
ShakeMap amps are %g and max horizontal
"""

import pyshake.shakelib.BeyerBommer2006 as BB

from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib import const

def shakemapToGMPE(sm_amps, gmpe, imt):

    if not isinstance(imt, PGV):
        gmpe_amps = sm_amps / 100

    gmpe_amps = BB.ampIMCtoIMC(gmpe_amps, const.IMC.GREATER_OF_TWO_HORIZONTAL, \
                               gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, imt)
    gmpe_amps = gmpe.to_distribution_values(gmpe_amps)

    return gmpe_amps

def gmpeToShakeMap(gmpe_amps, gmpe, imt):

    sm_amps = gmpe.to_imt_unit_values(gmpe_amps)
    sm_amps = BB.ampIMCtoIMC(sm_amps, gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, \
                             const.IMC.GREATER_OF_TWO_HORIZONTAL, imt)

    if not isinstance(imt, PGV):
        sm_amps *= 100

    return sm_amps
