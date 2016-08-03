"""
Routines to convert GMPE amps to/from ShakeMap values
GMPE amps are in ln(g), in an IMC of their choosing
ShakeMap amps are %g and max horizontal
"""

from shakemap.grind.conversions.imc.beyer_bommer_2006 import BeyerBommer2006 as bb06

from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib import const


def ampShakeMapToGMPE(sm_amps, gmpe, imt):
    """
    Convert Shakemap ground motion amplitudes to match the 
    intensity measure type 
    (`IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__)
    and intensity measure component
    (`IMC <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__)
    of a GMPE. 

    :param sa_amps:
        Numpy array of ground motion amplitudes. 
    :param gmpe:
        An Openquake `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__. 
    :param imt:
        An Openquake `IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__. 
    :returns:
        Numpy array of converted ground motion amplitudes. 
    """

    if not isinstance(imt, PGV):
        gmpe_amps = sm_amps / 100
    else:
        gmpe_amps = sm_amps.copy()

    gmpe_amps = bb06.ampIMCtoIMC(gmpe_amps, const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                 gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, imt)
    gmpe_amps = gmpe.to_distribution_values(gmpe_amps)

    return gmpe_amps


def ampGmpeToShakeMap(gmpe_amps, gmpe, imt):
    """
    Convert ground motion amplitudes from a GMPE to match the 
    intensity measure type 
    (`IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__)
    and intensity measure component
    (`IMC <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__)
    used by Shakemap. 

    :param gmpe_amps:
        Numpy array of ground motion amplitudes from a GMPE. 
    :param gmpe:
        An Openquake `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__. 
    :param imt:
        An Openquake `IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__. 
    :returns:
        Numpy array of converted ground motion amplitudes. 
    """

    sm_amps = gmpe.to_imt_unit_values(gmpe_amps)
    sm_amps = bb06.ampIMCtoIMC(sm_amps, gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT,
                               const.IMC.GREATER_OF_TWO_HORIZONTAL, imt)

    if not isinstance(imt, PGV):
        sm_amps *= 100

    return sm_amps


def sigmaShakeMapToGMPE(sm_sigmas, gmpe, imt):
    """
    Convert Shakemap ground motion standard deviations to match the 
    intensity measure type 
    (`IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__)
    and intensity measure component
    (`IMC <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__)
    of a GMPE. 

    :param sa_sigmas:
        Numpy array of ground motion standard deviations. 
    :param gmpe:
        An Openquake `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__. 
    :param imt:
        An Openquake `IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__. 
    :returns:
        Numpy array of converted ground motion standard deviations. 
    """

    gmpe_sigmas = gmpe.to_distribution_values(sm_sigmas)
    gmpe_sigmas = bb06.sigmaIMCtoIMC(gmpe_sigmas, const.IMC.GREATER_OF_TWO_HORIZONTAL,
                                     gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT, imt)

    return gmpe_sigmas


def sigmaGmpeToShakeMap(gmpe_sigmas, gmpe, imt):
    """
    Convert ground motion standard deviations from a GMPE to match the 
    intensity measure type 
    (`IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__)
    and intensity measure component
    (`IMC <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__)
    used by Shakemap. 

    :param gmpe_sigmas:
        Numpy array of ground motion standard deviations from a GMPE. 
    :param gmpe:
        An Openquake `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__. 
    :param imt:
        An Openquake `IMT <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__. 
    :returns:
        Numpy array of converted ground motion standard deviations. 
    """

    sm_sigmas = bb06.sigmaIMCtoIMC(gmpe_sigmas, gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT,
                                   const.IMC.GREATER_OF_TWO_HORIZONTAL, imt)
    sm_sigmas = gmpe.to_imt_unit_values(sm_sigmas)

    return sm_sigmas
