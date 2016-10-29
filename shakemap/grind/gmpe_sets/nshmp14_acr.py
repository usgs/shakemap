#!/usr/bin/env python

from openquake.hazardlib.gsim.abrahamson_2014 import AbrahamsonEtAl2014
from openquake.hazardlib.gsim.boore_2014 import BooreEtAl2014
from openquake.hazardlib.gsim.campbell_bozorgnia_2014 import CampbellBozorgnia2014
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014


def get_weights():
    """
    Method for getting the 2014 NSHMP ACR gmpe list and list of 
    their weights. Note that the large distance GMPE weights are
    None if they are not used (in which case there is no distance
    distinction in the GMPE weights). 

    TODO: Make a base class for gmpe_sets. 

    Returns:
        tuple: list of GMPEs, list of weights, list of large distance GMPE
            weights, distance cutoff (km), and list of GMPEs to use for site
            amplification. 
    """
    gmpes = [AbrahamsonEtAl2014(), BooreEtAl2014(),
             CampbellBozorgnia2014(), ChiouYoungs2014()]
    wts = [0.25, 0.25, 0.25, 0.25]

    # Notes:
    #    - No distance-depenent weights
    #    - All GMPEs have site term so no need for site GMPE list

    return gmpes, wts, None, None, None
