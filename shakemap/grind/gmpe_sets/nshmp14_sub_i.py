#!/usr/bin/env python

from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SInter
from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SInter
from openquake.hazardlib.gsim.atkinson_macias_2009 import AtkinsonMacias2009
from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SInter

# Note: I'm not using the "NSHMP" versions of the OQ GMPEs because they
#       fix hypocentral depth at 20 km.

def get_weights():
    """
    Method for getting the 2014 NSHMP SUB/interface gmpe list and list of 
    their weights. Note that the large distance GMPE weights are
    None if they are not used (in which case there is no distance
    distinction in the GMPE weights). 

    TODO: Make a base class for gmpe_sets. 

    Returns:
        tuple: list of GMPEs, list of weights, list of large distance GMPE
            weights, distance cutoff (km), and list of GMPEs to use for site
            amplification. 
    """

    gmpes = [AtkinsonBoore2003SInter(), ZhaoEtAl2006SInter(),
             AtkinsonMacias2009(), AbrahamsonEtAl2015SInter()]
    wts = [0.1, 0.3, 0.3, 0.3]

    # Note:
    #    - No distance-depenent weights
    # TODO:
    #    - Figure out site GMPEs. 

    return gmpes, wts, None, None, None
