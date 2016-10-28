#!/usr/bin/env python

from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlab
from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlabCascadia
from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlab
from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab


def get_weights():
    """
    Method for getting the 2014 NSHMP SUB/slab gmpe list and list of 
    their weights. Note that the large distance GMPE weights are
    None if they are not used (in which case there is no distance
    distinction in the GMPE weights). 

    TODO: Make a base class for gmpe_sets. 

    Returns:
        tuple: list of GMPEs, list of weights, list of large distance GMPE
            weights, distance cutoff (km), and list of GMPEs to use for site
            amplification. 
    """

    gmpes = [AtkinsonBoore2003SSlab(), AtkinsonBoore2003SSlabCascadia(),
             ZhaoEtAl2006SSlab(), AbrahamsonEtAl2015SSlab()]
    wts = [0.1667, 0.1667, 0.33, 0.33]

    # Note:
    #    - No distance-depenent weights
    # TODO:
    #    - Figure out site GMPEs. 

    return gmpes, wts, None, None, None
