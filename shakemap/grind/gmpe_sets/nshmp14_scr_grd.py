#!/usr/bin/env python

from openquake.hazardlib.gsim.frankel_1996 import FrankelEtAl1996MwNSHMP2008
# Note: the Toro implementation includes the 2002 update.
from openquake.hazardlib.gsim.toro_1997 import ToroEtAl1997MwNSHMP2008
from openquake.hazardlib.gsim.silva_2002 import SilvaEtAl2002MwNSHMP2008
from openquake.hazardlib.gsim.campbell_2003 import Campbell2003MwNSHMP2008
from openquake.hazardlib.gsim.tavakoli_pezeshk_2005 import TavakoliPezeshk2005MwNSHMP2008
from openquake.hazardlib.gsim.atkinson_boore_2006 import AtkinsonBoore2006Modified2011
from openquake.hazardlib.gsim.pezeshk_2011 import PezeshkEtAl2011
from openquake.hazardlib.gsim.boore_atkinson_2011 import Atkinson2008prime
from openquake.hazardlib.gsim.somerville_2001 import SomervilleEtAl2001NSHMP2008


def get_weights():
    """
    Method for getting the 2014 NSHMP SCR GRD gmpe list and list of
    their weights. Note that the large distance GMPE weights are
    None if they are not used (in which case there is no distance
    distinction in the GMPE weights). 

    TODO: Make a base class for gmpe_sets. 

    Returns:
        tuple: list of GMPEs, list of weights, list of large distance GMPE
            weights, distance cutoff (km), and list of GMPEs to use for site
            amplification. 
    """
    gmpes = [FrankelEtAl1996MwNSHMP2008(),
             ToroEtAl1997MwNSHMP2008(),
             SilvaEtAl2002MwNSHMP2008(),
             Campbell2003MwNSHMP2008(),
             TavakoliPezeshk2005MwNSHMP2008(),
             AtkinsonBoore2006Modified2011(),
             PezeshkEtAl2011(),
             Atkinson2008prime(),
             SomervilleEtAl2001NSHMP2008()]
    wts = [0.06, 0.13, 0.06, 0.13, 0.13, 0.25, 0.16, 0.08, 0.0]
    wts_large_dist = [0.16, 0.0, 0.0, 0.17, 0.17, 0.3, 0.2, 0.0, 0.0]
    dist_cutoff = 500
    site_gmpes = [AtkinsonBoore2006Modified2011()]

    return gmpes, wts, wts_large_dist, dist_cutoff, site_gmpes
