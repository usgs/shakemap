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

from shakemap.grind.gmpe_sets import nshmp14_scr_rlme

def test_nshmp14_scr_rlme():
    gmpes, wts, wts_large_dist, dist_cutoff, site_gmpes = \
        nshmp14_scr_rlme.get_weights()

    assert gmpes == \
        ['FrankelEtAl1996MwNSHMP2008()',
         'ToroEtAl1997MwNSHMP2008()',
         'SilvaEtAl2002MwNSHMP2008()',
         'Campbell2003MwNSHMP2008()',
         'TavakoliPezeshk2005MwNSHMP2008()',
         'AtkinsonBoore2006Modified2011()',
         'PezeshkEtAl2011()',
         'Atkinson2008prime()',
         'SomervilleEtAl2001NSHMP2008()']

    assert wts == [0.06, 0.11, 0.06, 0.11, 0.11, 0.22, 0.15, 0.08, 0.1]

    assert wts_large_dist == [0.16, 0.0, 0.0, 0.17, 0.17, 0.3, 0.2, 0.0, 0.0]

    assert dist_cutoff == 500

    assert site_gmpes == ['AtkinsonBoore2006Modified2011()']

