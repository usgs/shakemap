#!/usr/bin/env python

from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlab
from openquake.hazardlib.gsim.atkinson_boore_2003 import AtkinsonBoore2003SSlabCascadia
from openquake.hazardlib.gsim.zhao_2006 import ZhaoEtAl2006SSlab
from openquake.hazardlib.gsim.abrahamson_2015 import AbrahamsonEtAl2015SSlab


from shakemap.grind.gmpe_sets import nshmp14_sub_s


def test_nshmp14_sub_i():
    gmpes, wts, wts_large_dist, dist_cutoff, site_gmpes = \
        nshmp14_sub_s.get_weights()

    assert gmpes == \
        ['AtkinsonBoore2003SSlab()',
         'AtkinsonBoore2003SSlabCascadia()',
         'ZhaoEtAl2006SSlab()',
         'AbrahamsonEtAl2015SSlab()']

    assert wts == [0.1667, 0.1667, 0.33, 0.33]

    assert wts_large_dist is None

    assert dist_cutoff is None

#    assert site_gmpes == 

