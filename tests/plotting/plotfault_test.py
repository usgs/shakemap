import os.path
import sys

import matplotlib.pyplot as plt
import pytest

from shakemap.plotting.plotrupture import plot_rupture_wire3d
from shakemap.grind.rupture import read_rupture_file
from shakemap.grind.origin import Origin

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

@pytest.mark.mpl_image_compare(tolerance=4)
def test_plot_rupture_wire3d():
    origin = Origin({'id':'','lat':0,'lon':0,'depth':0,'mag':0})
    ff = os.path.join(shakedir,
        "tests/data/eventdata/hayward_RC_HN_HS_HE_Shaw09Mod_GEOL.txt")
    flt = read_rupture_file(origin, ff)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_rupture_wire3d(flt, ax)
    return fig

if __name__ == '__main__':
    test_plot_rupture_wire3d()
