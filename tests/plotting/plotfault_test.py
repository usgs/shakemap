import os.path
import sys

import matplotlib.pyplot as plt
import pytest

from shakemap.plotting.plotfault import plot_fault_wire3d
from shakemap.grind.fault import Fault

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

@pytest.mark.mpl_image_compare
def test_plot_fault_wire3d():
    ff = os.path.join(shakedir,
        "tests/data/eventdata/hayward_RC_HN_HS_HE_Shaw09Mod_GEOL.txt")
    flt = Fault.readFaultFile(ff)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_fault_wire3d(flt, ax)
    return fig

if __name__ == '__main__':
    test_plot_fault_wire3d()
