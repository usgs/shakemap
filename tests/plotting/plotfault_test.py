import os.path
import sys

import matplotlib.pyplot as plt
import pytest

from shakemap.plotting.plotrupture import plot_rupture_wire3d
from shakemap.plotting.plotrupture import map_rupture
from shakemap.grind.rupture import read_rupture_file
from shakemap.grind.origin import Origin

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

@pytest.mark.mpl_image_compare(tolerance=4)
def test_plot_rupture_wire3d_QuadRupture():
    origin = Origin({'id':'','lat':0,'lon':0,'depth':0,'mag':0})

    # A relatively complicated QuadRupture
    ff = os.path.join(shakedir,
        "tests/data/Barkaetal02_fault.txt")
    rup = read_rupture_file(origin, ff)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_rupture_wire3d(rup, ax)
    return fig

@pytest.mark.mpl_image_compare(tolerance=4)
def test_plot_rupture_wire3d_EdgeRupture():
    origin = Origin({'id':'','lat':0,'lon':0,'depth':0,'mag':0})

    # Cascadia
    ff = os.path.join(shakedir, "tests/data/cascadia.json")
    rup = read_rupture_file(origin, ff)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plot_rupture_wire3d(rup, ax)
    return fig

@pytest.mark.mpl_image_compare(tolerance=4)
def test_map_QuadRupture():
    origin = Origin({'id':'','lat':0,'lon':0,'depth':0,'mag':0})

    # Ismit
    ff = os.path.join(shakedir, 
        "tests/data/izmit.json")
    rup = read_rupture_file(origin, ff)
    fig = plt.figure()
    map_rupture(rup)
    return fig

@pytest.mark.mpl_image_compare(tolerance=4)
def test_map_EdgeRupture():
    origin = Origin({'id':'','lat':0,'lon':0,'depth':0,'mag':0})

    # Cascadia
    ff = os.path.join(shakedir, "tests/data/cascadia.json")
    rup = read_rupture_file(origin, ff)
    fig = plt.figure()
    map_rupture(rup)
    return fig


if __name__ == '__main__':
    test_plot_rupture_wire3d_QuadRupture()
    test_plot_rupture_wire3d_EdgeRupture()
