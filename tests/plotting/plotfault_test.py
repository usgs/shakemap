import os.path
import sys
import tempfile
import hashlib

import matplotlib.pyplot as plt

from shakemap.plotting.plotfault import plot_fault_wire3d
from shakemap.grind.fault import Fault

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '../..'))
sys.path.insert(0, shakedir)

def test_plot_fault_wire3d():
    # make a temporary directory
    p = tempfile.mkdtemp()
    testfile = os.path.join(p, 'fault.png')

    ff = os.path.join(shakedir,
        "tests/data/eventdata/hayward_RC_HN_HS_HE_Shaw09Mod_GEOL.txt")
    flt = Fault.readFaultFile(ff)
    ax = plot_fault_wire3d(flt)
    print(testfile)
    plt.savefig(testfile)
    data = open(testfile, 'rb').read()
    m = hashlib.md5()
    m.update(data)
    assert m.digest() == b'\xdd\r\x992\x9f\xdeN6\xa5Be\xc7\xa2*\xee\xe3'


if __name__ == '__main__':
    test_plot_fault_wire3d()
