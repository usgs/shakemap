
import os
import subprocess
import tempfile

import numpy as np

from shakelib.rupture.factory import get_rupture
from shakelib.rupture.origin import Origin


homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))

# Dummy origin
dummy = {
    'mag': np.nan,
    'id': 'dummy',
    'locstring': 'dummy',
    'mech': 'ALL',
    'lon': np.nan,
    'lat': np.nan,
    'depth': np.nan,
    'netid': "",
    'network': "",
    'time': ""
}
origin = Origin(dummy)

program = os.path.join(shakedir, 'bin', 'sm_rupture')


def test_rupture():
    with tempfile.TemporaryDirectory() as tmpdir:
        # Read in a fault file
        rup1_file = os.path.join(
            shakedir, 'tests', 'data', 'eventdata', 'northridge', 'current',
            'northridge_fault.txt')
        rup1 = get_rupture(origin, rup1_file)

        # Known point is p0
        dx = 0
        dy = 0
        p0 = rup1.getQuadrilaterals()[0][0]
        px = p0.x
        py = p0.y
        pz = p0.z
        length = rup1.getLength()
        width = rup1.getWidth()
        strike = rup1.getStrike()
        dip = rup1.getDip()

        outfile = os.path.join(tmpdir, 'test.json')

        op = subprocess.Popen(
            [program, outfile],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=False
        )
        responses = (
            'test\n' +
            '1\n' +
            str(px) + '\n' +
            str(py) + '\n' +
            str(pz) + '\n' +
            str(dx) + '\n' +
            str(dy) + '\n' +
            str(length) + '\n' +
            str(width) + '\n' +
            str(strike) + '\n' +
            str(dip) + '\n'
        )
        op.communicate(responses.encode('ascii'))
        rup2 = get_rupture(origin, outfile)

        # testing, note that some difference will occur since the original
        # points are not necessarily coplanar or even rectangular, which
        # are conditions enfored on the derived rupture and so this cannot
        # be a very precise comparison.
        rtol = 1e-4
        np.testing.assert_allclose(rup2.lats, rup1.lats, rtol=rtol)
        np.testing.assert_allclose(rup2.lons, rup1.lons, rtol=rtol)
        rtol = 2e-3
        np.testing.assert_allclose(rup2.depths, rup1.depths, rtol=rtol)


if __name__ == '__main__':
    test_rupture()
