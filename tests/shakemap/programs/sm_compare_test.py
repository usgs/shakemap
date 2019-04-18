#!/usr/bin/env python

import os
import subprocess

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))

########################################################################
# Test shake
########################################################################


def test_sm_compare():

    program = os.path.join(shakedir, 'bin', 'sm_compare')

    #
    # Run a real event
    #
    data_dir = os.path.join(shakedir, 'tests', 'data', 'eventdata',
                            'comparison_test')
    grid1 = os.path.join(data_dir, 'shakemap3.xml')
    grid2 = os.path.join(data_dir, 'shakemap4.xml')
    cp = subprocess.run([program, grid1, grid2, '--nocoasts'], shell=False)
    if os.path.isfile('compare.png'):
        os.remove('compare.png')
    assert not cp.returncode


########################################################################
# main program
########################################################################
if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_sm_compare()
