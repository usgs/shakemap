#!/usr/bin/env python

import os
import subprocess

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))

########################################################################
# Test shake
########################################################################


def test_shake():

    program = os.path.join(shakedir, 'bin', 'shake')
    #
    # Run a bogus event
    #
    cp = subprocess.run([program, 'not_an_event', 'assemble'], shell=False)
    assert cp.returncode
    #
    # Run a real event
    #
    cp = subprocess.run([program, 'nc72282711', 'assemble'], shell=False)
    assert not cp.returncode


########################################################################
# main program
########################################################################
if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_shake()
