#!/usr/bin/env python

import os
import subprocess

from shakemap.coremods.sleep import SleepModule


########################################################################
# Test assemble and augment
########################################################################

def test_sleep():

    amod = SleepModule('northridge', seconds=0)
    amod.execute()


def test_sleep_command_line():

    cp = subprocess.run(['shake', 'integration_test_0001',
                         'sleep', '-s', '0'], shell=False)
    assert not cp.returncode


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_sleep()
    test_sleep_command_line()
