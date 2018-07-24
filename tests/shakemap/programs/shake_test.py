#!/usr/bin/env python

import os
import subprocess
import shutil

from shakemap.utils.config import get_config_paths

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))

########################################################################
# Test shake
########################################################################


def clear_files(event_path):
    files = ['boat_fault.txt', 'dyfi_dat.xml', 'event.xml', 'model.conf',
             'model_select.conf', 'moment.xml', 'stationlist.xml',
             'shake_data.hdf']
    for fname in files:
        try:
            os.remove(os.path.join(event_path, fname))
        except OSError:
            pass


def set_files(event_path, files):
    for src, dst in files.items():
        shutil.copy(os.path.join(event_path, 'data', src),
                    os.path.join(event_path, dst))


def test_shake():

    installdir, datadir = get_config_paths()

    program = os.path.join(shakedir, 'bin', 'shake')
    #
    # Run a bogus event
    #
    cp = subprocess.run([program, 'not_an_event', 'assemble', '-c',
                         '"Test Comment"'], shell=False)
    assert cp.returncode
    #
    # Run a real event
    #
    event_path = os.path.join(datadir, 'nc72282711', 'current')
    clear_files(event_path)
    set_files(event_path, {'event.xml': 'event.xml'})
    cp = subprocess.run([program, 'nc72282711', 'assemble', '-c',
                         '"Test comment"'], shell=False)
    clear_files(event_path)
    assert not cp.returncode

    # run an event with the exception module, which also generates
    # an exception.  See what happens...
    cp = subprocess.run([program, 'nc72282711', 'exception'], shell=False)
    clear_files(event_path)
    assert cp.returncode


########################################################################
# main program
########################################################################
if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_shake()
