#!/usr/bin/env python

import os
import os.path

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.associate import AssociateModule


########################################################################
# Test sm_assemble
########################################################################
def test_associate():

    installpath, datapath = get_config_paths()

    amp_file = os.path.join(installpath, 'data', 'amps.db')
    if os.path.isfile(amp_file):
        os.remove(amp_file)

    # Process a non-existent event (should fail)
    amod = AssociateModule('not_an_event')
    with pytest.raises(NotADirectoryError):
        amod.execute()

    # Would succeed but we remove event.xml (should fail)
    event_file = os.path.join(datapath, 'integration_test_0001', 'current',
                              'event.xml')
    os.rename(event_file, event_file + '_safe')
    try:
        amod = AssociateModule('integration_test_0001')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        os.rename(event_file + '_safe', event_file)

    # Normal event (should succeed)
    try:
        amod = AssociateModule('integration_test_0001')
        amod.execute()
    finally:
        amp_file = os.path.join(installpath, 'data', 'amps.db')
        if os.path.isfile(amp_file):
            os.remove(amp_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_associate()
