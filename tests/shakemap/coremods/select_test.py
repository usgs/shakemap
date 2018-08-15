#!/usr/bin/env python
import os
import os.path

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.select import SelectModule
from common import clear_files, set_files


########################################################################
# Test select
########################################################################
def test_select():

    installpath, datapath = get_config_paths()

    # Process a non-existent event (should fail)
    smod = SelectModule('not_an_event')
    with pytest.raises(NotADirectoryError):
        smod.execute()

    # Normal event (should succeed)
    event_path = os.path.join(datapath, 'nc72282711', 'current')
    set_files(event_path, {'event.xml': 'event.xml'})
    conf_file = os.path.join(datapath, 'nc72282711', 'current',
                             'model_select.conf')
    smod = SelectModule('nc72282711')
    smod.execute()
    failed = False
    if not os.path.isfile(conf_file):
        failed = True
    clear_files(event_path)
    if failed:
        assert False

    # Subduction event (not over slab)
    conf_file = os.path.join(datapath, 'usp0004bxs', 'current',
                             'model_select.conf')
    if os.path.isfile(conf_file):
        os.remove(conf_file)
    try:
        smod = SelectModule('usp0004bxs')
        smod.execute()
    finally:
        if not os.path.isfile(conf_file):
            print('select failed!')
            assert False
        else:
            os.remove(conf_file)

    # Northridge, with moment tensor file
    conf_file = os.path.join(datapath, 'northridge2', 'current',
                             'model_select.conf')
    if os.path.isfile(conf_file):
        os.remove(conf_file)
    try:
        smod = SelectModule('northridge2')
        smod.execute()
    finally:
        if not os.path.isfile(conf_file):
            print('select failed!')
            assert False
        else:
            os.remove(conf_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_select()
