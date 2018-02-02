#!/usr/bin/env python
import os
import os.path
import shutil

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.select import SelectModule


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
    conf_file = os.path.join(datapath, 'nc72282711', 'current',
                             'model_zc.conf')
    if os.path.isfile(conf_file):
        os.remove(conf_file)
    try:
        smod = SelectModule('nc72282711')
        smod.execute()
    finally:
        if not os.path.isfile(conf_file):
            print('select failed!')
            assert False
        else:
            os.remove(conf_file)

    # Subduction event (should succeed)
    conf_file = os.path.join(datapath, 'official20041226005853450_30', 'current',
                             'model_zc.conf')
    if os.path.isfile(conf_file):
        os.remove(conf_file)
    try:
        smod = SelectModule('official20041226005853450_30')
        smod.execute()
    finally:
        if not os.path.isfile(conf_file):
            print('select failed!')
            assert False
        else:
            os.remove(conf_file)

    # Subduction event (near edge of subduction region)
    conf_file = os.path.join(datapath, 'usp000bnyr', 'current',
                             'model_zc.conf')
    if os.path.isfile(conf_file):
        os.remove(conf_file)
    try:
        smod = SelectModule('usp000bnyr')
        smod.execute()
    finally:
        if not os.path.isfile(conf_file):
            print('select failed!')
            assert False
        else:
            os.remove(conf_file)

    # Subduction event (not over slab)
    conf_file = os.path.join(datapath, 'usp0004bxs', 'current',
                             'model_zc.conf')
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


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_select()
