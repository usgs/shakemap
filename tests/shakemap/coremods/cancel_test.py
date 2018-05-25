#!/usr/bin/env python
import os
import shutil
# import subprocess

# import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.model import ModelModule
from shakemap.coremods.cancel import CancelModule

########################################################################
# Test cancel coremod
########################################################################


def test_cancel():
    installpath, datapath = get_config_paths()

    shutil.copytree(os.path.join(datapath, 'integration_test_0001'),
                    os.path.join(datapath, 'dummy_event'))

    # Process a non-existent event (should fail)
    try:
        mod = AssembleModule('dummy_event', comment='')
        mod.execute()
        mod = ModelModule('dummy_event')
        mod.execute()
        mod = CancelModule('dummy_event')
        mod.execute()
    finally:
        shutil.rmtree(os.path.join(datapath, 'dummy_event'))
        amp_file = os.path.join(installpath, 'data', 'amps.db')
        if os.path.isfile(amp_file):
            os.remove(amp_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_cancel()
