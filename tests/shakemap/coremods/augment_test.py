#!/usr/bin/env python
import os
import shutil
import subprocess

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.augment import AugmentModule
from shakemap.coremods.assemble import AssembleModule

########################################################################
# Test augment coremod
########################################################################


def test_augment():
    installpath, datapath = get_config_paths()

    # Process a non-existent event (should fail)
    augment = AugmentModule('not_an_event', comment='Test comment.')
    with pytest.raises(NotADirectoryError):
        augment.execute()

    # This would succeed, but we remove shake_data.hdf (should fail)
    hdf_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.rename(hdf_file, hdf_file + '_safe')
    try:
        augment = AugmentModule('wenchuan', comment='Test comment.')
        with pytest.raises(FileNotFoundError):
            augment.execute()
    finally:
        if os.path.isfile(hdf_file + '_safe'):
            os.rename(hdf_file + '_safe', hdf_file)

    # Normal event (should succeed)
    data_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        assemble = AssembleModule('wenchuan', comment='Test comment.')
        assemble.execute()
        augment = AugmentModule('wenchuan', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    # Do an event with model.conf (not model_zc.conf) and no zoneinfo
    # (should succeed)
    data_file = os.path.join(datapath, 'nc72282711',
                             'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        augment = AugmentModule('nc72282711', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    #
    # Make sure the location file substitutions work (should succeed)
    #
    data_file = os.path.join(
        datapath, 'northridge_points', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        assemble = AssembleModule('northridge_points', comment='Test comment.')
        assemble.execute()
        augment = AugmentModule('northridge_points', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    #
    # Try some bad config files
    #
    data_file = os.path.join(
        datapath, 'nc72282711_nodata_nofault', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        assemble = AssembleModule('nc72282711_nodata_nofault',
                                  comment='Test comment.')
        assemble.execute()
        # Should fail validation
        model_file = os.path.join(datapath, 'nc72282711_nodata_nofault',
                                  'current', 'model_zc.conf')
        os.rename(model_file, model_file + '_safe')
        shutil.copyfile(model_file + '.bad0', model_file)
        try:
            augment = AugmentModule('nc72282711_nodata_nofault',
                                    comment='Test comment.')
            with pytest.raises(RuntimeError):
                augment.execute()
        finally:
            os.rename(model_file + '_safe', model_file)

        # Should fail vs30 filename check
        model_file = os.path.join(datapath, 'nc72282711_nodata_nofault',
                                  'current', 'model_zc.conf')
        os.rename(model_file, model_file + '_safe')
        shutil.copyfile(model_file + '.bad1', model_file)
        try:
            augment = AugmentModule('nc72282711_nodata_nofault',
                                    comment='Test comment.')
            with pytest.raises(FileNotFoundError):
                augment.execute()
        finally:
            os.rename(model_file + '_safe', model_file)

        # Should fail prediction locations filename check
        model_file = os.path.join(datapath, 'nc72282711_nodata_nofault',
                                  'current', 'model_zc.conf')
        os.rename(model_file, model_file + '_safe')
        shutil.copyfile(model_file + '.bad2', model_file)
        try:
            augment = AugmentModule('nc72282711_nodata_nofault',
                                    comment='Test comment.')
            with pytest.raises(FileNotFoundError):
                augment.execute()
        finally:
            os.rename(model_file + '_safe', model_file)
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    #
    # Switch originators (should succeed)
    #
    model_file = os.path.join(datapath, 'nc72282711',
                              'current', 'model.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.cz', model_file)
    data_file = os.path.join(datapath, 'nc72282711',
                             'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        augment = AugmentModule('nc72282711', comment='Test comment.')
        augment.execute()
    finally:
        os.rename(model_file + '_safe', model_file)
        if os.path.isfile(data_file):
            os.remove(data_file)


def test_augment_command_line():
    installpath, datapath = get_config_paths()

    cp = subprocess.run(['shake', '--force', 'integration_test_0001',
                         'assemble', '-c', 'Assemble comment'], shell=False)
    assert not cp.returncode
    cp = subprocess.run(['shake', '--force', 'integration_test_0001',
                         'augment', '-c', 'Augment comment'], shell=False)
    assert not cp.returncode
    data_file = os.path.join(datapath, 'integration_test_0001', 'current',
                             'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_augment()
    test_augment_command_line()
