#!/usr/bin/env python

import os
import subprocess

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.augment import AugmentModule
from common import clear_files, set_files


########################################################################
# Test assemble and augment
########################################################################

def test_assemble_augment():

    installpath, datapath = get_config_paths()

    # Process a non-existent event (should fail)
    amod = AssembleModule('not_an_event', comment='Test comment.')
    with pytest.raises(NotADirectoryError):
        amod.execute()

    # Process a non-existent event (should fail)
    augment = AugmentModule('not_an_event', comment='Test comment.')
    with pytest.raises(NotADirectoryError):
        augment.execute()

    # Would succeed but we remove event.xml (should fail)
    event_file = os.path.join(datapath, 'wenchuan', 'current', 'event.xml')
    os.rename(event_file, event_file + '_safe')
    try:
        amod = AssembleModule('wenchuan', comment='Test comment.')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        os.rename(event_file + '_safe', event_file)

    # Normal event (should succeed)
    data_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        amod = AssembleModule('wenchuan', comment='Test comment.')
        amod.execute()
        #
        # Run a second time to exercise a different branch of the code
        #
        amod.execute()
    finally:
        pass

    # This would succeed, but we remove shake_data.hdf (should fail)
    hdf_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.rename(hdf_file, hdf_file + '_safe')
    try:
        augment = AugmentModule('wenchuan', comment='Test comment.')
        with pytest.raises(FileNotFoundError):
            augment.execute()
    finally:
        # Put shake_data.hdf back for the next test
        if os.path.isfile(hdf_file + '_safe'):
            os.rename(hdf_file + '_safe', hdf_file)

    # Normal event (should succeed)
    try:
        augment = AugmentModule('wenchuan', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    # Do an event with model.conf (not model_select.conf) and no zoneinfo
    # (should succeed)
    event_path = os.path.join(datapath, 'nc72282711', 'current')
    set_files(event_path, {'event.xml': 'event.xml',
                           'stationlist.xml': 'stationlist.xml',
                           'boat_fault.txt': 'boat_fault.txt',
                           'model.conf': 'model.conf'})
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
        clear_files(event_path)

    # Try not having an event-specific config (should succeed)
    set_files(event_path, {'event.xml': 'event.xml',
                           'stationlist.xml': 'stationlist.xml',
                           'boat_fault.txt': 'boat_fault.txt'})
    try:
        amod = AssembleModule('nc72282711', comment='Test comment.')
        amod.execute()
        augment = AugmentModule('nc72282711', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)
        clear_files(event_path)

    # Do an event with DYFI data (should succeed)
    set_files(event_path, {'event.xml': 'event.xml',
                           'stationlist.xml': 'stationlist.xml',
                           'boat_fault.txt': 'boat_fault.txt',
                           'dyfi_dat.xml': 'dyfi_dat.xml'})
    try:
        amod = AssembleModule('nc72282711', comment='Test comment.')
        amod.execute()
        augment = AugmentModule('nc72282711', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)
        clear_files(event_path)

    #
    # Try some bad config files
    #
    # Should fail validation
    set_files(event_path, {'event.xml': 'event.xml',
                           'model_select.conf.bad0': 'model_select.conf'})
    try:
        amod = AssembleModule('nc72282711', comment='Test comment.')
        with pytest.raises(RuntimeError):
            amod.execute()
    finally:
        clear_files(event_path)

    # Should fail vs30 filename check
    set_files(event_path, {'event.xml': 'event.xml',
                           'model_select.conf.bad1': 'model_select.conf'})
    try:
        amod = AssembleModule('nc72282711', comment='Test comment.')
        with pytest.raises(RuntimeError):
            amod.execute()
    finally:
        clear_files(event_path)

    # Should fail prediction locations filename check
    set_files(event_path, {'event.xml': 'event.xml',
                           'model_select.conf.bad2': 'model_select.conf'})
    try:
        amod = AssembleModule('nc72282711', comment='Test comment.')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        clear_files(event_path)
    #
    # Make sure the location file substitutions work (should succeed)
    #
    data_file = os.path.join(datapath, 'northridge_points', 'current',
                             'shake_data.hdf')
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
    try:
        set_files(event_path, {'event.xml': 'event.xml'})
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        # Should fail validation
        set_files(event_path, {'model_select.conf.bad0': 'model_select.conf'})
        augment = AugmentModule('nc72282711', comment='Test comment.')
        with pytest.raises(RuntimeError):
            augment.execute()
        clear_files(event_path)

        # Should fail vs30 filename check
        set_files(event_path, {'event.xml': 'event.xml'})
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        set_files(event_path, {'event.xml': 'event.xml',
                               'model_select.conf.bad1': 'model_select.conf'})
        augment = AugmentModule('nc72282711', comment='Test comment.')
        with pytest.raises(RuntimeError):
            augment.execute()
        clear_files(event_path)

        # Should fail prediction locations filename check
        set_files(event_path, {'event.xml': 'event.xml'})
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        set_files(event_path, {'event.xml': 'event.xml',
                               'model_select.conf.bad2': 'model_select.conf'})
        augment = AugmentModule('nc72282711', comment='Test comment.')
        with pytest.raises(FileNotFoundError):
            augment.execute()
    finally:
        clear_files(event_path)

    #
    # Switch originators (should succeed)
    #
    set_files(event_path, {'event.xml': 'event.xml',
                           'model.conf.cz': 'model.conf'})
    try:
        assemble = AssembleModule('nc72282711', comment='Test comment.')
        assemble.execute()
        augment = AugmentModule('nc72282711', comment='Test comment.')
        augment.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)
        clear_files(event_path)


def test_assemble_augment_command_line():
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
    test_assemble_augment()
    test_assemble_augment_command_line()
