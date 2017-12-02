#!/usr/bin/env python

import os
import shutil

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.assemble import AssembleModule


########################################################################
# Test sm_assemble
########################################################################
def test_assemble():

    installpath, datapath = get_config_paths()

    # Process a non-existent event (should fail)
    amod = AssembleModule('not_an_event')
    with pytest.raises(NotADirectoryError):
        amod.execute()

    # Would succeed but we remove event.xml (should fail)
    event_file = os.path.join(datapath, 'wenchuan', 'current', 'event.xml')
    os.rename(event_file, event_file + '_safe')
    try:
        amod = AssembleModule('wenchuan')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        os.rename(event_file + '_safe', event_file)

    # Normal event (should succeed)
    data_file = os.path.join(datapath, 'wenchuan', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        amod = AssembleModule('wenchuan')
        amod.execute()
        #
        # Run a second time to exercise a different branch of the code
        #
        amod.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    # Do an event with model.conf (not model_zc.conf) and no zoneinfo
    # (should succeed)
    data_file = os.path.join(datapath, 'nc72282711', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        amod = AssembleModule('nc72282711')
        amod.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

    # Try not having an event-specific config (should succeed)
    model_file = os.path.join(datapath, 'nc72282711', 'current',
                              'model.conf')
    os.rename(model_file, model_file + '_safe')
    data_file = os.path.join(datapath, 'nc72282711', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        amod = AssembleModule('nc72282711')
        amod.execute()
    finally:
        os.rename(model_file + '_safe', model_file)
        if os.path.isfile(data_file):
            os.remove(data_file)

    # Do an event with DYFI data (should succeed)
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current',
                            'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.rename(hdf_file, hdf_file + '_safe')
    try:
        amod = AssembleModule('nc72282711_dyfi')
        amod.execute()
    finally:
        if os.path.isfile(hdf_file + '_safe'):
            os.rename(hdf_file + '_safe', hdf_file)
    
    #
    # Try some bad config files
    #
    # Should fail validation
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad0', model_file)
    try:
        amod = AssembleModule('nc72282711_nodata_nofault')
        with pytest.raises(RuntimeError):
            amod.execute()
    finally:
        os.rename(model_file + '_safe', model_file)

    # Should fail vs30 filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad1', model_file)
    try:
        amod = AssembleModule('nc72282711_nodata_nofault')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        os.rename(model_file + '_safe', model_file)

    # Should fail prediction locations filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad2', model_file)
    try:
        amod = AssembleModule('nc72282711_nodata_nofault')
        with pytest.raises(FileNotFoundError):
            amod.execute()
    finally:
        os.rename(model_file + '_safe', model_file)
    #
    # Make sure the location file substitutions work (should succeed)
    #
    data_file = os.path.join(datapath, 'northridge_points', 'current', 'shake_data.hdf')
    if os.path.isfile(data_file):
        os.remove(data_file)
    try:
        amod = AssembleModule('northridge_points')
        amod.execute()
    finally:
        if os.path.isfile(data_file):
            os.remove(data_file)

if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_assemble()
