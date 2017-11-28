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
    amod = AssembleModule('wenchuan')
    with pytest.raises(FileNotFoundError):
        amod.execute()
    os.rename(event_file + '_safe', event_file)

    # Normal event (should succeed)
    amod = AssembleModule('wenchuan')
    amod.execute()

    # Do an event with model.conf (not model_zc.conf) and no zoneinfo
    # (should succeed)
    amod = AssembleModule('nc72282711')
    amod.execute()

    # Try not having an event-specific config (should succeed)
    model_file = os.path.join(datapath, 'nc72282711', 'current',
                              'model.conf')
    os.rename(model_file, model_file + '_safe')
    amod = AssembleModule('nc72282711')
    amod.execute()
    os.rename(model_file + '_safe', model_file)

    # Remove the existing hdf file from a no-backup event (should succeed)
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current',
                            'shake_data.hdf')
    if os.path.isfile(hdf_file):
        os.rename(hdf_file, hdf_file + '_safe')
    amod = AssembleModule('nc72282711_dyfi')
    amod.execute()
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
    amod = AssembleModule('nc72282711_nodata_nofault')
    with pytest.raises(RuntimeError):
        amod.execute()
    os.rename(model_file + '_safe', model_file)

    # Should fail vs30 filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad1', model_file)
    amod = AssembleModule('nc72282711_nodata_nofault')
    with pytest.raises(FileNotFoundError):
        amod.execute()
    os.rename(model_file + '_safe', model_file)

    # Should fail prediction locations filename check
    model_file = os.path.join(datapath, 'nc72282711_nodata_nofault', 
                              'current', 'model_zc.conf')
    os.rename(model_file, model_file + '_safe')
    shutil.copyfile(model_file + '.bad2', model_file)
    amod = AssembleModule('nc72282711_nodata_nofault')
    with pytest.raises(FileNotFoundError):
        amod.execute()
    os.rename(model_file + '_safe', model_file)
    #
    # Make sure the location file substitutions work (should succeed)
    #
    amod = AssembleModule('northridge_points')
    amod.execute()

