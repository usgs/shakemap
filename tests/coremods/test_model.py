import os
import shutil

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule

########################################################################
# Test sm_model
########################################################################
def test_model():

    installpath, datapath = get_config_paths()
    #
    # This is Northridge for a set of output points (not a grid)
    # Remove the products directory to hit that code
    # (should succeed)
    #
    assemble = AssembleModule('northridge_points')
    assemble.execute()
    products_dir = os.path.join(datapath, 'northridge_points', 'current', 
                            'products')
    if os.path.isdir(products_dir):
        shutil.rmtree(products_dir)
    model = ModelModule('northridge_points')
    model.execute()

    #
    # This is a small grid with station data only (should succeed)
    #
    assemble = AssembleModule('nc72282711')
    assemble.execute()
    model = ModelModule('nc72282711')
    model.execute()

    #
    # This is a small grid with DYFI data only (should succeed)
    #
    assemble = AssembleModule('nc72282711_dyfi')
    assemble.execute()
    model = ModelModule('nc72282711_dyfi')
    model.execute()

    #
    # Run with no data and no fault, and use the default extent.
    #
    assemble = AssembleModule('nc72282711_nodata_nofault')
    assemble.execute()
    model = ModelModule('nc72282711_nodata_nofault')
    model.execute()

    #
    # Set the bias and outlier magnitude limits low to test additional
    # code branches
    #
    assemble = AssembleModule('nc72282711_nofault')
    assemble.execute()
    model = ModelModule('nc72282711_nofault')
    model.execute()

    #
    # This event exists, but we hide the hdf file (should fail)
    #
    hdf_file = os.path.join(datapath, 'nc72282711_dyfi', 'current', 
                            'shake_data.hdf')
    os.rename(hdf_file, hdf_file + '_safe')
    model = ModelModule('nc72282711_dyfi')
    with pytest.raises(FileNotFoundError):
        model.execute()
    os.rename(hdf_file + '_safe', hdf_file)

    #
    # This event doesn't exist (should fail)
    #
    model = ModelModule('not_an_event')
    with pytest.raises(NotADirectoryError):
        model.execute()

