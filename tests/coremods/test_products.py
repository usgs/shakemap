#!/usr/bin/env python
import os
import os.path

import matplotlib
matplotlib.use('Agg')

import pytest

from shakemap.utils.config import get_config_paths
from shakemap.coremods.model import ModelModule
from shakemap.coremods.assemble import AssembleModule
from shakemap.coremods.contour import ContourModule
from shakemap.coremods.gridxml import GridXMLModule
from shakemap.coremods.info import InfoModule
from shakemap.coremods.raster import RasterModule
from shakemap.coremods.rupture import RuptureModule
from shakemap.coremods.stations import StationModule
from shakemap.coremods.mapping import MappingModule

########################################################################
# Test sm_model
########################################################################
def test_products():

    installpath, datapath = get_config_paths()
    try:
        #
        # Make sure an output file exists
        #
        assemble = AssembleModule('nc72282711')
        assemble.execute()
        model = ModelModule('nc72282711')
        model.execute()

        #
        # Test the creation of products -- currently not checking results
        # for validity or consistency, but probably should
        #
        mod = ContourModule('nc72282711')
        mod.execute()
        mod = GridXMLModule('nc72282711')
        mod.execute()
        mod = InfoModule('nc72282711')
        mod.execute()
        mod = RasterModule('nc72282711')
        mod.execute()
        mod = RuptureModule('nc72282711')
        mod.execute()
        mod = StationModule('nc72282711')
        mod.execute()
        mod = MappingModule('nc72282711')
        mod.execute()
    finally:
        data_file = os.path.join(datapath, 'nc72282711', 'current',
                                 'shake_data.hdf')
        if os.path.isfile(data_file):
            os.remove(data_file)
        res_file = os.path.join(datapath, 'nc72282711', 'current',
                                'products', 'shake_results.hdf')
        if os.path.isfile(res_file):
            os.remove(res_file)

if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_products()
