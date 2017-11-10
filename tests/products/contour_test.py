#!/usr/bin/env python

#stdlib imports
import os.path
import tempfile
import shutil

#third party imports
from configobj import ConfigObj
import numpy as np

#neic imports
from shakemap.utils.config import get_config_paths
from shakelib.utils.containers import OutputContainer
from shakemap.products.contour import (contour,
                                       contour_to_files,
                                       _get_default_intervals)

def test_intervals():
    fgrid = np.linspace(0.007,11.7,50)
    intervals = _get_default_intervals(fgrid)
    np.testing.assert_almost_equal(intervals[0],0.01)
    np.testing.assert_almost_equal(intervals[-1],10)

    fgrid = np.linspace(0.1,70,50)
    intervals = _get_default_intervals(fgrid)
    np.testing.assert_almost_equal(intervals[0],0.3)
    np.testing.assert_almost_equal(intervals[-1],30)

    fgrid = np.array([3.5,5.6,8.2])
    intervals = _get_default_intervals(fgrid,interval_type='linear')
    np.testing.assert_almost_equal(intervals[0],4.0)
    np.testing.assert_almost_equal(intervals[-1],8.0)


def test_contour():
    install_path, data_path = get_config_paths()
    homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
    datafile = os.path.join(homedir,'..','data','shake_result.hdf')
    container = OutputContainer.load(datafile)
    config_file = os.path.join(install_path, 'config', 'products.conf')
    config = ConfigObj(config_file)
    line_strings = contour(container,'MMI','Larger',None)
    assert line_strings[0]['properties']['value'] == 4.0
    assert line_strings[0]['geometry']['type'] == 'MultiLineString'
    assert len(line_strings[0]['geometry']['coordinates'][0]) == 81

    
#TODO - test with variety of shakemaps, incl those spanning 180 meridian.
# def test_contour_files():
#     outfolder = tempfile.mkdtemp()
#     try:
#         install_path, data_path = get_config_paths()
#         homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
#         datafile = os.path.join(homedir,'..','data','shake_result.hdf')
#         container = OutputContainer.load(datafile)
#         config_file = os.path.join(install_path, 'config', 'products.conf')
#         config = ConfigObj(config_file)
#         contour_to_files(container,config,outfolder,None)
#         jsonfiles = os.listdir(outfolder)
#         cmpfiles = ['us19940117123055_MMI_Larger.json',
#                     'us19940117123055_PGA_Larger.json',
#                     'us19940117123055_PGV_Larger.json',
#                     'us19940117123055_PSA0p3_Larger.json',
#                     'us19940117123055_PSA1p0_Larger.json',
#                     'us19940117123055_PSA3p0_Larger.json']
#         assert sorted(jsonfiles) == cmpfiles
#     except Exception as e:
#         raise(e)
#     finally:
#         shutil.rmtree(outfolder) #remove outfolder and all files in it
                    

if __name__ == '__main__':
    test_intervals()
    test_contour()
    #test_contour_files()
