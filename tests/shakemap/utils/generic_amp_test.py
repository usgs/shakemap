#!/usr/bin/env python

import os
import shutil

import numpy as np

from shakemap.utils.generic_amp import get_generic_amp_factors
from shakemap.utils.config import get_config_paths
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from mapio.gridcontainer import GridHDFContainer

do_test = True


class Dummy(object):
    pass


def make_generic_amps():
    imts = ['PGA', 'SA(0.3)', 'SA(1.0)', 'SA(3.0)']
    install_path, _ = get_config_paths()
    geodict = {'dx': 0.016666666666666666,
               'dy': 0.016666666666666666,
               'nx': 301,
               'ny': 151,
               'xmax': -116.0,
               'xmin': -121.0,
               'ymax': 35.5,
               'ymin': 33.0}
    gd = GeoDict(geodict)

    # make east-west file (1s on the left, 0s on the right)
    data = np.ones((gd.ny, gd.nx))
    data[:, 151:] = 0
    outfolder = os.path.join(install_path, 'data', 'GenericAmpFactors')
    east_west_file = os.path.join(outfolder, 'Test_basin_east_west.hdf')
    east_west = GridHDFContainer.create(east_west_file)
    for imt in imts:
        grid = Grid2D(data, gd)
        east_west.setGrid(imt, grid)
    east_west.close()

    # make east-west file (1s on the left, 0s on the right)
    data = np.ones((gd.ny, gd.nx))
    data[76:151, :] = 0
    outfolder = os.path.join(install_path, 'data', 'GenericAmpFactors')
    north_south_file = os.path.join(outfolder, 'Test_basin_north_south.hdf')
    north_south = GridHDFContainer.create(north_south_file)
    for imt in imts:
        grid = Grid2D(data, gd)
        north_south.setGrid(imt, grid)
    north_south.close()

    return (east_west_file, north_south_file)


def test_generic_amp():
    try:
        #
        # Using the LA basin test data set, make a set of lons and lats
        # Bounds of data set: -119.284/-117.284/32.9/34.6
        # But we want to go somewhat beyond the bounds in order to test
        # that functionality
        #
        lons = np.linspace(-121.0, -118.5, 25)
        lons = np.append(lons, np.linspace(-118.5, -116.0, 25))
        lats = np.linspace(33.0, 35.8, 25)
        lats = np.append(lats, np.linspace(33.0, 35.8, 25))

        sx = Dummy()
        sx.lons = lons
        sx.lats = lats

        install_path, _ = get_config_paths()
        outfolder = os.path.join(install_path, 'data', 'GenericAmpFactors')
        tmpout = os.path.join(install_path, 'data', 'ga_tmp')
        shutil.move(outfolder, tmpout)
        gaf = get_generic_amp_factors(sx, 'PGA')
        assert gaf is None
        shutil.move(tmpout, outfolder)
        gaf = get_generic_amp_factors(sx, 'PGA')
        assert gaf is None

        east_west_file, north_south_file = make_generic_amps()

        gaf = get_generic_amp_factors(sx, 'PGA')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('PGA:')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'PGV')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('PGV:')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'SA(0.3)')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('SA(0.3):')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'SA(2.0)')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('SA(2.0):')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'SA(4.0)')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('SA(4.0):')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'SA(0.1)')
        gaf_target = np.array(
          [1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  2.,  2.,
           2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  0.,  0.,  0.,  1.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  1.,  1.,
           1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  0.,  0.,  0.])
        if do_test is True:
            assert np.allclose(gaf, gaf_target)
        else:
            print('SA(0.1):')
            print(repr(gaf))

        gaf = get_generic_amp_factors(sx, 'XXXX')
        gaf_target = np.array(
          [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
           0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])
        assert np.allclose(gaf, gaf_target)

    except Exception:
        assert 1 == 2
    finally:
        os.remove(east_west_file)
        os.remove(north_south_file)


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_generic_amp()
