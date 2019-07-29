#!/usr/bin/env python

# stdlib imports
import datetime
import io
import numpy as np
import random
import string
import sys
import tempfile
import os.path

# third party imports
import pytest

# local imports
from shakelib.utils.containers import ShakeMapInputContainer
from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakelib.rupture.point_rupture import PointRupture


homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
sys.path.insert(0, shakedir)


def randomword(length):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(length))


def dict_equal(d1, d2):
    s1 = sorted(set(d1.keys()))
    s2 = sorted(set(d2.keys()))
    return s1 == s2


def test_input_container():
    f, datafile = tempfile.mkstemp()
    os.close(f)
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9"
time="2008-05-12T06:28:01Z"
depth="19.0" locstring="EASTERN SICHUAN, CHINA" productcode="us2008ryan"
mech="" netid="us" network="" />"""
    try:
        config = {
            'alliance': 'chaotic neutral',
            'race': 'Elf',
            'armor': 5,
            'class': 'Warrior',
            'intelligence': 10
        }
        rupturefile = os.path.join(homedir, 'container_data',
                                   'Barkaetal02_fault.txt')
        eventfile = io.StringIO(event_text)
        datafiles = [os.path.join(
            homedir, 'container_data/northridge_stations_dat.xml')]

        timestamp = datetime.datetime.utcnow().strftime('%FT%TZ')
        originator = 'us'
        version = 1
        history = {'history': [[timestamp, originator, version]]}

        container = ShakeMapInputContainer.createFromInput(
            datafile,
            config,
            eventfile,
            history,
            datafiles=datafiles,
            rupturefile=rupturefile)
        cfile = container.getFileName()
        assert datafile == cfile
        config = container.getConfig()
        station = container.getStationList()
        rupture = container.getRuptureObject()
        history = container.getVersionHistory()
        container.close()

        container2 = ShakeMapInputContainer.load(datafile)
        config2 = container2.getConfig()
        station2 = container2.getStationList()  # noqa
        rupture2 = container2.getRuptureObject()  # noqa
        history2 = container2.getVersionHistory()  # noqa

        assert dict_equal(config, config2)
        df1, _ = station.getStationDictionary(instrumented=False)
        df2, _ = station2.getStationDictionary(instrumented=False)
        assert dict_equal(df1, df2)
        df1, _ = station.getStationDictionary(instrumented=True)
        df2, _ = station2.getStationDictionary(instrumented=True)
        assert dict_equal(df1, df2)
        assert history['history'][-1][0] == history['history'][-1][0]
        assert history['history'][-1][1] == history['history'][-1][1]
        assert history['history'][-1][2] == history['history'][-1][2]

        container2.close()

        eventfile.seek(0)
        container3 = ShakeMapInputContainer.createFromInput(
            datafile,
            config,
            eventfile,
            {})
        try:
            # this should fail, because we haven't set any station data yet
            station = container3.getStationList()
        except AttributeError:
            assert 1 == 1
        rupture = container3.getRuptureObject()
        history = container3.getVersionHistory()
        assert len(history) == 0
        assert isinstance(rupture, PointRupture)

        container3.setStationData(datafiles)

        #
        # Test the getStationDict() and setStationDict() functions with
        # some dummy data
        #
        config = {
            'alliance': 'chaotic neutral',
            'race': 'Elf',
            'armor': 5,
            'class': 'Warrior',
            'intelligence': 10
        }
        with pytest.raises(AttributeError):
            container3.getStationDict()
        with pytest.raises(TypeError):
            container3.setStationDict(None)
        container3.setStationDict(config)
        config2 = container3.getStationDict()
        assert dict_equal(config, config2)
        container3.close()

    except Exception:
        assert 1 == 2
    finally:
        os.remove(datafile)


def test_output_container():
    nrows = 400
    ncols = 230
    # create MMI mean data for maximum component
    mean_mmi_maximum_data = np.random.rand(nrows, ncols)
    mean_mmi_maximum_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Gandalf',
        'color': 'white',
        'powers': 'magic'
    }

    # create MMI std data for maximum component
    std_mmi_maximum_data = mean_mmi_maximum_data / 10
    std_mmi_maximum_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Legolas',
        'color': 'green',
        'powers': 'good hair'
    }

    # create MMI mean data for rotd50 component
    mean_mmi_rotd50_data = np.random.rand(nrows, ncols)
    mean_mmi_rotd50_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Gimli',
        'color': 'brown',
        'powers': 'axing'
    }

    # create MMI std data for rotd50 component
    std_mmi_rotd50_data = mean_mmi_rotd50_data / 10
    std_mmi_rotd50_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Aragorn',
        'color': 'white',
        'powers': 'scruffiness'
    }

    # create PGA mean data for maximum component
    mean_pga_maximum_data = np.random.rand(nrows, ncols)
    mean_pga_maximum_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Pippin',
        'color': 'purple',
        'powers': 'rashness'
    }

    # create PGA std data for maximum component
    std_pga_maximum_data = mean_pga_maximum_data / 10
    std_pga_maximum_metadata = {
        'xmin': -118.5,
        'xmax': -114.5,
        'ymin': 32.1,
        'ymax': 36.7,
        'dx': 0.01,
        'dy': 0.02,
        'nx': 400,
        'ny': 230,
        'name': 'Merry',
        'color': 'grey',
        'powers': 'hunger'
    }

    f, datafile = tempfile.mkstemp()
    os.close(f)

    try:
        container = ShakeMapOutputContainer.create(datafile)
        # LookupError raised if trying to dropIMTs if there are none
        with pytest.raises(LookupError):
            container.dropIMT('mmi')

        # Add imts
        container.setIMTGrids('mmi',
                              mean_mmi_maximum_data, mean_mmi_maximum_metadata,
                              std_mmi_maximum_data, std_mmi_maximum_metadata,
                              component='maximum')
        container.setIMTGrids('mmi',
                              mean_mmi_rotd50_data, mean_mmi_rotd50_metadata,
                              std_mmi_rotd50_data, std_mmi_rotd50_metadata,
                              component='rotd50')
        container.setIMTGrids('pga',
                              mean_pga_maximum_data, mean_pga_maximum_metadata,
                              std_pga_maximum_data, std_pga_maximum_metadata,
                              component='maximum')

        # get the maximum MMI imt data
        mmi_max_dict = container.getIMTGrids('mmi', component='maximum')
        np.testing.assert_array_equal(mmi_max_dict['mean'],
                                      mean_mmi_maximum_data)
        np.testing.assert_array_equal(mmi_max_dict['std'],
                                      std_mmi_maximum_data)
        assert mmi_max_dict['mean_metadata'] == mean_mmi_maximum_metadata
        assert mmi_max_dict['std_metadata'] == std_mmi_maximum_metadata

        # get the rotd50 MMI imt data
        mmi_rot_dict = container.getIMTGrids('mmi', component='rotd50')
        np.testing.assert_array_equal(mmi_rot_dict['mean'],
                                      mean_mmi_rotd50_data)
        np.testing.assert_array_equal(mmi_rot_dict['std'],
                                      std_mmi_rotd50_data)
        assert mmi_rot_dict['mean_metadata'] == mean_mmi_rotd50_metadata
        assert mmi_rot_dict['std_metadata'] == std_mmi_rotd50_metadata

        # get list of all imts
        imts = container.getIMTs()

        # get list of maximum imts
        max_imts = container.getIMTs(component='maximum')
        assert sorted(max_imts) == ['mmi', 'pga']

        # get list of components for mmi
        mmi_comps = container.getComponents('mmi')
        assert sorted(mmi_comps) == ['maximum', 'rotd50']

        # Test dropIMT
        imts = container.getIMTs('maximum')
        assert imts == ['mmi', 'pga']
        container.dropIMT('mmi')
        imts = container.getIMTs('maximum')
        assert imts == ['pga']
        container.close()

    except Exception as e:
        raise(e)
    finally:
        os.remove(datafile)


def test_output_arrays():

    f, datafile = tempfile.mkstemp()
    os.close(f)

    try:
        container = ShakeMapOutputContainer.create(datafile)
        #
        # Test that no data type is set
        #
        assert container.getDataType() is None

        #
        # Make some array data and metadata
        #
        mean = np.random.rand(100)
        std = np.random.rand(100)
        lats = np.random.rand(100)
        lons = np.random.rand(100)
        ids = np.array([randomword(4).encode('ascii') for x in range(100)])
        metadata = {
            'units': '%g',
            'digits': 4
        }
        #
        # Put the data in the container
        #
        container.setIMTArrays('PGA', lons, lats, ids,
                               mean, metadata,
                               std, metadata, 'Larger')
        #
        # Now extract it and compare it to what we put in there
        #
        dout = container.getIMTArrays('PGA', 'Larger')
        assert all(dout['lons'] == lons)
        assert all(dout['lats'] == lats)
        assert all(dout['ids'] == ids)
        assert all(dout['mean'] == mean)
        assert all(dout['std'] == std)
        #
        # Check the data type
        #
        assert container.getDataType() == 'points'
        #
        # Try raising some exceptions
        #
        # Shouldn't be able to find this IMT
        with pytest.raises(LookupError):
            container.getIMTArrays('JUNK', 'Larger')
        # Shapes of inputs not the same
        with pytest.raises(ValueError):
            empty = np.array([])
            container.setIMTArrays('PGV', empty, lats, ids,
                                   mean, metadata,
                                   std, metadata, 'Larger')
        # IMT already exists
        with pytest.raises(LookupError):
            container.setIMTArrays('PGA', lons, lats, ids,
                                   mean, metadata,
                                   std, metadata, 'Larger')
        # Trying to set a grid in a file with points
        with pytest.raises(TypeError):
            container.setIMTGrids('PGV', mean, metadata,
                                  std, metadata, 'Larger')
        # Trying to get a grid in a file with points
        with pytest.raises(TypeError):
            container.getIMTGrids('PGA', 'Larger')

        container.close()

    except Exception as e:
        if os.path.isfile(datafile):
            os.remove(datafile)
        raise(e)
    finally:
        if os.path.isfile(datafile):
            os.remove(datafile)


if __name__ == '__main__':
    test_input_container()
    test_output_container()
    test_output_arrays()
