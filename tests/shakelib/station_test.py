#!/usr/bin/env python

# stdlib modules
import sys
import os.path
import pickle

# third party modules
import numpy as np

# local imports
from shakelib.station import StationList

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)

#
# Set SAVE to True to write new versions of the output to disk,
# set it to False to actually run the tests.
#
SAVE = False


def test_station():

    homedir = os.path.dirname(os.path.abspath(__file__))

    #
    # First test the Calexico data on its own
    #
    event = 'Calexico'

    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    dyfifile = os.path.join(datadir, 'ciim3_dat.xml')
    xmlfiles = [inputfile, dyfifile]

    stations = StationList.loadFromXML(xmlfiles, ":memory:")

    df1 = stations.getStationDictionary(instrumented=True)
    df2 = stations.getStationDictionary(instrumented=False)

    ppath = os.path.abspath(os.path.join(datadir, '..', 'database',
                                         'test1.pickle'))
    if SAVE:
        ldf = [df1, df2]
        with open(ppath, 'wb') as f:
            pickle.dump(ldf, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open(ppath, 'rb') as f:
            ldf = pickle.load(f)

        saved_df1 = ldf[0]
        saved_df2 = ldf[1]

        compare_dataframes(saved_df1, df1)
        compare_dataframes(saved_df2, df2)

    #
    # Should at least hit this code
    #
    imtlist = stations.getIMTtypes()
    assert 'PGA' in imtlist
    assert 'PGV' in imtlist

    #
    # Add the Northridge data to the Calexico data to test
    # addData()
    #
    event = 'northridge'
    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    inputfile = os.path.join(datadir, 'hist_dat.xml')
    dyfifile = os.path.join(datadir, 'dyfi_dat.xml')
    xmlfiles = [inputfile, dyfifile]

    stations = stations.addData(xmlfiles)

    df1 = stations.getStationDictionary(instrumented=True)
    df2 = stations.getStationDictionary(instrumented=False)

    ppath = os.path.abspath(os.path.join(datadir, '..', 'database',
                                         'test2.pickle'))
    if SAVE:
        ldf = [df1, df2]
        with open(ppath, 'wb') as f:
            pickle.dump(ldf, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open(ppath, 'rb') as f:
            ldf = pickle.load(f)

        saved_df1 = ldf[0]
        saved_df2 = ldf[1]

        compare_dataframes(saved_df1, df1)
        compare_dataframes(saved_df2, df2)


def test_station2():

    #
    # Test the wenchuan data
    #
    homedir = os.path.dirname(os.path.abspath(__file__))

    event = 'wenchuan'

    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    inputfile = os.path.join(datadir, 'stationlist.xml')
    xmlfiles = [inputfile]

    stations = StationList.loadFromXML(xmlfiles, ":memory:")

    df1 = stations.getStationDictionary(instrumented=True)
    df2 = stations.getStationDictionary(instrumented=False)

    ppath = os.path.abspath(os.path.join(datadir, '..', 'database',
                                         'test3.pickle'))
    if SAVE:
        ldf = [df1, df2]
        with open(ppath, 'wb') as f:
            pickle.dump(ldf, f, pickle.HIGHEST_PROTOCOL)
    else:
        with open(ppath, 'rb') as f:
            ldf = pickle.load(f)

        saved_df1 = ldf[0]
        saved_df2 = ldf[1]

        compare_dataframes(saved_df1, df1)
        compare_dataframes(saved_df2, df2)

        #
        # Dump the database to SQL and then restore it to a new
        # StationList object. Compare dataframes.
        #
        sql = stations.dumpToSQL()

        stations2 = StationList.loadFromSQL(sql)

        df1 = stations2.getStationDictionary(instrumented=True)
        df2 = stations2.getStationDictionary(instrumented=False)

        compare_dataframes(saved_df1, df1)
        compare_dataframes(saved_df2, df2)


def test_station3():

    #
    # Exercise the geojson code. Can't really compare it to anything
    # because the ordering is scrambled by the hashes in the XML
    # parsing stage. Once (if) we institute a loadFromJSON() method, we
    # can do a comparison.
    #
    homedir = os.path.dirname(os.path.abspath(__file__))

    event = 'wenchuan'

    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    inputfile = os.path.join(datadir, 'stationlist.xml')
    xmlfiles = [inputfile]

    stations = StationList.loadFromXML(xmlfiles, ":memory:")

    myjson = stations.getGeoJson()  # noqa


def test_station4():

    homedir = os.path.dirname(os.path.abspath(__file__))

    event = 'northridge'
    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    dyfifile = os.path.join(datadir, 'dyfi_dat.xml')
    xmlfiles = [dyfifile]

    stations = StationList.loadFromXML(xmlfiles, ":memory:")

    df1 = stations.getStationDictionary(instrumented=True)  # noqa
    df2 = stations.getStationDictionary(instrumented=False)  # noqa
    assert df1 is None


def test_station5():

    homedir = os.path.dirname(os.path.abspath(__file__))

    event = 'Calexico'

    datadir = os.path.abspath(os.path.join(homedir, 'station_data'))
    datadir = os.path.abspath(os.path.join(datadir, event, 'input'))

    inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    dyfifile = os.path.join(datadir, 'ciim3_dat.xml')

    xmlfiles = [inputfile, dyfifile]
    stations1 = StationList.loadFromXML(xmlfiles, ":memory:")
    #
    # Load the data more than once to exercise the code that handles
    # repeated entries.
    #
    xmlfiles = [inputfile, inputfile, dyfifile, dyfifile]
    stations2 = StationList.loadFromXML(xmlfiles, ":memory:")

    df1 = stations1.getStationDictionary(instrumented=True)
    df2 = stations2.getStationDictionary(instrumented=True)

    compare_dataframes(df1, df2)


def compare_dataframes(df1, df2):

    assert sorted(list(df1.keys())) == sorted(list(df2.keys()))

    idx1 = np.argsort(df1['id'])
    idx2 = np.argsort(df2['id'])

    for key in df1.keys():
        if df1[key].dtype == np.float:
            assert np.allclose(df1[key][idx1], df2[key][idx2], equal_nan=True)
        else:
            assert (df1[key][idx1] == df2[key][idx2]).all()


if __name__ == '__main__':
    test_station()
    test_station2()
    test_station3()
    test_station4()
    test_station5()
