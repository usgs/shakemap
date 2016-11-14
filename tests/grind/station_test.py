#!/usr/bin/env python

# stdlib modules
import sys
import os.path
import time as time
import shutil

# third party modules
import pandas.util.testing as pdt

# local imports
from shakemap.grind.station import StationList
from shakemap.grind.source import Source
from shakemap.grind.multigmpe import MultiGMPE
from shakemap.grind.sites import Sites
from shakemap.grind.gmice.wgrw12 import WGRW12
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.gsim.allen_2012_ipe import AllenEtAl2012

homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
sys.path.insert(0, shakedir)


def test_station(tmpdir):

    homedir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.join(homedir, '..', 'data', 
            'eventdata', 'Calexico', 'input'))

    #
    # Read the event, source, and rupture files and produce a Source object
    #
    inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    dyfifile = os.path.join(datadir, 'ciim3_dat.xml')
    eventfile = os.path.join(datadir, 'event.xml')
    rupturefile = os.path.join(datadir, 'wei_fault.txt')

    source_obj = Source.fromFile(eventfile, rupturefile=rupturefile)

    #
    # Set up the GMPE, IPE, and GMICE
    #
    gmpe_cy14 = ChiouYoungs2014()

    gmpe = MultiGMPE.from_list([gmpe_cy14], [1.0])

    gmice = WGRW12()

    ipe = AllenEtAl2012()

    #
    # 
    #
    rupture_ctx = source_obj.getRuptureContext([gmpe])

    smdx = 0.0083333333
    smdy = 0.0083333333
    lonspan = 6.0
    latspan = 4.0
    vs30filename = os.path.join(datadir, '..', 'vs30', 'vs30.grd')

    sites_obj_grid = Sites.fromCenter(
            rupture_ctx.hypo_lon, rupture_ctx.hypo_lat, lonspan, latspan, 
            smdx, smdy, defaultVs30=760.0, vs30File=vs30filename, 
            vs30measured_grid=None, padding=False, resample=False
        )

    xmlfiles = [inputfile, dyfifile]
#    dbfile = str(tmpdir.join('stations.db'))
    dbfile = os.path.join(str(tmpdir), 'stations.db')

    stations = StationList.fromXML(xmlfiles, dbfile, source_obj, 
            sites_obj_grid, gmpe, ipe, gmice)

    df1 = stations.getStationDataframe(1, sort=True)
    df2 = stations.getStationDataframe(0, sort=True)


    #
    # In case the test starts failing because of some minor change
    # in one of the prediction or conversion equations (or roundoff
    # or whatever), but the code is running correctly, uncomment 
    # these lines and re-run the test. Then, copy the new stations.db
    # file into tests/data/eventdata/Calexico/database/. Then
    # recomment these lines and rerun the test. It should succeed.
    #
    #shutil.copy(dbfile,'./stations.db')
    #print(os.getcwd())

    #
    # We should probably check these dataframes against some established
    # set, and also check the database against a known database. 
    #

    ref_dbfile = os.path.join(datadir, '..', 'database', 'stations.db')

    stations2 = StationList(ref_dbfile)

    ref_df1 = stations2.getStationDataframe(1, sort=True)
    ref_df2 = stations2.getStationDataframe(0, sort=True)

#    assert ref_df1.equals(df1)
#    assert ref_df2.equals(df2)

    pdt.assert_frame_equal(df1, ref_df1)
    pdt.assert_frame_equal(df2, ref_df2)

if __name__ == '__main__':
    import tempfile

    td = tempfile.TemporaryDirectory()
    test_station(td.name)

