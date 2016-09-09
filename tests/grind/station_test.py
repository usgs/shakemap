#!/usr/bin/env python

# stdlib modules
import sys
import os.path
#import tempfile
import time as time

# third party modules

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
# put this at the front of the system path, ignoring any installed
# shakemap stuff
sys.path.insert(0, shakedir)

# local imports
from shakemap.grind.station import StationList
from shakemap.grind.source import Source
from shakemap.grind.multigmpe import MultiGMPE
from shakemap.grind.sites import Sites
from shakemap.grind.gmice.wgrw12 import WGRW12
from openquake.hazardlib.gsim.chiou_youngs_2014 import ChiouYoungs2014
from openquake.hazardlib.gsim.allen_2012_ipe import AllenEtAl2012


def test_station(tmpdir):

    homedir = os.path.dirname(os.path.abspath(__file__))
    datadir = os.path.abspath(os.path.join(homedir, '..', 'data', 
            'eventdata', 'Calexico', 'input'))

    #
    # Read the event, source, and fault files and produce a Source object
    #
    inputfile = os.path.join(datadir, 'stationlist_dat.xml')
    dyfifile = os.path.join(datadir, 'ciim3_dat.xml')
    eventfile = os.path.join(datadir, 'event.xml')
    faultfile = os.path.join(datadir, 'wei_fault.txt')

    source_obj = Source.readFromFile(eventfile, faultfile=faultfile)

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

    sites_obj_grid = Sites.createFromCenter(
            rupture_ctx.hypo_lon, rupture_ctx.hypo_lat, lonspan, latspan, 
            smdx, smdy, defaultVs30=760.0, vs30File=vs30filename, 
            vs30measured_grid=None, padding=False, resample=False
        )

    xmlfiles = [inputfile, dyfifile]
    dbfile = str(tmpdir.join('stations.db'))

    stations = StationList.fromXML(xmlfiles, dbfile, source_obj, 
            sites_obj_grid, gmpe, ipe, gmice)

    df1 = stations.getStationDataframe(1, sort=True)

    df2 = stations.getStationDataframe(0, sort=False)

    #
    # We should probably check these dataframes against some established
    # set, and also check the database against a known database. 
    #

