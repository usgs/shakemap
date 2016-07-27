#!/usr/bin/env python

#stdlib modules
import sys
import os.path

import tempfile
import time
from datetime import datetime

#third party modules

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..','..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed shakemap stuff

#local imports
from shakemap.grind.station import StationList
from shakemap.grind.source import Source

def _test():

    tmp,dbfile = tempfile.mkstemp()
    os.close(tmp)
    os.remove(dbfile)

    homedir = os.path.dirname(os.path.abspath(__file__))
    xmlfile = os.path.abspath(os.path.join(homedir,'..','data','eventdata',
                                           'northridge','northridge_stations.xml'))
    stationfile = os.path.abspath(os.path.join(homedir,'..','data','eventdata',
                                               'northridge','northridge_stations.db'))
    eventdict = {'lat':34.213,'lon':-118.537,'depth':18.2,
                 'mag':6.7,'time':datetime(1994,1,17,12,30,55),
                 'mech':'ALL','dip':45,'rake':90}

    try:
        print('Testing load from XML format...')
        t1 = time.time()
        stations1 = StationList.loadFromXML([xmlfile],dbfile)
        t2 = time.time()
        print('Passed load from XML format %i stations in %.2f seconds.' % (len(stations1),t2-t1))

        print('Testing filling in distance and derived MMI/PGM values...')
        source = Source(eventdict)
        stations1.fillTables(source)
        print('Passed filling in distance and derived MMI/PGM values...')
        
        print('Testing retrieval of MMI data from StationList object...')
        t1 = time.time()
        mmidf1 = stations1.getMMIStations()
        t2 = time.time()
        print('Passed retrieval of %i MMI data in %.2f seconds from StationList object.' % (len(mmidf1),t2-t1))

        print('Testing retrieval of instrumented data from StationList object...')
        t1 = time.time()
        imtdf1 = stations1.getInstrumentedStations()
        t2 = time.time()
        print('Passed retrieval of %i instrumented data in %.2f seconds from StationList object.' % (len(imtdf1),t2-t1))


        print('Testing load from sqlite format...')
        t1 = time.time()
        stations2 = StationList(stationfile)
        t2 = time.time()
        print('Passed load from sqlite format %i stations in %.2f seconds.' % (len(stations1),t2-t1))

        print('Testing retrieval of MMI data from StationList object...')
        t1 = time.time()
        mmidf2 = stations2.getMMIStations()
        t2 = time.time()
        print('Passed retrieval of %i MMI data in %.2f seconds from StationList object.' % (len(mmidf2),t2-t1))

        print('Testing retrieval of instrumented data from StationList object...')
        t1 = time.time()
        imtdf2 = stations2.getInstrumentedStations()
        t2 = time.time()
        print('Passed retrieval of %i instrumented data in %.2f seconds from StationList object.' % (len(imtdf1),t2-t1))

        assert(len(stations1) == len(stations2))

        
    except Exception as msg:
        print('Error caught: %s' % str(msg))
    if os.path.isfile(dbfile):
        os.remove(dbfile)
    
