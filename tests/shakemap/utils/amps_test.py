#!/usr/bin/env python

import os.path
from datetime import datetime

import numpy as np

from shakemap.utils.config import get_config_paths
from shakemap.utils.amps import AmplitudeHandler, dt_to_timestamp
from shakemap.coremods.select import validate_config
from configobj import ConfigObj


def test_amps():
    try:
        install_path, data_path = get_config_paths()

        # dbfile location
        homedir = os.path.dirname(os.path.abspath(__file__))
        dbfile = os.path.join(homedir,'..','..','data','install','data','amps.db')

        if os.path.isfile(dbfile):
            os.remove(dbfile)
        handler = AmplitudeHandler(install_path,data_path)

        # test inserting events into the database
        event = {'id':'ci37889959',
                 'time':datetime(2018,3,7,18,5,0 ),
                 'lat':35.487,
                 'lon':-120.027,
                 'depth':8.0,
                 'mag':3.7}
        handler.insertEvent(event)
        info = handler.getStats()
        assert info['events'] == 1
        
        homedir = os.path.dirname(os.path.abspath(__file__))
        xmlfile = os.path.join(homedir,'..','..',
                               'data','ampdata','USR_100416_20180307_180450.xml')

        handler.insertAmps(xmlfile)
        info = handler.getStats()
        assert info['stations'] == 1
        assert info['stations'] == 1
        assert info['station_min'] == datetime(2018, 3, 7, 18, 4, 49)
        assert info['station_max'] == datetime(2018, 3, 7, 18, 4, 49)
        assert info['channels'] == 3
        assert info['pgms'] == 15
        eqtime = dt_to_timestamp(event['time'])
        eqlat = event['lat']
        eqlon = event['lon']
        df = handler.associate(eqtime,eqlat,eqlon)
        df_pga = df[df['imt']=='pga']
        np.testing.assert_almost_equal(df_pga['value'].sum(),0.10420000000000001)

        del handler
        os.remove(dbfile)

        # test global associator
        handler = AmplitudeHandler(install_path,data_path)
        handler.insertEvent(event)
        handler.insertAmps(xmlfile)
        nassociated = handler.associateAll(write_data=False)
        assert nassociated == 1

        del handler
        os.remove(dbfile)

        #test clean methods
        handler = AmplitudeHandler(install_path,data_path)
        handler.insertEvent(event)
        handler.insertAmps(xmlfile)
        handler.cleanEvents(threshold=1)
        handler.cleanAmps(threshold=1)
        info = handler.getStats()
        assert info['events'] == 0
        assert info['stations'] == 0
        assert info['channels'] == 0
        assert info['pgms'] == 0
    except:
        assert 1==2
    finally:
        if os.path.isfile(dbfile):
            os.remove(dbfile)
    
if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_amps()
