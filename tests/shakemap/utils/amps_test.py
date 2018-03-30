#!/usr/bin/env python

import os.path
from datetime import datetime

import numpy as np

from shakemap.utils.config import get_config_paths
from shakemap.utils.amps import AmplitudeHandler
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
        event = {'id':'nc72981481',
                 'time':datetime(2018,3,9,6,1,28),
                 'lat':40.2918,
                 'lon':-124.5423,
                 'depth':15.1,
                 'mag':4.5}
        handler.insertEvent(event)
        info = handler.getStats()
        assert info['events'] == 1
        
        homedir = os.path.dirname(os.path.abspath(__file__))
        xmlfile = os.path.join(homedir,'..','..',
                               'data','ampdata','USR_nw0147_20180309_060100.xml')

        handler.insertAmps(xmlfile)
        info = handler.getStats()
        assert info['stations'] == 1
        assert info['stations'] == 1
        assert info['station_min'] == datetime(2018, 3, 9, 13, 1, 20)
        assert info['station_max'] == datetime(2018, 3, 9, 13, 1, 20)
        assert info['channels'] == 3
        assert info['pgms'] == 15
        eqtime = int(event['time'].timestamp())
        eqlat = event['lat']
        eqlon = event['lon']
        df = handler.associate(eqtime,eqlat,eqlon)
        df_pga = df[df['imt']=='pga']
        np.testing.assert_almost_equal(df_pga['value'].sum(),0.34539999999999998)

        del handler
        os.remove(dbfile)
        
        handler = AmplitudeHandler(install_path,data_path)
        handler.insertEvent(event)
        handler.insertAmps(xmlfile)
        nassociated = handler.associateAll(write_data=False)
        assert nassociated == 1
    except:
        assert 1==2
    finally:
        if os.path.isfile(dbfile):
            os.remove(dbfile)
    
if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_amps()
