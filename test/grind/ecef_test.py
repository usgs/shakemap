#!/usr/bin/env python

#stdlib imports
import os.path
import sys

#third party
import numpy as np

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.ecef import latlon2ecef,ecef2latlon

def test():
    print('Testing ECEF conversion code...')
    lat = 32.1
    lon = 118.5
    dep = 10.0

    x,y,z = latlon2ecef(lat,lon,dep)

    TESTX = -2576515.419
    TESTY = 4745351.087
    TESTZ = 3364516.124
    
    np.testing.assert_almost_equal(x,TESTX,decimal=2)
    np.testing.assert_almost_equal(y,TESTY,decimal=2)
    np.testing.assert_almost_equal(z,TESTZ,decimal=2)

    lat2,lon2,dep2 = ecef2latlon(x,y,z)
    np.testing.assert_almost_equal(lat2,lat,decimal=2)
    np.testing.assert_almost_equal(lon2,lon,decimal=2)
    np.testing.assert_almost_equal(dep2,dep,decimal=2)
    print('Passed tests of ECEF conversion code.')
    
if __name__ == '__main__':
    test()
    
