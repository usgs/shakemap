#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import io

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..','..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed shakemap stuff

#local imports
from shakemap.utils.exception import ShakeMapException
from shakemap.grind.source import Source

#third party
from openquake.hazardlib.gsim import base,abrahamson_2014
from openquake.hazardlib.geo.mesh import Mesh
import numpy as np

def test():
    fault_text = """30.979788       103.454422      1
31.691615       104.419160      1
31.723569       104.374760      1
32.532213       105.220821      1
32.641450       105.135050      20
31.846790       104.246202      20
31.942158       104.205286      20
31.290105       103.284388      20
30.979788       103.454422      1"""
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01" timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681" type="" />
    """
    print('Testing creation of source object...')
    source_text = """mech=RS"""
    ffile = io.StringIO(fault_text)
    efile = io.StringIO(event_text)
    sfile = io.StringIO(source_text)
    source = Source.readFromFile(efile,faultfile=ffile,sourcefile=sfile)
    print('Passed creation of source object.')

    print('Testing creation of RuptureContext object...')
    gmpe = abrahamson_2014.AbrahamsonEtAl2014()
    rupture = source.getRuptureContext([gmpe])
    testdict = {'mag':7.9,
                'strike': -133.083550974,
                'dip': 49.8524115024,
                'rake': 45.0,
                'ztor':0.999999999995,
                'hypo_lon':103.3639,
                'hypo_lat':30.9858,
                'hypo_depth':19.0,
                'width':27.8623813381}
    for key in testdict.keys():
        value = eval('rupture.%s' % key)
        np.testing.assert_almost_equal(testdict[key],value)
    print('Passed creation of RuptureContext object...')
    
    print('Test setting mechanism and rake/dip...')
    mech = 'RS'
    exp_dip = 40
    exp_rake = 90
    source.setMechanism(mech)
    assert source.getEventParam('dip') == exp_dip
    assert source.getEventParam('rake') == exp_rake
    source.setMechanism('ALL',dip=45,rake=315)
    assert source.getEventParam('rake') == -45
    #this should raise an exception
    try:
        source.setMechanism('ALL',dip=110)
    except ShakeMapException as sme:
        print('Exception raised appropriately for dip greater than 90.')
    #this should raise an exception
    try:
        source.setMechanism('ALL',rake=370)
    except ShakeMapException as sme:
        print('Exception raised appropriately for rake greater than 360.')
    print('Test setting mechanism and rake/dip...')

def _test_northridge():
    fault_text = """
    # Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """
    event_text = """<?xml version="1.0" encoding="US-ASCII" standalone="yes"?>
<earthquake id="blah" lat="34.213" lon="-118.537" mag="7.9" year="1994" month="01" day="17" hour="12" minute="30" second="55" timezone="GMT" depth="18.4" locstring="NORTHRIDGE" created="1211173621" otime="1210573681" type="" />
    """
    source_text = """mech=RS"""
    ffile = io.StringIO(fault_text)
    efile = io.StringIO(event_text)
    sfile = io.StringIO(source_text)
    source = Source.readFromFile(efile,faultfile=ffile,sourcefile=sfile)
    gmpe = abrahamson_2014.AbrahamsonEtAl2014()
    rupture = source.getRuptureContext(gmpe)
    mapwidth = 2.0
    latmin = rupture.hypo_lat - mapwidth
    latmax = rupture.hypo_lat + mapwidth
    lonmin = rupture.hypo_lon - mapwidth
    lonmax = rupture.hypo_lon + mapwidth
    dim = 0.02
    lats = np.arange(latmin,latmax,dim)
    lons = np.arange(lonmin,lonmax,dim)
    lon,lat = np.meshgrid(lons,lats)
    dep = np.zeros_like(lon)
    mesh = Mesh(lon,lat,dep)
    distances = source.getDistanceContext(gmpe,mesh)
    rupture = source.getRuptureContext(gmpe)
    for key in rupture._slots_:
        try:
            value = eval('rupture.%s' % key)
        except:
            print('No value set for %s' % key)
            continue
        print('%s = %s' % (key,str(value)))    

    cbuf = io.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
        
if __name__ == '__main__':
    test()
            
        
