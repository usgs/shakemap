#!/usr/bin/env python

#stdlib imports
from xml.dom import minidom
import io
import sys

#third party
import pytz
import numpy as np
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo import point
from openquake.hazardlib.geo import Mesh
from openquake.hazardlib.geo.surface import multi
from openquake.hazardlib.geo.surface import planar
from openquake.hazardlib.geo import utils
from openquake.hazardlib.gsim import base,abrahamson_2014
from openquake.hazardlib.const import TRT
from shapely import geometry
from .timeutils import ShakeDateTime
from .ecef import latlon2ecef,ecef2latlon
from openquake.hazardlib.gsim.base import GMPE
from .fault import Fault
from .distance import getDistance

#local imports
from shakemap.utils.exception import ShakeMapException

REQUIRED_KEYS = ['lat','lon','depth','time','mag']
OPTIONAL_KEYS = ['id','timezone','locstring','type','created']
DEG2RAD = 180./np.pi
RAKEDICT = {'SS':0.0,'NM':-90.0,'RS':90.0,'ALL':None}
DEFAULT_STRIKE = 0.0
DEFAULT_DIP = 90.0
DEFAULT_WIDTH = 0.0
DEFAULT_ZTOR = 0.0


def readEventFile(eventxml):
    """
    Read event.xml file from disk, returning a dictionary of attributes.
    Input XML format looks like this:
    ******
    <earthquake id="2008ryan" lat="30.9858" lon="103.3639" mag="7.9" year="2008" month="05" day="12" hour="06" minute="28" second="01" timezone="GMT" depth="19.0" locstring="EASTERN SICHUAN, CHINA" created="1211173621" otime="1210573681" type="" />
    ******
    :param eventxml:
        Path to event XML file OR file-like object
    :returns:
       Dictionary with keys as indicated above in earthquake element attributes.
    """
    event = {}
    if isinstance(eventxml,str):
        root = minidom.parse(eventxml)
    else:
        data = eventxml.read()
        root = minidom.parseString(data)
        
    eq = root.getElementsByTagName('earthquake')[0]
    event['lat'] = float(eq.getAttribute('lat'))
    event['lon'] = float(eq.getAttribute('lon'))
    event['depth'] = float(eq.getAttribute('depth'))
    event['mag'] = float(eq.getAttribute('mag'))
    year = int(eq.getAttribute('year'))
    month = int(eq.getAttribute('month'))
    day = int(eq.getAttribute('day'))
    hour = int(eq.getAttribute('hour'))
    minute = int(eq.getAttribute('minute'))
    second = float(eq.getAttribute('second'))
    msec = int((second - int(second))*1e6)
    second = int(second)

    #TODO: Handle:
    # 1) timezones other than UTC (use pytz)
    # 2) times < year 1900 (strftime fails) - subclass datetime, use Obspy utcdatetime
    # 3) daylight savings time (use pytz)
    
    event['time'] = ShakeDateTime(year,month,day,hour,minute,second,msec)
    for key in ['id','locstring','type','timezone']:
        if eq.hasAttribute(key):
            event[key] = eq.getAttribute(key)
    if 'created' in event:
        event['created'] = ShakeDateTime.utcfromtimestamp(int(eq.getAttribute('created')))

    root.unlink()
    return event

def readSource(sourcefile):
    """
    Read source.txt file, which has lines like key=value.
    :param sourcefile:
        Path to source.txt file OR file-like object
    :returns:
        Dictionary containing key/value pairs from file.
    """
    isFile = False
    if isinstance(sourcefile,str):
        f = open(sourcefile,'rt')
    else:
        isFile = True
        f = sourcefile
    params = {}
    for line in f.readlines():
        if line.startswith('#'):
            continue
        parts = line.split('=')
        key = parts[0].strip()
        value = parts[1].strip()
        #see if this is some sort of number
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass
        params[key] = value
    if not isFile:
        f.close()
    return params

class Source(object):
    """
    Encapsulate everything we need to know about an earthquake source.
    """
    def __init__(self,event,fault=None,sourcedict=None):
        """
        Construct a Source object.
        :param event:
            dictionary of values (see readEventFile())
        :param fault:
            a Fault object
        :param sourcedict:
            Dictionary containing values from source.txt file (see readSource())
        :returns:
            Source object.
        """
        missing = []
        for req in REQUIRED_KEYS:
            if req not in list(event.keys()):
                missing.append(req)
        if len(missing):
           raise NameError('Input event dictionary is missing the following required keys: "%s"' % (','.join(missing)))
        self.Fault = fault
        self.EventDict = event.copy()
        if isinstance(sourcedict,dict):
            for key,value in sourcedict.items():
                self.setEventParam(key,value)
        

    @classmethod
    def readFromFile(cls,eventxmlfile,faultfile=None,sourcefile=None):
        """
        Class method to create a Source object by specifying an event.xml file, a fault file, and a source.txt file.
        :param eventxmlfile:
            Event xml file (see readEventFile())
        :param faultfile:
            Fault text file (see fault.py)
        :param sourcefile:
            source.txt file (see readSource())
        :returns:
            Source object
        """
        event = readEventFile(eventxmlfile)
        if faultfile is not None:
            fault = Fault.readFaultFile(faultfile)
        else:
            fault=None
        params = None
        if sourcefile is not None:
            params = readSource(sourcefile)
        return cls(event,fault=fault,sourcedict=params)

    def getDistanceContext(self,gmpe,mesh):
        """
        Create a DistancesContext object
        :param gmpe: 
            Concrete subclass of GMPE (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
        :param mesh:
            Mesh object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        :returns:
            DistancesContext object with distance grids required by input gmpe.
        :raises TypeError:
            if gmpe is not a subclass of GMPE
        """
        if not isinstance(gmpe,GMPE):
            raise TypeError('getDistanceContext() cannot work with objects of type "%s"' % type(gmpe))
        context = base.DistancesContext()
        for method in gmpe.REQUIRES_DISTANCES:
            print('Calculating method %s' % method)
            dist = self.calcDistance(mesh,method)
            #can't figure out how to dynamically assign dist to a property on context
            #this doesn't work:
            #context.__dict__[method] = dist
            if method == 'rx':
                context.rx = dist.reshape(mesh.lons.shape)
            if method == 'rrup':
                context.rrup = dist.reshape(mesh.lons.shape)
            if method == 'rjb':
                context.rjb = dist.reshape(mesh.lons.shape)
            if method == 'ry0':
                context.ry0 = dist.reshape(mesh.lons.shape)
            if method == 'rcdpp':
                context.rcdpp = dist.reshape(mesh.lons.shape)
            if method == 'repi':
                context.repi = dist.reshape(mesh.lons.shape)
            if method == 'rhypo':
                context.rhypo = dist.reshape(mesh.lons.shape)

        return context

    def getRuptureContext(self,gmpe):
        """
        Return a GEM RuptureContext https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py
        :returns:
            RuptureContext object with all known parameters filled in.
        If Source does not contain a Fault, strike, dip, ztor, and width will be filled with default values.
        Rake may not be known, or may be estimated from a focal mechanism.
        """
        #rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'
        reqset = gmpe.REQUIRES_RUPTURE_PARAMETERS
        if 'hypo_loc' in reqset:
            raise ShakeMapException('Rupture parameter "hypo_loc" is not supported!')
        rup = base.RuptureContext()
        rup.mag = self.getEventParam('mag')
        if self.Fault is not None:
            rup.strike = self.Fault.getStrike()
            rup.dip = self.Fault.getDip()
            rup.ztor = self.Fault.getTopOfRupture()
            rup.width = self.Fault.getWidth()
        else:
            rup.strike = DEFAULT_STRIKE
            rup.dip = DEFAULT_DIP
            rup.ztor = DEFAULT_ZTOR
            rup.width = DEFAULT_WIDTH

        if 'rake' in self.EventDict:
            rup.rake = self.getEventParam('rake')
        elif 'mech' in self.EventDict:
            mech = self.EventDict['mech']
            rup.rake = RAKEDICT[mech]
        else:
            rup.rake = RAKEDICT['ALL']
        
        rup.hypo_lat = self.getEventParam('lat')
        rup.hypo_lon = self.getEventParam('lon')
        rup.hypo_depth = self.getEventParam('depth')
        return rup

    def setTectonicRegion(self,region):
        """
        Set tectonic region.
        :param region:
            TRT object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/const.py)
        """
        if not isinstance(region,TRT):
            raise ValueError('Input tectonic region must be of type openquake.hazardlib.const.TRT')
        self.TectonicRegion = region

    def getTectonicRegion(self):
        """
        Get tectonic region.
        :returns:
            TRT object (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/const.py)
        """
        return self.TectonicRegion
        
    def getEventDict(self):
        """
        Return the event dictionary.
        :returns:
           Copy of dictionary of event keys/values
        """
        return self.EventDict.copy()

    def setEventParam(self,key,value):
        """
        Set a parameter in the internal event dictionary
        :param key:
            string key
        :param value:
            value (any object type)
        """
        self.EventDict[key] = value

    def getEventParam(self,key):
        """
        Get a parameter from the internal event dictionary
        :param key:
            string key
        :returns:
            value (any object type)
        """
        if key not in list(self.EventDict.keys()):
            raise NameError('Key "%s" not found in event dictionary' % key)
        return self.EventDict[key]

    def setFaultReference(self,reference):
        """
        Set the citeable reference for the fault
        :param reference:
            string citeable reference
        """
        self.Fault.setReference(reference)

    def getFaultReference(self):
        """
        Get the citeable reference for the fault
        :returns:
            string citeable reference
        """
        return self.Fault.getReference()

    def getQuadrilaterals(self):
        """
        Get the list of quadrilaterals defining the fault.
        :returns:
            List of quadrilateral tuples (4 Point objects each)
        """
        return self.Fault.getQuadrilaterals()
    
    def calcDistance(self,mesh,method='rjb'):
        """
        Calculate distance from the source (fault or point source)
        :param mesh:
          Instance of an openquake.hazardlib.geo.Mesh object.
        :param method:
          One of: 'rjb','rx','rrup','ry0','rcdbp','epi','hypo'. NB - If no fault is provided, 'rjb' and 'rx' => 'repi',
          anything else => 'rhypo'.
        """
        if self.Fault is not None:
            quadlist = self.Fault.getQuadrilaterals()
        else:
            quadlist = None
        lat = self.getEventParam('lat')
        lon = self.getEventParam('lon')
        depth = self.getEventParam('depth')
        oqpoint = point.Point(lon,lat,depth)
        #if user has specified a distance measure dependent on a fault but we don't have a fault
        #override their choice with the appropriate point source distance.
        if quadlist is None and (method != 'repi' or method != 'rhypo'):
            if method in ['rjb','rx']:
                method = 'repi'
            else:
                method = 'rhypo'
            
        return getDistance(method,mesh,quadlist=quadlist,mypoint=oqpoint)

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
    for key in rupture.__slots__:
        try:
            value = eval('rupture.%s' % key)
        except:
            print('No value set for %s' % key)
            continue
        print('%s = %s' % (key,str(value)))

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
    for key in rupture.__slots__:
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
            
        
