#!/usr/bin/env python

#stdlib imports
from xml.dom import minidom
import sys

#third party
import pytz
import numpy as np
from openquake.hazardlib.geo import geodetic
from openquake.hazardlib.geo import point
from openquake.hazardlib.geo import Mesh
from openquake.hazardlib.geo.geodetic import distance
from openquake.hazardlib.geo.surface import multi
from openquake.hazardlib.geo.surface import planar
from openquake.hazardlib.geo import utils
from openquake.hazardlib.gsim import base
from openquake.hazardlib.const import TRT
from shapely import geometry
from timeutils import ShakeDateTime
from ecef import latlon2ecef,ecef2latlon
from openquake.hazardlib.gsim.base import GMPE

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
    if isinstance(eventxml,str) or isinstance(eventxml,unicode):
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
        if event.has_key(key):
            event[key] = eq.getAttribute(key)
    if event.has_key('created'):
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
    if isinstance(sourcefile,str) or isinstance(sourcefile,unicode):
        f = open(sourcefile,'rt')
    else:
        isFile = True
        f = sourcefile
    params = {}
    for line in f.readlines:
        parts = line.split('=')
        key = parts[0].strip()
        value = params[1].strip()
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
            if req not in event.keys():
                missing.append(req)
        if len(missing):
           raise NameError('Input event dictionary is missing the following required keys: "%s"' % (','.join(missing)))
        self.Fault = fault
        self.EventDict = event.copy()
        self.TectonicRegion = tectonicRegion
        if isinstance(sourcedict,dict):
            for key,value in sourcedict.iteritems():
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
        for method in gmpe.REQUIRES_DISTANCES:
            context = base.DistancesContext
            dist = self.calcDistance(mesh,method)
            eval('context.%s = dist' % method)

        return context

    def getRuptureContext(self):
        """
        Return a GEM RuptureContext https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py
        :returns:
            RuptureContext object with all known parameters filled in.
        If Source does not contain a Fault, strike, dip, ztor, and width will be filled with default values.
        Rake may not be known, or may be estimated from a focal mechanism.
        """
        #gmpe is a subclass of hazardlib GMPE
        params = gmpe.REQUIRES_RUPTURE_PARAMETERS
        #rupturecontext constructor inputs:
        # 'mag', 'strike', 'dip', 'rake', 'ztor', 'hypo_lon', 'hypo_lat',
        # 'hypo_depth', 'width', 'hypo_loc'
        mag = self.getEventParam('mag')
        if self.Fault is not None:
            strike = self.Fault.getStrike()
            dip = self.Fault.getDip()
            ztor = self.Fault.getTopOfRupture()
            width = self.Fault.getWidth()
        else:
            strike = DEFAULT_STRIKE
            dip = DEFAULT_DIP
            ztor = DEFAULT_ZTOR
            width = DEFAULT_WIDTH

        if self.EventDict.has_key('rake'):
            rake = self.getEventParam('rake')
        elif self.EventDict.has_key('mech'):
            mech = self.EventDict['mech']
            rake = RAKEDICT[mech]
        else:
            rake = RAKEDICT['ALL']
        
        region = self.TectonicRegion
        lat = self.getEventParam('lat')
        lon = self.getEventParam('lon')
        depth = self.getEventParam('depth')
        locstr = self.getEventParam('locstr')
        rup = base.RuptureContext(mag,strike,dip,rake,ztor,lon,lat,depth,width,locstr)

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
        if key not in self.EventDict.keys():
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
        
        """
        if self.Fault is not None:
            quadlist = self.Fault.getQuadrilaterals()
        else:
            quadlist = None
        lat = self.getEventParam('lat')
        lon = self.getEventParam('lon')
        depth = self.getEventParam('depth')
        point = point.Point(lat,lon,depth)
        #if user has specified a distance measure dependent on a fault but we don't have a fault
        #override their choice with the appropriate point source distance.
        if quadlist is None and (method != 'epi' or method != 'hypo'):
            if method in ['rjb','rx']:
                method = 'epi'
            else:
                method = 'hypo'
            
        return distance.getDistance(method,mesh,quadlist=quadlist,point=point)

def test(eventxml=None,faultfile=None):
    event = {'lat':28.230,'lon':84.731,'depth':8.2,
             'mag':7.8,'time':ShakeDateTime.utcnow()}
    xmin = event['lon'] - (2.5/2.0)
    xmax = event['lon'] + (2.5/2.0)
    ymin = event['lat'] - (2.5/2.0)
    ymax = event['lat'] + (2.5/2.0)
    xdim = 0.0083
    ydim = 0.0083

    fault = np.array([[28.427,84.614,20],
                      [27.986,86.179,20],
                      [27.469,85.880,13],
                      [27.956,84.461,13],
                      [28.427,84.614,20]])
    lons = fault[:,1]
    lats = fault[:,0]
    depths = fault[:,2]

    reference = 'Schmidt, J.J.J 1999 Journal of GeoTectonical Science'
    source = Source(event,lon=lons,lat=lats,depth=depths,reference=reference)
    x,y = np.meshgrid(np.arange(xmin,xmax,xdim),np.arange(ymin,ymax,ydim))
    shakemap = Mesh(x,y)
    print shakemap.shape
    dist = source.calcDistance(shakemap,method='rjb')
    #test reading of files
    lon = None
    lat = None
    depth = None
    
    if eventxml is not None:
        source = Source.readFromFile(eventxml,faultfile=faultfile)
        dist = source.calcDistance(shakemap,method='rjb')
        
if __name__ == '__main__':
    if len(sys.argv) > 1:
        eventfile = sys.argv[1]
    if len(sys.argv) > 2:
        faulttxt = sys.argv[2]
    test(eventxml=eventfile,faultfile=faulttxt)
            
        
