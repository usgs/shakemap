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
from shapely import geometry
from timeutils import ShakeDateTime
from ecef import latlon2ecef,ecef2latlon

REQUIRED_KEYS = ['lat','lon','depth','time','mag']
OPTIONAL_KEYS = ['id','timezone','locstring','type','created']
DEG2RAD = 180./np.pi

def readEventFile(eventxml):
    event = {}
    root = minidom.parse(eventxml)
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

def readFaultFile(faultfile):
    """
    Read fault file format as defined in ShakeMap Software Guide.  
    Input: 
    Fault file in GMT psxy format, where

     * Fault vertices are space separated lat,lon,depth triplets on a single line.
     * Fault segments are separated by lines containing ">"
     * Fault segments must be closed.
     * Fault segments must be all clockwise or all counter-clockwise.

    Raises Exception when above conditions are not met.
    """
    x = []
    y = []
    z = []
    if isinstance(faultfile,str) or isinstance(faultfile,unicode):
        faultfile = open(faultfile,'rt')
        faultlines = faultfile.readlines()
    else:
        faultlines = faultfile.readlines()
    reference = ''
    for line in faultlines:
        sline = line.strip()
        if sline.startswith('#'):
            reference += sline
            continue
        if sline.startswith('>'):
            if len(x): #start of new line segment
                x.append(np.nan)
                y.append(np.nan)
                z.append(np.nan)
                continue
            else: #start of file
                continue
        parts = sline.split()
        if len(parts) < 3:
            raise Exception('Finite fault file %s has no depth values.' % faultfile)
        y.append(float(parts[0]))
        x.append(float(parts[1]))
        z.append(float(parts[2]))
    faultfile.close()
    if np.isnan(x[-1]):
        x = x[0:-1]
        y = y[0:-1]
        z = z[0:-1]
        
    return (x,y,z,reference)
    

class Source(object):
    def __init__(self,event,lon=None,lat=None,depth=None,reference=None):
        missing = []
        for req in REQUIRED_KEYS:
            if req not in event.keys():
                missing.append(req)
        if len(missing):
           raise NameError('Input event dictionary is missing the following required keys: "%s"' % (','.join(missing)))
        allNone = lon is None and lat is None and depth is None
        noNone = lon is not None and lat is not None and depth is not None
        if not allNone and not noNone:
            raise ValueError('lon,lat,depth input parameters must be all valid sequence objects or all None.')
        self.EventDict = event.copy()
        self.reference = reference
        if noNone:
            self.setCoordinates(lon,lat,depth)
        else:
            self.setCoordinates([event['lon']],[event['lat']],[event['depth']])

    @classmethod
    def readFromFile(cls,eventxmlfile,faultfile=None):
        event = readEventFile(eventxmlfile)
        if faultfile is not None:
            lon,lat,depth,reference = readFaultFile(faultfile)
        else:
            lon = None
            lat = None
            depth = None
        return cls(event,lon=lon,lat=lat,depth=depth,reference=reference)
            
    def getEventDict(self):
        return self.EventDict

    def setEventParam(self,key,value):
        self.EventDict[key] = value

    def getEventParam(self,key):
        if key not in self.EventDict.keys():
            raise NameError('Key "%s" not found in event dictionary' % key)
        return self.EventDict[key]

    def setFaultReference(self,reference):
        self.reference = reference

    def getFaultReference(self):
        return self.reference
    
    def setCoordinates(self,lon,lat,depth):
        lon = np.array(lon)
        lat = np.array(lat)
        depth = np.array(depth)
        inan = np.isnan(lon)
        numnans = len(lon[inan])
        numsegments = numnans + 1
        #requirements:
        # 1) Coordinate arrays must be same length
        # 2) Polygons must be quadrilaterals
        # 3) Quads must be closed
        # 4) Quads must be planar
        if len(lon) != len(lat) != len(depth):
            raise IndexError('Length of input lon,lat,depth arrays must be equal')
        
        istart = 0
        endpoints = list(np.where(np.isnan(lon))[0])
        endpoints.append(len(lon))
        surfaces = []
        for iend in endpoints:
            lonseg = lon[istart:iend]
            latseg = lat[istart:iend]
            depthseg = depth[istart:iend]
            #each segment can have many contiguous quadrilaterals defined in it
            #separations (nans) between segments mean that segments are not contiguous.
            npoints = len(lonseg)
            nquads = ((npoints - 5)/2) + 1
            ioff = 0
            for i in range(0,nquads):
                endidx = ioff-1 #we have the closing vertex that we're not interested in here
                topLeft = point.Point(lonseg[ioff],latseg[ioff],depthseg[ioff])
                topRight = point.Point(lonseg[ioff+1],latseg[ioff+1],depthseg[ioff])
                bottomRight = point.Point(lonseg[endidx-2],latseg[endidx-2],depthseg[endidx-2])
                bottomLeft = point.Point(lonseg[endidx-1],latseg[endidx-1],depthseg[endidx-1])
                rect = self.getPlanarRect(topLeft,topRight,bottomRight,bottomLeft)
                surface = planar.PlanarSurface.from_corner_points(0.01,rect[0],rect[1],rect[2],rect[3])
                surfaces.append(surface)
                ioff += 1
            istart = iend+1

        self.MultiSurface = multi.MultiSurface(surfaces)
        self.lon = lon
        self.lat = lat
        self.depth = depth

    def getPlanarRect(self,topLeft,topRight,bottomRight,bottomLeft):
        lon1 = topLeft.longitude
        lat1 = topLeft.latitude
        dep1 = topLeft.depth
        lon2 = topRight.longitude
        lat2 = topRight.latitude
        dep2 = topRight.depth
        lon3 = bottomRight.longitude
        lat3 = bottomRight.latitude
        dep3 = bottomRight.depth
        lon4 = bottomLeft.longitude
        lat4 = bottomLeft.latitude
        dep4 = bottomLeft.depth

        lats = np.array([lat1,lat2,lat3,lat4])
        lons = np.array([lon1,lon2,lon3,lon4])
        depths = np.array([dep1,dep2,dep3,dep4])
        
        points = np.vstack((lons,lats,depths)).T
        ctrPoint,orthogonal = utils.plane_fit(points)
        toplen = distance(lons[0],lats[0],depths[0],lons[1],lats[1],depths[1])
        bottomlen = distance(lons[2],lats[2],depths[2],lons[3],lats[3],depths[3])
        rightheight = distance(lons[1],lats[1],depths[1],lons[2],lats[2],depths[2])
        leftheight = distance(lons[3],lats[3],depths[3],lons[0],lats[0],depths[0])
        width = max(toplen,bottomlen)
        height = max(leftheight,rightheight)
        mlat = np.mean(lats)*DEG2RAD
        width = np.cos(mlat) * (width / 111.191)
        height = height / 111.191
        dx = width/2
        dy = height/2
        i,j,k = orthogonal
        a,b,c = ctrPoint
        z1 = (-i*dx - j*dy)/k + c
        z2 = (-i*dx - j*dy)/k + c
        z3 = (i*dx - j*dy)/k + c
        z4 = (i*dx - j*dy)/k + c

        xplane = np.array([a-dx,a+dx,a+dx,a-dx])
        yplane = np.array([b+dy,b+dy,b-dy,b-dy])
        zplane = np.array([z1,z2,z3,z4])

        topLeft = point.Point(xplane[0],yplane[0],-zplane[0])
        topRight = point.Point(xplane[1],yplane[1],-zplane[1])
        bottomRight = point.Point(xplane[2],yplane[2],-zplane[2])
        bottomLeft = point.Point(xplane[3],yplane[3],-zplane[3])
        return (topLeft,topRight,bottomRight,bottomLeft)

    def getSurface(self):
        return self.MultiSurface
    
    def calcDistance(self,mesh,method='rjb'):
        if method == 'rjb':
            if len(self.lon) == 1:
                p = point.Point(self.lon[0],self.lat[0],0.0)
                return p.distance_to_mesh(mesh)
            else:
                return self.MultiSurface.get_joyner_boore_distance(mesh)
        elif method == 'rrup':
            if len(self.lon) == 1:
                p = point.Point(self.lon[0],self.lat[0],self.depths[0])
                return p.distance_to_mesh(mesh)
            else:
                return self.MultiSurface.get_min_distance(mesh)
        elif method == 'rx':
            if len(self.lon) == 1:
                return np.zeros(mesh.shape)
            else:
                return self.MultiSurface.get_rx_distance(mesh)
        else:
            raise NotImplementedError('No distance measure exists for method type "%s"' % method)
        


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
            
        
