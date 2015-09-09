#!/usr/bin/env python

from neicio.gmt import GMTGrid
from openquake.hazardlib.gsim.base import SitesContext

class SitesException(Exception):
    """
    Class to represent errors in the Sites class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def getGeodict(bounds,xdim,ydim):
    lons = np.arange(bounds[2],bounds[3]+ydim,ydim)
    lats = np.arange(bounds[0],bounds[1]+xdim,xdim)
    geodict = {'xmin':bounds[0],
               'xmax':bounds[1],
               'ymin':bounds[2],
               'ymax':bounds[3],
               'xdim':xdim,
               'ydim':ydim,
               'nrows':len(lons),
               'ncols':len(lats)}
    return (geodict,lons,lats)

class Sites(object):
    def __init__(self,hypocenter,bounds,xdim,ydim,vs30File=None,defaultVs30=760):
        self.GeoDict,lons,lats = getGeodict(bounds,xdim,ydim)
        if vs30File is not None:
            vsgrid = GMTGrid(vs30File)
            bigbounds = (bounds[0]-xdim*4,bounds[1]+xdim*4,bounds[2]-ydim*4,bounds[3]+ydim*4)
            vsgrid.load(bigbounds)
            vsgrid.interpolateToGrid(self.GeoDict)
            self.Vs30 = vsgrid.griddata
        else:
            self.Vs30 = np.ones((geodict['nrows'],geodict['ncols']))*defaultVs30
        self._calculateZ1P0()
        self._calculateZ2P5()
        self.SitesContext = SitesContext()
        self.SitesContext.vs30 = self.Vs30
        self.SitesContext.z1pt0 = self.Z1Pt0
        self.SitesContext.z2pt5 = self.Z2Pt5
        self.SitesContext.backarc = False #zoneconfig might have this info
        self.SitesContext.vs30measured = False #no idea where this might ever come from
        self.SitesContext.lons = lons
        self.SitesContext.lats = lats
        
    def getSitesContext(self):
        return self.SitesContext

    def _calculateZ2P5(self):
        self.Z2Pt5 = self.Vs30*2.5 #fill in with real equation here

    def _calculateZ1P5(self):
        self.Z2Pt5 = self.Vs30*1.0 #fill in with real equation here
            
        
