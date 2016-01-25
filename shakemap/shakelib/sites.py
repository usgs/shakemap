#!/usr/bin/env python

#stdlib imports
import sys

#third party imports
from neicio.gmt import GMTGrid
from openquake.hazardlib.gsim.base import SitesContext
import numpy as np


class SitesException(Exception):
    """
    Class to represent errors in the Sites class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def getGeodict(bounds,xdim,ydim):
    """Return a geodict dictionary as defined here: https://github.com/usgs/neicio/blob/master/neicio/grid.py
    :param bounds:
        Tuple of floats containing (lonmin,lonmax,latmin,latmax)
    :param xdim:
        Float width of desired cells in decimal degrees.
    :param ydim:
        Float height of desired cells in decimal degrees.
    :returns:
       Dictionary containing fields as referenced above:
        - xmin
        - xmax
        - ymin
        - ymax
        - xdim
        - ydim
        - nrows
        - ncols
    """
    lons = np.arange(bounds[0],bounds[1]+xdim,xdim)
    lats = np.arange(bounds[2],bounds[3]+ydim,ydim)
    geodict = {'xmin':bounds[0],
               'xmax':bounds[1],
               'ymin':bounds[2],
               'ymax':bounds[3],
               'xdim':xdim,
               'ydim':ydim,
               'nrows':len(lats),
               'ncols':len(lons)}
    return (geodict,lons,lats)

class Sites(object):
    """An object to encapsulate information used to generate a GEM SitesContext.
    (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
    """
    def __init__(self,bounds,xdim,ydim,vs30File=None,defaultVs30=686,vs30measured=False,backarc=False):
        """
        Construct a Sites object.
        :param bounds:
            Tuple of floats containing (lonmin,lonmax,latmin,latmax)
        :param xdim:
            Float width of desired cells in decimal degrees.
        :param ydim:
            Float height of desired cells in decimal degrees.
        :param vs30File:
            Path to NetCDF file containing Vs30 values (m/s).
        :param defaultVs30:
            Default Vs30 value.
        :param vs30measured:
            Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
            Boolean indicating whether event is on the backarc as defined here: 
            http://earthquake.usgs.gov/learn/glossary/?term=backarc
        """
        self.GeoDict,lons,lats = getGeodict(bounds,xdim,ydim)
        if vs30File is not None:
            bigbounds = (bounds[0]-xdim*4,bounds[1]+xdim*4,bounds[2]-ydim*4,bounds[3]+ydim*4)
            vsgrid = GMTGrid(grdfile=vs30File,bounds=bigbounds)
            vsgrid.interpolateToGrid(self.GeoDict)
            self.Vs30 = vsgrid.griddata
        else:
            self.Vs30 = np.ones((self.GeoDict['nrows'],self.GeoDict['ncols']))*defaultVs30
        self._calculateZ1P0()
        self._calculateZ2P5()
        self.SitesContext = SitesContext()
        self.SitesContext.vs30 = self.Vs30
        self.SitesContext.z1pt0 = self.Z1Pt0
        self.SitesContext.z2pt5 = self.Z2Pt5
        self.SitesContext.backarc = backarc #zoneconfig might have this info
        self.SitesContext.vs30measured = vs30measured #no idea where this might ever come from
        self.SitesContext.lons = lons
        self.SitesContext.lats = lats
        
    def getSitesContext(self):
        """Return a SitesContext object appropriate for GMPE use.
        :returns:
           SitesContext object.
        """
        return self.SitesContext

    def _calculateZ1P0(self):
        c1 = 6.745
        c2 = 1.35
        c3 = 5.394
        c4 = 4.48
        self.Z1Pt0 = np.zeros_like(self.Vs30)
        self.Z1Pt0[self.Vs30 < 180] = np.exp(c1)
        idx = (self.Vs30 >= 180) & (self.Vs30 <= 500)
        self.Z1Pt0[idx] = np.exp(c1 - c2*np.log(self.Vs30[idx]/180.0))
        idx = self.Vs30 > 500
        self.Z1Pt0[idx] = np.exp(c3 -  c4 * np.log(self.Vs30[idx]/500.0))

    def _calculateZ2P5(self):
        c1 = 519
        c2 = 3.595
        self.Z2Pt5 = c1 + self.Z1Pt0*c2
            
        
def _test(vs30file=None):
    hypo = {'lat':34.1,'lon':-118.2,'depth':23.5}
    bounds = (hypo['lon']-2.5,hypo['lon']+2.5,hypo['lat']-2.5,hypo['lat']+2.5)
    xdim = 0.0083
    ydim = 0.0083
    mysite = Sites(bounds,xdim,ydim,vs30File=vs30file,defaultVs30=760)
    sc = mysite.getSitesContext()
    pass

if __name__ == '__main__':
    if len(sys.argv) > 1:
        vs30file = sys.argv[1]
    _test(vs30file=vs30file)
