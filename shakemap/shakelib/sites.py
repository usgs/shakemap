#!/usr/bin/env python

#stdlib imports
import sys

#third party imports
from mapio.gmt import GMTGrid
from mapio.geodict import GeoDict
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

def calculateZ1P0(vs30)
    c1 = 6.745
    c2 = 1.35
    c3 = 5.394
    c4 = 4.48
    Z1Pt0 = np.zeros_like(vs30)
    Z1Pt0[vs30 < 180] = np.exp(c1)
    idx = (vs30 >= 180) & (vs30 <= 500)
    Z1Pt0[idx] = np.exp(c1 - c2*np.log(vs30[idx]/180.0))
    idx = vs30 > 500
    Z1Pt0[idx] = np.exp(c3 -  c4 * np.log(vs30[idx]/500.0))
    return Z1Pt0

def calculateZ2P5(z1pt0):
    c1 = 519
    c2 = 3.595
    Z2Pt5 = c1 + z1pt0*c2
    return Z2Pt5

class Sites(object):
    """An object to encapsulate information used to generate a GEM SitesContext.
    (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/gsim/base.py)
    """
    def __init__(self,vs30grid,vs30measured=False,backarc=false):
        """
        Construct a Sites object.
        :param vs30grid:
            MapIO Grid2D object containing Vs30 values.
        :param vs30measured:
            Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
            Boolean indicating whether event is on the backarc as defined here: 
            http://earthquake.usgs.gov/learn/glossary/?term=backarc
        """
        self.Vs30 = vs30grid
        self.GeoDict = vs30grid.getGeoDict().copy()
        lons = np.arange(self.GeoDict.xmin,self.GeoDict.xmax,self.GeoDict.dx)
        lats = np.arange(self.GeoDict.ymin,self.GeoDict.ymax,self.GeoDict.dy)
        self.Z1Pt0 = calculateZ1P0(self.Vs30)
        self.Z2Pt5 = calculateZ2P5(self.Z1Pt0)
        self.SitesContext = SitesContext()
        self.SitesContext.vs30 = self.Vs30.getData().copy()
        self.SitesContext.z1pt0 = self.Z1Pt0
        self.SitesContext.z2pt5 = self.Z2Pt5
        self.SitesContext.backarc = backarc #zoneconfig might have this info
        self.SitesContext.vs30measured = vs30measured #no idea where this might ever come from
        self.SitesContext.lons = lons
        self.SitesContext.lats = lats

    @classmethod
    def createFromCenter(cls,cx,cy,dx,dy,xspan,yspan,defaultVs30=686.0,vs30File=None,vs30measured=False,backarc=False):
        """Create a Sites object by defining a center point, resolution, extent, and Vs30 values.

        :param cx:
          X coordinate of desired center point.
        :param cy:
          X coordinate of desired center point.
        :param dx:
          Resolution of desired grid in X direction.
        :param dy:
          Resolution of desired grid in Y direction.
        :param xspan:
          Width of desired grid.
        :param yspan:
          Height of desired grid.
        :param defaultVs30:
          Default Vs30 value to use if vs30File not specified.
        :param vs30File:
          Name of GMT or GDAL format grid file containing Vs30 values.
        :param vs30measured:
          Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
            Boolean indicating whether event is on the backarc as defined here: 
            http://earthquake.usgs.gov/learn/glossary/?term=backarc
        """
        geodict = GeoDict.createDictFromCenter(cx,cy,dx,dy,xspan,yspan)
        if vs30File is not None:
            try:
                vs30grid = GMTGrid.load(vs30File,samplegeodict=geodict,resample=True,method='linear')
            except Exception as msg1:
                try:
                    vs30grid = GDALGrid.load(vs30File,samplegeodict=geodict,resample=True,method='linear')
                except Exception as msg2:
                    msg = 'Load failure of %s - error messages: "%s"\n "%s"' % (vs30File,str(msg1),str(msg2))
                    raise SitesException(msg)
        else:
            griddata = np.ones((geodict.ny,geodict.nx))*defaultVs30
            vs30grid = Grid2D(griddata,geodict)
        return cls(vs30grid,vs30measured=vs30measured,backarc=backarc)

    def sampleFromSites(lats,lons):
        """Create a SitesContext object by sampling the current Sites object.

        :param lats:
          Sequence of latitudes.
        :param lons:
          Sequence of longitudes.
        :returns:
          SitesContext object where data are sampled from the current Sites object.
        :raises SitesException:
           When lat/lon input sequences do not share dimensionality.
        """
        lats = np.array(lats)
        lons = np.array(lons)
        latshape = lats.shape
        lonshape = lons.shape
        if latshape != lonshape:
            msg = 'Input lat/lon arrays must have the same dimensions'
            raise SitesException(msg)
        
        site = SitesContext()
        site.vs30 = self.Vs30.getValue(lats,lons)
        site.lats = lats
        site.lons = lons
        site.z1pt0 = calculateZ1P0(site.vs30)
        site.z2pt5 = calculateZ2P5(site.z1pt0)
        site.vs30measured = self.SitesContext.vs30measured
        site.backarc = self.SitesContext.backarc

        return site

    def getVs30Grid(self):
        """
        Return the Grid2D object containing Vs30 values for this Sites object.
        """
        return self.Vs30
    
    def getSitesContext(self):
        """Return a SitesContext object appropriate for GMPE use.
        :returns:
           SitesContext object.
        """
        return self.SitesContext
        
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
