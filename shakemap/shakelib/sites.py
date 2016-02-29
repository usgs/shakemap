#!/usr/bin/env python

#stdlib imports
import sys

#third party imports
from mapio.gmt import GMTGrid
from mapio.gdal import GDALGrid
from mapio.grid2d import Grid2D
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

def _load(vs30File,samplegeodict=None,resample=False,method='linear',doPadding=False,padValue=np.nan):
    try:
        vs30grid = GMTGrid.load(vs30File,samplegeodict=samplegeodict,resample=resample,
                                method=method,doPadding=doPadding,padValue=padValue)
    except Exception as msg1:
        try:
            vs30grid = GDALGrid.load(vs30File,samplegeodict=samplegeodict,resample=resample,
                                     method=method,doPadding=doPadding,padValue=padValue)
        except Exception as msg2:
            msg = 'Load failure of %s - error messages: "%s"\n "%s"' % (vs30File,str(msg1),str(msg2))
            raise SitesException(msg)
    return vs30grid

def _getFileGeoDict(fname):
    geodict = None
    try:
        geodict = GMTGrid.getFileGeoDict(fname)
    except Exception as msg1:
        try:
            geodict = GDALGrid.getFileGeoDict(fname)
        except Exception as msg2:
            msg = 'File geodict failure with %s - error messages: "%s"\n "%s"' % (fname,str(msg1),str(msg2))
            raise SitesException(msg)
    return geodict
    
def calculateZ1P0(vs30):
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
    def __init__(self,vs30grid,vs30measured=False,backarc=False,defaultVs30=686.0):
        """
        Construct a Sites object.
        :param vs30grid:
            MapIO Grid2D object containing Vs30 values.
        :param vs30measured:
            Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
            Boolean indicating whether event is on the backarc as defined here: 
            http://earthquake.usgs.gov/learn/glossary/?term=backarc
        :param defaultVs30:
          Default Vs30 value to use in locations where Vs30Grid is not specified.
        """
        self.Vs30 = vs30grid
        self.defaultVs30 = defaultVs30
        self.GeoDict = vs30grid.getGeoDict().copy()
        lons = np.arange(self.GeoDict.xmin,self.GeoDict.xmax,self.GeoDict.dx)
        lats = np.arange(self.GeoDict.ymin,self.GeoDict.ymax,self.GeoDict.dy)
        self.Z1Pt0 = calculateZ1P0(self.Vs30.getData())
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
    def _create(cls,geodict,defaultVs30,vs30File,padding,resample):
        if vs30File is not None:
            fgeodict = _getFileGeoDict(vs30File)
            if not resample:
                if not padding:
                    geodict = fgeodict.getBoundsWithin(geodict) #we want something that is within and aligned
                else:
                    geodict = fgeodict.getAligned(geodict) #we want something that is just aligned, since we're padding edges
            vs30grid = _load(vs30File,samplegeodict=geodict,resample=resample,
                             method='linear',doPadding=padding,padValue=defaultVs30)

        return vs30grid
        
    @classmethod
    def createFromBounds(cls,xmin,xmax,ymin,ymax,dx,dy,defaultVs30=686.0,vs30File=None,
                         vs30measured=False,backarc=False,padding=False,resample=False):
        """Create a Sites object by defining a center point, resolution, extent, and Vs30 values.

        :param xmin:
          X coordinate of left edge of bounds.
        :param xmax:
          X coordinate of right edge of bounds.
        :param ymin:
          Y coordinate of bottom edge of bounds.
        :param ymax:
          Y coordinate of top edge of bounds.
        :param dx:
          Resolution of desired grid in X direction.
        :param dy:
          Resolution of desired grid in Y direction.
        :param defaultVs30:
          Default Vs30 value to use if vs30File not specified.
        :param vs30File:
          Name of GMT or GDAL format grid file containing Vs30 values.
        :param vs30measured:
          Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
          Boolean indicating whether event is on the backarc as defined here: 
          http://earthquake.usgs.gov/learn/glossary/?term=backarc
        :param padding:
          Boolean indicating whether or not to pad resulting Vs30 grid out to edges of input
          bounds.  If False, grid will be clipped to the extent of the input file.
        :param resample:
          Boolean indicating whether or not the grid should be resampled.
        """
        geodict = GeoDict.createDictFromBox(xmin,xmax,ymin,ymax,dx,dy)
        if vs30File is not None:
            vs30grid = cls._create(geodict,defaultVs30,vs30File,padding,resample)
        else:
            griddata = np.ones((geodict.ny,geodict.nx))*defaultVs30
            vs30grid = Grid2D(griddata,geodict)
        return cls(vs30grid,vs30measured=vs30measured,backarc=backarc,defaultVs30=defaultVs30)
        
    @classmethod
    def createFromCenter(cls,cx,cy,xspan,yspan,dx,dy,defaultVs30=686.0,vs30File=None,
                         vs30measured=False,backarc=False,padding=False,resample=False):
        """Create a Sites object by defining a center point, resolution, extent, and Vs30 values.

        :param cx:
          X coordinate of desired center point.
        :param cy:
          X coordinate of desired center point.
        :param xspan:
          Width of desired grid.
        :param yspan:
          Height of desired grid.
        :param dx:
          Resolution of desired grid in X direction.
        :param dy:
          Resolution of desired grid in Y direction.
        :param defaultVs30:
          Default Vs30 value to use if vs30File not specified.
        :param vs30File:
          Name of GMT or GDAL format grid file containing Vs30 values.
        :param vs30measured:
          Boolean indicating whether Vs30 values were measured or derived (i.e., from slope)
        :param backarc:
          Boolean indicating whether event is on the backarc as defined here: 
          http://earthquake.usgs.gov/learn/glossary/?term=backarc
        :param padding:
          Boolean indicating whether or not to pad resulting Vs30 grid out to edges of input
          bounds.  If False, grid will be clipped to the extent of the input file.
        :param resample:
          Boolean indicating whether or not the grid should be resampled.
        """
        geodict = GeoDict.createDictFromCenter(cx,cy,dx,dy,xspan,yspan)
        if vs30File is not None:
            vs30grid = cls._create(geodict,defaultVs30,vs30File,padding,resample)
        else:
            griddata = np.ones((geodict.ny,geodict.nx))*defaultVs30
            vs30grid = Grid2D(griddata,geodict)
        return cls(vs30grid,vs30measured=vs30measured,backarc=backarc,defaultVs30=defaultVs30)

    def sampleFromSites(self,lats,lons):
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
        site.vs30 = self.Vs30.getValue(lats,lons,default=self.defaultVs30) #use default vs30 if outside grid
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
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 3.0
    yspan = 3.0
    mysite = Sites.createFromCenter(cx,cy,xspan,yspan,dx,dy,vs30File=vs30file,
                                    padding=True,resample=False)
    sc = mysite.getSitesContext()
    
    cx = -118.2
    cy = 83
    dx = 0.0083
    dy = 0.0083
    xspan = 3.0
    yspan = 3.0
    mysite = Sites.createFromCenter(cx,cy,xspan,yspan,dx,dy,vs30File=vs30file,padding=True,resample=False)

    xmin = 116.234
    xmax = 120.876
    ymin = 20.12345
    ymax = 24.75435
    dx = 0.0083
    dy = 0.0083
    mysite = Sites.createFromBounds(xmin,xmax,ymin,ymax,dx,dy,vs30File=vs30file,padding=False,resample=False)

if __name__ == '__main__':
    vs30file = None
    if len(sys.argv) > 1:
        vs30file = sys.argv[1]
    _test(vs30file=vs30file)
