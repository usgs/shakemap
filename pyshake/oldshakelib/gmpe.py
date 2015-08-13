#!/usr/bin/env python

#third party imports
import numpy as np

DEFAULT_LATSPAN = 2.5 #dd
DEFAULT_LONSPAN = 4.0 #dd
DEFAULT_XRES = 0.008
DEFAULT_YRES = 0.008

class GMPE(object):
    """
    Base class for all GMPE objects.
    """
    #public methods
    def __init__(self,sourceparams,fault,sitecorrect,doPSA=True):
        """Instantiate the base class.
        Among other as-yet unspecified functionality, the constructor will create an 
        empty 3-D grid object using a set of rules to define the extent and resolution.
         
        @param sourceparams: A sourceparams object (from event.xml file)
        @param fault: A fault object (from fault file).  Can be None.
        @param sitecorrect: A site correction object.  Can be None.
        @keyword doPSA: If True, calculate 0.3 s, 1.0 s, and 3.0 s spectral accelerations.
        """
        self.lat = sourceparams['lat']
        self.lon = sourceparams['lon']
        self.depth = sourceparams['depth']
        self.mag = sourceparams['mag']

        self.fault = fault

        self.siteCorrection = sitecorrect

        self.doPSA = doPSA
        
        #figure out grid reference information
        self.xmin = self.lon - DEFAULT_LONSPAN/2.0
        self.ymax = self.lat + DEFAULT_LATSPAN/2.0
        self.xmax = self.lon + DEFAULT_LONSPAN/2.0
        self.ymin = self.lat - DEFAULT_LATSPAN/2.0
        self.ncols = int((self.xmax - self.xmin)/DEFAULT_XRES)
        self.nrows = int((self.ymax - self.ymin)/DEFAULT_YRES)
        self.xres = DEFAULT_XRES
        self.yres = DEFAULT_YRES

        #create empty grids for our modeled values
        self.pga = np.zeros((self.nrows,self.ncols))
        self.pgv = np.zeros((self.nrows,self.ncols))
        self.mmi = np.zeros((self.nrows,self.ncols))
        if doPSA:
            self.psa03 = np.zeros((self.nrows,self.ncols))
            self.psa10 = np.zeros((self.nrows,self.ncols))
            self.psa30 = np.zeros((self.nrows,self.ncols))

    def getGridDef(self): 
        """Return the grid definition tuple as
        determined from the input source params. 
        Return: (xmin,ymax, xres, yres,nrows,ncols) where:
                xmin: Upper left X coordinate of the grid.
                ymax: Upper left Y coordinate of the grid.
                xres: X resolution (dd)
                yres: Y resolution (dd)
                nrows: Number of cells in the Y direction.
                ncols: Number of cells in the Y direction.
        """
        return (self.xmin,self.ymax,self.xres,self.yres,self.nrows,self.ncols)

    def setGridDef(self,xmin,ymax,xres,yres,nrows,ncols):
        """Override the grid definition tuple as
        determined from the input source params. 
        
        @param xmin: Tuple of (xmin,ymax)
        @param xres: X resolution (dd)
        @param yres: Y resolution (dd)
        @param nrows: Number of cells in the Y direction.
        @param ncols: Number of cells in the Y direction.
        """
        self.xmin = xmin
        self.ymax = ymax
        self.nrows = nrows
        self.ncols = ncols
        self.xres = xres
        self.yres = yres
        self.xmax = self.xmin + self.ncols*self.xres
        self.ymin = self.ymax - self.nrows*self.yres

    def getRowCol(self,lat,lon):
        """Return data row and column from given geographic coordinates (lat/lon decimal degrees).
        @param lat: Input latitude.
        @param lon: Input longitude.
        @return: Tuple of row and column.
        """
        ulx = self.xmin
        uly = self.ymax
        xdim = self.xres
        ydim = self.yres
        col = np.floor((lon-ulx)/xdim)
        row = np.floor((uly-lat)/ydim)
        return (row,col)

    def getLatLon(self,row,col):
        """Return geographic coordinates (lat/lon decimal degrees) for given data row and column.
        @param row: Row dimension index into internal data array.
        @param col: Column dimension index into internal data array.
        @return: Tuple of latitude and longitude.
        """
        ulx = self.xmin
        uly = self.ymax
        xdim = self.xres
        ydim = self.yres
        lon = ulx + col*xdim
        lat = uly - row*ydim
        return (lat,lon)
        
    def setLonOffset(self,lonoffset):
        """
        Modify the longitude offset of the center of the map relative to epicenter.
        @param lonoffset: Longitude offset.
        """
        pass

    def setLatOffset(self,latoffset):
        """
        Modify the latitude offset of the center of the map relative to epicenter.
        @param latoffset: Latitude offset.
        """
        pass

    def setLonSpan(self,lonspan):
        """
        Modify the longitude extent of the map.
        @param lonspan: Longitude extent.
        """
        

    def setLatSpan(self,latspan):
        """
        Modify the latitude extent of the map.
        @param latspan: Latitude extent.
        """
        pass

    def setLatLonRation(self,latlonratio):
        """
        Modify the longitude extent of the map based on latitude span and new lonratio.
        @param latlonratio: Latitude/Longitude ratio.
        """
        pass
        
    def setAmplitudes(self,ampdata):
        """
        Set ground truth amplitude data in GMPE object.
        @param ampdata: Data object. 
        """
        pass

    def setIntensities(self,mmidata):
        """
        Set intensity-derived ground motion data in GMPE object.
        @param mmidata: Data object. 
        """
        pass

    def calculateAmplitudes(self,Vs30):
        """
        Calculate amplitudes for every cell in grid given GMPE definition in paper. 
        Implemented in sub-classes.
        """
        pass

    def applySiteCorrection(self):
        """Apply site correction data as supplied in constructor.
        Implemented in sub-classes (?)
        """
        pass

    def adjustModelToData(self,ampdata):
        """
        Adjust the modeled amplitudes by a weighted average of data values (raw pga/pgv and mmi-derived both).
        Implemented in base class (?) 
        """
        

    def calculateSigma(self):
        """
        Calculate uncertainties for modeled data.
        Implemented in sub-class.
        """
        pass

    def calculateDirectivity(self):
        """
        Calculate directivity for modeled data.
        Implemented in sub-class (?)
        """
        pass

    def getPGA(self):
        """
        Return PGA layer.
        @return: 2-D numpy array containing the PGA modeled/data values.
        """

    def getPGV(self):
        """
        Return PGV layer.
        @return: 2-D numpy array containing the PGV modeled/data values.
        """

    def getMMI(self):
        """
        Return MMI layer.
        @return: 2-D numpy array containing the MMI modeled/data values.
        """

    def getPSA03(self):
        """
        Return PSA03 layer.
        @return: 2-D numpy array containing the PSA03 modeled/data values.
        @raise Exception: If GMPE was created without PSA values.
        """

    def getPSA10(self):
        """
        Return PSA10 layer.
        @return: 2-D numpy array containing the PSA10 modeled/data values.
        @raise Exception: If GMPE was created without PSA values.
        """

    def getPSA30(self):
        """
        Return PSA30 layer.
        @return: 2-D numpy array containing the PSA30 modeled/data values.
        @raise Exception: If GMPE was created without PSA values.
        """
    
    #private methods
    def _calculateDistanceToSource(self):
        """
        Calculate distance from each cell to the source, in manner appropriate for GMPE
        Implemented in sub-classes, using global vectorized functions for calculating distances.
        """
        pass

    def _calculateIntensity(self,gmice):
        """
        Calculate intensity for all grid cells, using input GMICE object.
        Implemented in base class.
        @param gmice: GMICE object used to calculate MMI from PGA/PGV.
        """
        pass

    #static methods - that is, methods that are associated with the class but do not require an 
    #instantiated object.
    @staticmethod
    def uncertaintyFunction():
        """
        A static class method implemented in sub-classes, providing a function for calculating uncertainties.
        This will be used both by the sub-class itself, and the Data object for calculating uncertainties.
        NB: Not entirely sure about the syntax of the @staticmethod decorator... putting it in the base class
        may not work.  Experimentation required.
        """
        pass
    
