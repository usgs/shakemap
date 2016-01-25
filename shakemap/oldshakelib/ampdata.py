#!/usr/bin/env python

#stdlib imports
from xml.dom import minidom
import sys

#third party imports
import numpy as np

#local imports
from distance import getCartesianDistance
from fault import Fault

class AmpData(object):
    def __init__(self,datafile,intensity=False,gmice=None):
        """
        Read in data XML file, turn into internal data structure.
        Ignore intensity values unless intensity flag=True, in which case
        ground motions will be ignored, and intensities will be converted to 
        ground motions using a supplied gmice object.  The AmpData object will
        be flagged as an "intensity" AmpData.
        @param datafile: ShakeMap XML data file.
        @keyword intensity: Treat this data file as an intensity data set.
        @keyword gmice: If intensity==True, use this GMICE object to convert intensities to ground motions.
        """
        root = minidom.parse(datafile)
        stations = root.getElementsByTagName('station')
        self.stationlist = []
        for station in stations:
            sdict = {}
            sdict['lat'] = float(station.getAttribute('lat'))
            sdict['lon'] = float(station.getAttribute('lat'))
            sdict['code'] = station.getAttribute('code')
            comps = station.getElementsByTagName('comp')
            pga = []
            pgv = []
            psa03 = []
            psa10 = []
            psa30 = []
            for comp in comps:
                cname = comp.getAttribute('name').strip()
                if cname.lower() == 'z' or cname.lower() == 'ud':
                    continue
                #anything else we find is presumed to be horizontal data - we want to grab all components, 
                #then take the maximum values for any ground motion
                for child in comp.childNodes:
                    if child.nodeType != child.ELEMENT_NODE:
                        continue
                    nodename = child.nodeName
                    if nodename == 'acc':
                        pga.append(float(child.getAttribute('value')))
                    if nodename == 'vel':
                        pgv.append(float(child.getAttribute('value')))
                    if nodename == 'psa03':
                        psa03.append(float(child.getAttribute('value')))
                    if nodename == 'psa10':
                        psa10.append(float(child.getAttribute('value')))
                    if nodename == 'psa30':
                        psa30.append(float(child.getAttribute('value')))
                    
            sdict['pga'] = np.max(np.abs(np.array(pga)))
            sdict['pgv'] = np.max(np.abs(np.array(pgv)))
            sdict['psa03'] = np.max(np.abs(np.array(psa03)))
            sdict['psa10'] = np.max(np.abs(np.array(psa10)))
            sdict['psa30'] = np.max(np.abs(np.array(psa30)))
            self.stationlist.append(sdict.copy())
            
        root.unlink()

    def calculateDistances(self,source):
        """
        Calculate distance from source (epicenter, fault, etc.) to all stations
        using supplied distance function.
        @param source: Either a dictionary of {'lat','lon','depth'} or a fault object. (??)
        @param distfunction: The distance function to use to calculate distances.
        """
        if isinstance(source,Fault):
            raise Exception,'Faults not currently supported'
        lat = source['lat']
        lon = source['lon']
        for station in self.stationlist:
            station['distance'] = getCartesianDistance(lat,lon,0,station['lat'],station['lon'],0)/1000.0
        self.stationlist.sort(key=lambda station: station['distance'])

    def correctToRock(self,doBasinCorrection=False):
        """
        Correct data to rock, remove site amplifications,optionally perform basin correction.
        @keyword doBasinCorrection: Apply basin correction.
        """
        pass

    def flagStations(self,stationlist):
        """
        Flag the stations that users have told us are bad.
        @param stationlist: List of station codes that should be flagged.
        """
        pass

    def calculateGMPESigma(self,sigmafunction):
        """
        Calculate GMPE sigma for each station, using uncertainty function defined by GMPE.
        @param sigmafunction: Static GMPE uncertainty function.
        """
        pass

    def calculateGMICESigma(self,sigmafunction):
        """
        Calculate GMICE sigma for each station, using uncertainty function defined by GMICE.
        @param sigmafunction: Static GMICE uncertainty function.
        """
        pass

    def calculateBias(self,optfunction,errfunction,threshold,params):
        """
        Determine the magnitude which minimizes differences between data and model, flag those stations that
        are more than threshold std away from GMPE model.
        @param optfunction: Optimization function (Simplex, etc.) to explore parameter space.
        @param errfunction: Function used to determine minimum.
        @param threshold: Number of standard deviations above which data should be flagged.
        @param params: Dictionary containing bias calculation parameters (??)
        @return: Magnitude offset representing bias.
        """
        pass

if __name__ == '__main__':
    ampfile = sys.argv[1]
    lat = float(sys.argv[2])
    lon = float(sys.argv[3])
    ampdata = AmpData(ampfile)
    source = {'lat':lat,'lon':lon}
    ampdata.calculateDistances(source)
    print ampdata.stationlist[0]['distance']
    print ampdata.stationlist[1]['distance']
    print ampdata.stationlist[-1]['distance']
