#!/usr/bin/env python

# stdlib imports
import os.path
import sys

# third party

# hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))
# put this at the front of the system path, ignoring any installed mapio stuff
sys.path.insert(0, shakedir)

from shakemap.mapping.mapmaker import MapMaker
from impactutils.colors.cpalette import ColorPalette

from mapio.basemapcity import BasemapCities
from mapio.shake import ShakeGrid

# local imports
from shakemap.grind.station import StationList
from shakemap.grind.fault import Fault
from shakemap.grind.source import Source


def _test_intensity():

    datadir = os.path.abspath(os.path.join(
        homedir, '..', 'data', 'eventdata', 'northridge'))
    shakefile = os.path.join(datadir, 'northridge_grid.xml')
    topofile = os.path.join(datadir, 'northridge_topo.grd')
    faultfile = os.path.join(datadir, 'northridge_fault.txt')
    cityfile = os.path.join(datadir, 'northridge_cities.txt')
    coastfile = os.path.join(datadir, 'northridge_coastline.json')
    countryfile = os.path.join(datadir, 'northridge_countries.json')
    statefile = os.path.join(datadir, 'northridge_states.json')
    lakefile = os.path.join(datadir, 'northridge_lakes.json')
    oceanfile = os.path.join(datadir, 'northridge_ocean.json')
    stationfile = os.path.join(datadir, 'northridge_stations.db')
    roadfile = os.path.join(datadir, 'northridge_roads.json')
    tancptfile = os.path.join(shakedir, 'shakemap', 'mapping', 'tan.cpt')
    shakecptfile = os.path.join(
        shakedir, 'shakemap', 'mapping', 'shakecpt.cpt')

    layerdict = {'coast': coastfile,
                 'ocean': oceanfile,
                 'lake': lakefile,
                 'country': countryfile,
                 'roads': roadfile,
                 'state': statefile}

    tancolormap = ColorPalette.fromPreset('shaketopo')
    shakecolormap = ColorPalette.fromPreset('mmi')
    cities = BasemapCities.loadFromCSV(cityfile)
    shakemap = ShakeGrid.load(shakefile, adjust='res')
    stations = StationList(stationfile)
    fault = Fault.readFaultFile(faultfile)
    edict = shakemap.getEventDict()
    eventdict = {'lat': edict['lat'],
                 'lon': edict['lon'],
                 'depth': edict['depth'],
                 'mag': edict['magnitude'],
                 'time': edict['event_timestamp']}
    source = Source(eventdict, fault)
    maker = MapMaker(shakemap, topofile, stations,
                     fault, layerdict, source, cities)

    # draw intensity map
    outfolder = os.path.expanduser('~')
    maker.setIntensityLayer('mmi')
    maker.setIntensityGMTColorMap(shakecolormap)
    intensity_map = maker.drawIntensityMap(outfolder)
    print('Intensity map saved as: %s' % intensity_map)

    # draw contour maps
    maker.setContourGMTColorMap(tancolormap)

    # Draw pgv contours
    maker.setContourLayer('pgv')
    contour_pgv_map = maker.drawContourMap(outfolder)
    print('PGV contour map saved as: %s' % contour_pgv_map)

    # Draw pga contours
    maker.setContourLayer('pga')
    contour_pga_map = maker.drawContourMap(outfolder)
    print('PGA contour map saved as: %s' % contour_pga_map)

    # Draw psa0.3 contours
    maker.setContourLayer('psa03')
    contour_psa03_map = maker.drawContourMap(outfolder)
    print('PSA0.3 contour map saved as: %s' % contour_psa03_map)

    # Draw psa1.0 contours
    maker.setContourLayer('psa10')
    contour_psa10_map = maker.drawContourMap(outfolder)
    print('PSA1.0 contour map saved as: %s' % contour_psa10_map)

    # Draw psa3.0 contours
    maker.setContourLayer('psa30')
    contour_psa30_map = maker.drawContourMap(outfolder)
    print('PSA3.0 contour map saved as: %s' % contour_psa30_map)

if __name__ == '__main__':
    _test_intensity()
