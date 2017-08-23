#!/usr/bin/env python

# stdlib
import os.path
import time
from functools import partial
from collections import OrderedDict

# third party
#from mapio.basemapcity import BasemapCities
from mapio.gmt import GMTGrid
import fiona
from matplotlib.patches import Polygon
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
from shapely.ops import transform
from shapely.geometry import MultiPolygon
from shapely.geometry import MultiLineString
from shapely.geometry import GeometryCollection
from shapely.geometry import Polygon as sPolygon
from shapely.geometry import shape as sShape
from shapely.geometry import mapping
from mpl_toolkits.basemap import Basemap
import numpy as np
from descartes import PolygonPatch
from scipy.ndimage import gaussian_filter
import pyproj

# local imports
from shakelib.utils.exception import ShakeLibException
from shakelib.gmice.wgrw12 import WGRW12


CITY_COLS = 2
CITY_ROWS = 2
CITIES_PER_GRID = 10
FIG_WIDTH = 8
FIG_HEIGHT = 8
BASEMAP_RESOLUTION = 'c'
WATERCOLOR = '#7AA1DA'
NCONTOURS = 6
VERT_EXAG = 0.1

# all of the zorder values for different plotted parameters
IMG_ZORDER = 1
STATIONS_ZORDER = 250
CITIES_ZORDER = 100
FAULT_ZORDER = 500
EPICENTER_ZORDER = 500
CONTOUR_ZORDER = 800
ROAD_ZORDER = 5
SCALE_ZORDER = 1500
GRATICULE_ZORDER = 1200
OCEAN_ZORDER = 1000
BORDER_ZORDER = 1001


def getProjectedPolygon(polygon, m):
    extlon, extlat = zip(*polygon.exterior.coords[:])
    extx, exty = m(extlon, extlat)
    extpts = list(zip(extx, exty))
    ints = []
    for interior in polygon.interiors:
        try:
            intlon, intlat = zip(*interior.coords[:])
        except:
            x = 1
        intx, inty = m(intlon, intlat)
        ints.append(list(zip(intx, inty)))
    ppolygon = sPolygon(extpts, ints)
    return ppolygon


def getProjectedPatches(polygon, m, edgecolor=WATERCOLOR):
    patches = []
    if isinstance(polygon, MultiPolygon):
        for p in polygon:
            ppolygon = getProjectedPolygon(p, m)
            patch = PolygonPatch(ppolygon, facecolor=WATERCOLOR, edgecolor=edgecolor,
                                 zorder=OCEAN_ZORDER, linewidth=1, fill=True, visible=True)
            patches.append(patch)
    else:
        ppolygon = getProjectedPolygon(polygon, m)
        patch = PolygonPatch(ppolygon, facecolor=WATERCOLOR, edgecolor=edgecolor,
                             zorder=OCEAN_ZORDER, linewidth=1, fill=True, visible=True)
        patches.append(patch)

    return patches


class MapMaker(object):

    def __init__(self, shakemap, topofile, stations, fault, layerdict, source, cities=None):
        req_keys = set(['coast', 'ocean', 'lake', 'country', 'state', 'roads'])
        if len(set(layerdict.keys()).intersection(req_keys)) != len(req_keys):
            raise ShakeMapException(
                'layerdict input must have all keys from %s' % str(req_keys))
        self.shakemap = shakemap
        self.topofile = topofile
        self.layerdict = layerdict
        self.cities = cities
        self.city_cols = CITY_COLS
        self.city_rows = CITY_ROWS
        self.cities_per_grid = CITIES_PER_GRID
        self.imt_layer = None
        self.contour_layer = None
        self.intensity_colormap = None
        self.contour_colormap = None
        self.stations = stations
        self.fault = fault
        self.source = source
        self.fig_width = FIG_WIDTH
        self.fig_height = FIG_HEIGHT

        # clip all the vector data now so that map rendering will be fast
        t1 = time.time()
        self._clipBounds()
        t2 = time.time()
        print('%.1f seconds to clip vectors.' % (t2 - t1))

    def _clipBounds(self):
        # returns a list of GeoJSON-like mapping objects
        xmin, xmax, ymin, ymax = self.shakemap.getBounds()
        bbox = (xmin, ymin, xmax, ymax)
        bboxpoly = sPolygon([(xmin, ymax), (xmax, ymax),
                             (xmax, ymin), (xmin, ymin), (xmin, ymax)])
        self.vectors = {}
        for key, value in self.layerdict.items():
            vshapes = []
            f = fiona.open(value, 'r')
            shapes = f.items(bbox=bbox)
            for shapeidx, shape in shapes:
                tshape = sShape(shape['geometry'])
                intshape = tshape.intersection(bboxpoly)
                vshapes.append(intshape)
            print('Filename is %s' % value)
            f.close()
            self.vectors[key] = vshapes

    def setCityGrid(self, nx=2, ny=2, cities_per_grid=10):
        self.city_cols = nx
        self.city_rows = ny
        self.cities_per_grid = cities_per_grid

    def setFigureSize(self, figwidth, figheight):
        self.fig_width = figwidth
        self.fig_height = figheight

    def setCityList(self, dataframe):
        self.cities = BaseMapCities(dataframe)  # may raise exception
        self.city_rows = None
        self.city_cols = None
        self.cities_per_grid = None

    def setIntensityLayer(self, imt_layer):
        self.imt_layer = imt_layer

    def setContourLayer(self, contour_layer):
        self.contour_layer = contour_layer

    def setIntensityGMTColorMap(self, colormap):
        self.intensity_colormap = colormap

    def setContourGMTColorMap(self, colormap):
        self.contour_colormap = colormap

    def _setMap(self, gd):
        clon = gd.xmin + (gd.xmax - gd.xmin) / 2.0
        clat = gd.ymin + (gd.ymax - gd.ymin) / 2.0
        f = plt.figure(figsize=(self.fig_width, self.fig_height))
        ax = f.add_axes([0.1, 0.1, 0.8, 0.8])

        m = Basemap(llcrnrlon=gd.xmin, llcrnrlat=gd.ymin, urcrnrlon=gd.xmax, urcrnrlat=gd.ymax,
                    rsphere=(6378137.00, 6356752.3142),
                    resolution=BASEMAP_RESOLUTION, projection='merc',
                    lat_0=clat, lon_0=clon, lat_ts=clat, ax=ax, suppress_ticks=True)
        return m

    def _projectGrid(self, data, m, gd):
        # set up meshgrid to project topo and mmi data
        xmin = gd.xmin
        if gd.xmax < gd.xmin:
            xmin -= 360
        lons = np.linspace(xmin, gd.xmax, gd.nx)
        # backwards so it plots right side up
        lats = np.linspace(gd.ymax, gd.ymin, gd.ny)
        llons1, llats1 = np.meshgrid(lons, lats)
        pdata = m.transform_scalar(np.flipud(data), lons, lats[::-1], gd.nx, gd.ny, returnxy=False,
                                   checkbounds=False, order=1, masked=False)
        return pdata

    def _getDraped(self, data, topodata):
        maxvalue = self.intensity_colormap.vmax
        mmisc = data / maxvalue
        rgba_img = self.intensity_colormap.cmap(mmisc)
        rgb = np.squeeze(rgba_img[:, :, 0:3])
        # use lightsource class to make our shaded topography
        ls = LightSource(azdeg=135, altdeg=45)
        # intensity = ls.hillshade(ptopo,fraction=0.25,vert_exag=1.0)

        ls1 = LightSource(azdeg=120, altdeg=45)
        ls2 = LightSource(azdeg=225, altdeg=45)
        intensity1 = ls1.hillshade(
            topodata, fraction=0.25, vert_exag=VERT_EXAG)
        intensity2 = ls2.hillshade(
            topodata, fraction=0.25, vert_exag=VERT_EXAG)
        intensity = intensity1 * 0.5 + intensity2 * 0.5

        draped_hsv = ls.blend_hsv(rgb, np.expand_dims(intensity, 2))

        return draped_hsv

    def _drawBoundaries(self, m):
        allshapes = self.vectors['country'] + self.vectors['state']
        for shape in allshapes:
            # shape is a geojson-like mapping thing
            blon, blat = zip(*shape.coords[:])
            bx, by = m(blon, blat)
            m.plot(bx, by, 'k', zorder=BORDER_ZORDER)

    def _drawRoads(self, m):
        allshapes = self.vectors['roads']
        xmin = 9999999
        ymin = xmin
        xmax = -99999999
        ymax = xmax
        for shape in allshapes:
            # shape is a shapely geometry
            if isinstance(shape, MultiLineString):
                blon = []
                blat = []
                for mshape in shape:
                    tlon, tlat = zip(*mshape.coords[:])
                    blon += tlon
                    blat += tlat
            else:
                blon, blat = zip(*shape.coords[:])

            if min(blon) < xmin:
                xmin = min(blon)
            if min(blat) < ymin:
                ymin = min(blat)
            if max(blon) > xmax:
                xmax = max(blon)
            if max(blat) > ymax:
                ymax = max(blat)
            bx, by = m(blon, blat)
            m.plot(bx, by, '#808080', zorder=ROAD_ZORDER)
        x = 1

    def _drawLakes(self, m, gd):
        lakes = self.vectors['lake']
        for lake in lakes:
            ppatches = getProjectedPatches(lake, m, edgecolor='k')
            for ppatch in ppatches:
                m.ax.add_patch(ppatch)

    def _drawOceans(self, m, gd):
        ocean = self.vectors['ocean'][0]  # this is one shapely polygon
        ppatches = getProjectedPatches(ocean, m)
        for ppatch in ppatches:
            m.ax.add_patch(ppatch)

    def _drawCoastlines(self, m, gd):
        coasts = self.vectors['coast']
        for coast in coasts:  # these are polygons?
            clon, clat = zip(*coast.exterior.coords[:])
            cx, cy = m(clon, clat)
            m.plot(cx, cy, 'k', zorder=BORDER_ZORDER)

    def _drawGraticules(self, m, gd):
        par = np.arange(np.ceil(gd.ymin), np.floor(gd.ymax) + 1, 1.0)
        mer = np.arange(np.ceil(gd.xmin), np.floor(gd.xmax) + 1, 1.0)
        merdict = m.drawmeridians(mer, labels=[0, 0, 0, 1], fontsize=10,
                                  linewidth=0.5, color='gray', zorder=GRATICULE_ZORDER)
        pardict = m.drawparallels(par, labels=[1, 0, 0, 0], fontsize=10,
                                  linewidth=0.5, color='gray', zorder=GRATICULE_ZORDER)

        # loop over meridian and parallel dicts, change/increase font, draw
        # ticks
        xticks = []
        for merkey, mervalue in merdict.items():
            merline, merlablist = mervalue
            merlabel = merlablist[0]
            merlabel.set_family('sans-serif')
            merlabel.set_fontsize(12.0)
            xticks.append(merline[0].get_xdata()[0])

        yticks = []
        for parkey, parvalue in pardict.items():
            parline, parlablist = parvalue
            parlabel = parlablist[0]
            parlabel.set_family('sans-serif')
            parlabel.set_fontsize(12.0)
            yticks.append(parline[0].get_ydata()[0])

        # plt.tick_params(axis='both',color='k',direction='in')
        plt.xticks(xticks, ())
        plt.yticks(yticks, ())
        m.ax.tick_params(direction='out')

    def _drawTitle(self, isContour=False):
        # Add a title
        hlon = self.shakemap.getEventDict()['lon']
        hlat = self.shakemap.getEventDict()['lat']
        edict = self.shakemap.getEventDict()
        eloc = edict['event_description']
        timestr = edict['event_timestamp'].strftime('%b %d, %Y %H:%M:%S')
        mag = edict['magnitude']
        if hlon < 0:
            lonstr = 'W%.2f' % np.abs(hlon)
        else:
            lonstr = 'E%.2f' % hlon
        if hlat < 0:
            latstr = 'S%.2f' % np.abs(hlat)
        else:
            latstr = 'N%.2f' % hlat
        dep = edict['depth']
        eid = edict['event_id']
        net = edict['event_network']
        if not eid.startswith(net):
            eid = net + eid
        tpl = (timestr, mag, latstr, lonstr, dep, eid)
        layername = 'MMI'
        if isContour:
            layername = self.contour_layer.upper()
        plt.suptitle('USGS ShakeMap (%s): %s' % (layername, eloc),
                     fontsize=14, verticalalignment='top', y=0.95)
        plt.title('%s UTC M%.1f %s %s Depth: %.1fkm ID:%s' %
                  tpl, fontsize=10, verticalalignment='bottom')
        return eid

    def _drawStations(self, m, fill=False, imt='pga'):
        # get the locations and values of the MMI observations
        mmidf = self.stations.getStationDataframe(0)
        # get the locations and values of the instrumented observations
        instdf = self.stations.getStationDataframe(1)
        if imt in ('mmi', 'pga', 'pgv'):
            dimt = imt.upper()
        else:
            perstr = imt.replace('psa', '')
            dimt = 'SA(' + perstr[0] + '.' + perstr[1] + ')'
        if not fill:
            # plot MMI as small circles
            mmilat = mmidf['lat']
            mmilon = mmidf['lon']
            m.plot(mmilon, mmilat, 'ko', latlon=True, fillstyle='none',
                   markersize=4, zorder=STATIONS_ZORDER)

            # plot MMI as slightly larger triangles
            instlat = instdf['lat']
            instlon = instdf['lon']
            m.plot(instlon, instlat, 'k^', latlon=True, fillstyle='none',
                   markersize=6, zorder=STATIONS_ZORDER)
        else:
            for idx, value in enumerate(mmidf['lat']):
                mlat = mmidf['lat'][idx]
                mlon = mmidf['lon'][idx]
                mmi = mmidf['MMI'][idx]
                mcolor = self.intensity_colormap.getDataColor(mmi)
                m.plot(mlon, mlat, 'o', latlon=True,
                       markerfacecolor=mcolor, markeredgecolor='k', 
                       markersize=4, zorder=STATIONS_ZORDER)

            for idx, value in enumerate(instdf['lat']):
                mlat = instdf['lat'][idx]
                mlon = instdf['lon'][idx]
                dmmi = instdf[dimt][idx]
                mcolor = self.intensity_colormap.getDataColor(dmmi)
                m.plot(mlon, mlat, '^', latlon=True,
                       markerfacecolor=mcolor, markeredgecolor='k', 
                       markersize=6, zorder=STATIONS_ZORDER)

    def _drawFault(self, m):
        lats = self.fault.lats
        lons = self.fault.lons
        x, y = m(lons, lats)
        m.plot(x, y, 'k', lw=2, zorder=FAULT_ZORDER)

    def drawIntensityMap(self, outfolder):
        if self.intensity_colormap is None:
            raise ShakeMapException(
                'MapMaker.setGMTColormap() has not been called.')
        t0 = time.time()
        # resample shakemap to topogrid
        # get the geodict for the topo file
        topodict = GMTGrid.getFileGeoDict(self.topofile)[0]
        # get the geodict for the ShakeMap
        smdict = self.shakemap.getGeoDict()
        # get a geodict that is aligned with topo, but inside shakemap
        sampledict = topodict.getBoundsWithin(smdict)

        self.shakemap = self.shakemap.interpolateToGrid(sampledict)

        gd = self.shakemap.getGeoDict()

        # establish the basemap object
        m = self._setMap(gd)

        # get topo layer and project it
        topogrid = GMTGrid.load(
            self.topofile, samplegeodict=sampledict, resample=False)
        topodata = topogrid.getData().copy()
        ptopo = self._projectGrid(topodata, m, gd)

        # get intensity layer and project it
        imtdata = self.shakemap.getLayer(self.imt_layer).getData().copy()
        pimt = self._projectGrid(imtdata, m, gd)

        # get the draped intensity data
        draped_hsv = self._getDraped(pimt, ptopo)  # where will 10.0 come from

        # draw the draped intensity data
        m.imshow(draped_hsv, interpolation='none', zorder=IMG_ZORDER)

        # draw country/state boundaries
        self._drawBoundaries(m)

        # draw whatever road data is available
        self._drawRoads(m)

        # draw lakes
        self._drawLakes(m, gd)

        # draw oceans (pre-processed with islands taken out)
        t1 = time.time()
        self._drawOceans(m, gd)
        t2 = time.time()
        print('%.1f seconds to render oceans.' % (t2 - t1))

        # draw coastlines
        self._drawCoastlines(m, gd)

        # draw meridians, parallels, labels, ticks
        self._drawGraticules(m, gd)

        # draw map scale
        scalex = gd.xmin + (gd.xmax - gd.xmin) / 5.0
        scaley = gd.ymin + (gd.ymax - gd.ymin) / 10.0
        yoff = (0.007 * (m.ymax - m.ymin))
        clon = (gd.xmin + gd.xmax) / 2.0
        clat = (gd.ymin + gd.ymax) / 2.0
        m.drawmapscale(scalex, scaley, clon, clat, length=100,
                       barstyle='fancy', yoffset=yoff, zorder=SCALE_ZORDER)

        # draw fault polygon, if present
        self._drawFault(m)  # get the fault loaded

        # draw epicenter
        hlon = self.shakemap.getEventDict()['lon']
        hlat = self.shakemap.getEventDict()['lat']
        m.plot(hlon, hlat, 'k*', latlon=True, fillstyle='none',
               markersize=22, mew=1.2, zorder=EPICENTER_ZORDER)

        # draw cities
        # reduce the number of cities to those whose labels don't collide
        # set up cities
        if self.city_cols is not None:
            self.cities = self.cities.limitByBounds(
                (gd.xmin, gd.xmax, gd.ymin, gd.ymax))
            self.cities = self.cities.limitByGrid(nx=self.city_cols, ny=self.city_rows,
                                                  cities_per_grid=self.cities_per_grid)
            print("Available fonts: ", self.cities._fontlist)
            self.cities = self.cities.limitByMapCollision(m, fontname='Times New Roman')
        self.cities.renderToMap(m.ax, zorder=CITIES_ZORDER)

        # draw title and supertitle
        eventid = self._drawTitle()

        # draw station and macroseismic locations
        self._drawStations(m)  # need stationlist object

        # save plot to file
        plt.draw()
        outfile = os.path.join(outfolder, 'intensity_%s.pdf' % eventid)
        plt.savefig(outfile)
        tn = time.time()
        print('%.1f seconds to render entire map.' % (tn - t0))
        return outfile

    def _getShaded(self, ptopo):
        maxvalue = self.contour_colormap.vmax
        ls1 = LightSource(azdeg=120, altdeg=45)
        ls2 = LightSource(azdeg=225, altdeg=45)
        intensity1 = ls1.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
        intensity2 = ls2.hillshade(ptopo, fraction=0.25, vert_exag=VERT_EXAG)
        intensity = intensity1 * 0.5 + intensity2 * 0.5

        ptoposc = ptopo / maxvalue
        rgba = self.contour_colormap.cmap(ptoposc)
        rgb = np.squeeze(rgba)

        draped_hsv = ls1.blend_hsv(rgb, np.expand_dims(intensity, 2))

        return draped_hsv

    def round_to(self, n, precision):
        correction = 0.5 if n >= 0 else -0.5
        return int(n / precision + correction) * precision

    def getContourLevels(self, dmin, dmax, imt):
        # groupings taken from table on
        # https://en.wikipedia.org/wiki/Peak_ground_acceleration
        if imt == 'pgv':
            # table of minimum dmax and dinc levels
            dmax_dinc = OrderedDict([(1.1, 0.1),
                                     (3.4, 0.25),
                                     (8.1, 0.5),
                                     (16.0, 2.0),
                                     (31.0, 5.0),
                                     (60.0, 10.0),
                                     (116.0, 10.0),
                                     (200.0, 25.0)])
            keys = np.array(list(dmax_dinc.keys()))
            didx = np.where(keys < dmax)[0].max()
            dinc = dmax_dinc[keys[didx]]
            newdmin = self.round_to(dmin, dinc)
            newdmax = self.round_to(dmax, dinc)
        else:
            dmax_dinc = OrderedDict([(0.014 * 100, 0.1),
                                     (0.039 * 100, 0.5),
                                     (0.092 * 100, 1.0),
                                     (0.18 * 100, 2.5),
                                     (0.34 * 100, 5.0),
                                     (0.65 * 100, 10.0),
                                     (1.24 * 100, 15.0),
                                     (3.0 * 100, 37.5)])
            keys = np.array(list(dmax_dinc.keys()))
            didx = np.where(keys < dmax)[0].max()
            dinc = dmax_dinc[keys[didx]]
            newdmin = self.round_to(dmin, dinc)
            newdmax = self.round_to(dmax, dinc)
        levels = np.arange(newdmin, newdmax + dinc, dinc)
        return levels

    def drawContourMap(self, outfolder, cmin=None, cmax=None):
        if self.contour_colormap is None:
            raise ShakeMapException(
                'MapMaker.setGMTColormap() has not been called.')
        t0 = time.time()
        # resample shakemap to topogrid
        # get the geodict for the topo file
        topodict = GMTGrid.getFileGeoDict(self.topofile)[0]
        # get the geodict for the ShakeMap
        smdict = self.shakemap.getGeoDict()
        # get a geodict that is aligned with topo, but inside shakemap
        sampledict = topodict.getBoundsWithin(smdict)

        self.shakemap = self.shakemap.interpolateToGrid(sampledict)

        gd = self.shakemap.getGeoDict()

        # establish the basemap object
        m = self._setMap(gd)

        # get topo layer and project it
        topogrid = GMTGrid.load(
            self.topofile, samplegeodict=sampledict, resample=False)
        topodata = topogrid.getData().copy()
        ptopo = self._projectGrid(topodata, m, gd)

        # get contour layer and project it1
        imtdata = self.shakemap.getLayer(self.contour_layer).getData().copy()
        pimt = self._projectGrid(imtdata, m, gd)

        # get the draped intensity data
        hillshade = self._getShaded(ptopo)

        # draw the draped intensity data
        m.imshow(hillshade, interpolation='none', zorder=IMG_ZORDER)

        # draw the contours of imt data
        xmin = gd.xmin
        if gd.xmax < gd.xmin:
            xmin -= 360
        lons = np.linspace(xmin, gd.xmax, gd.nx)
        # backwards so it plots right side up
        lats = np.linspace(gd.ymax, gd.ymin, gd.ny)
        x, y = m(*np.meshgrid(lons, lats))
        pimt = gaussian_filter(pimt, 5.0)
        dmin = pimt.min()
        dmax = pimt.max()
        levels = self.getContourLevels(dmin, dmax, self.contour_layer)
        cs = m.contour(x, y, np.flipud(pimt), colors='w',
                       cmap=None, levels=levels, zorder=CONTOUR_ZORDER)
        clabels = plt.clabel(cs, colors='k', fmt='%.1f',
                             fontsize=8.0, zorder=CONTOUR_ZORDER)
        for cl in clabels:
            bbox = dict(boxstyle="round", facecolor='white', edgecolor='w')
            cl.set_bbox(bbox)
            cl.set_zorder(CONTOUR_ZORDER)

        # draw country/state boundaries
        self._drawBoundaries(m)

        # draw lakes
        self._drawLakes(m, gd)

        # draw oceans (pre-processed with islands taken out)
        t1 = time.time()
        self._drawOceans(m, gd)
        t2 = time.time()
        print('%.1f seconds to render oceans.' % (t2 - t1))

        # draw coastlines
        self._drawCoastlines(m, gd)

        # draw meridians, parallels, labels, ticks
        self._drawGraticules(m, gd)

        # draw filled symbols for MMI and instrumented measures
        self._drawStations(m, fill=True, imt=self.contour_layer)

        # draw map scale
        scalex = gd.xmin + (gd.xmax - gd.xmin) / 5.0
        scaley = gd.ymin + (gd.ymax - gd.ymin) / 10.0
        yoff = (0.007 * (m.ymax - m.ymin))
        clon = (gd.xmin + gd.xmax) / 2.0
        clat = (gd.ymin + gd.ymax) / 2.0
        m.drawmapscale(scalex, scaley, clon, clat, length=100,
                       barstyle='fancy', yoffset=yoff, zorder=SCALE_ZORDER)

        # draw fault polygon, if present
        self._drawFault(m)  # get the fault loaded

        # draw epicenter
        hlon = self.shakemap.getEventDict()['lon']
        hlat = self.shakemap.getEventDict()['lat']
        m.plot(hlon, hlat, 'k*', latlon=True, fillstyle='none',
               markersize=22, mew=1.2, zorder=EPICENTER_ZORDER)

        # draw cities
        # reduce the number of cities to those whose labels don't collide
        # set up cities
        if self.city_cols is not None:
            self.cities = self.cities.limitByBounds(
                (gd.xmin, gd.xmax, gd.ymin, gd.ymax))
            self.cities = self.cities.limitByGrid(nx=self.city_cols, ny=self.city_rows,
                                                  cities_per_grid=self.cities_per_grid)
            self.cities = self.cities.limitByMapCollision(m, fontname='Times New Roman')
        self.cities.renderToMap(m.ax, zorder=CITIES_ZORDER)

        # draw title and supertitle
        eventid = self._drawTitle(isContour=True)

        # draw whatever road data is available
        # self._drawRoads(m)

        # save plot to file
        plt.draw()
        outfile = os.path.join(outfolder, 'contour_%s_%s.pdf' %
                               (self.contour_layer, eventid))
        plt.savefig(outfile)
        tn = time.time()
        print('%.1f seconds to render entire map.' % (tn - t0))
        return outfile
