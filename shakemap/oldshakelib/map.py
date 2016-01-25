#!/usr/bin/env python

#third party imports
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.colors import ListedColormap,LinearSegmentedColormap,Normalize,BoundaryNorm
from matplotlib.colors import LightSource

from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np

from neicio.gmt import GMTGrid
from neicio.shake import ShakeGrid
from neicutil.colormap import GMTColormap

def main(gridfile,topofile):
    WATER = [ 0.50980392,  0.60784314,  0.88235294]

    cmap = GMTColormap('shakecpt.cpt')
    clist = cmap.getColorList()
    boundaries = cmap.getZValues()
    palette = ListedColormap(clist,'my_colormap')

    mmi = ShakeGrid('grid.xml',variable='MMI')
    topo = GMTGrid('w100n40.grd')

    topo.interpolateToGrid(mmi.getGeoDict())

    mmidata = flipud(mmi.griddata.copy())
    topodata = topo.griddata.copy()

    bounds = mmi.getRange()

    bounds = list(bounds)

    if bounds[1] < 0 and bounds[0] > bounds[1]:
        bounds[1] = bounds[1] + 360

    clat = bounds[2] + (bounds[3] - bounds[2])/2
    clon = bounds[0] + (bounds[1] - bounds[0])/2
    dx = (bounds[1] - bounds[0])*111191 * np.cos(np.radians(clat))
    dy = (bounds[3] - bounds[2])*111191
    figwidth=6
    aspect = dy/dx
    figheight = aspect * figwidth
    fig = figure(figsize=(figwidth,figheight),edgecolor='g',facecolor='g')
    ax1 = fig.add_axes([0,0,1.0,1.0])
    m = Basemap(llcrnrlon=bounds[0],llcrnrlat=bounds[2],
                         urcrnrlon=bounds[1],urcrnrlat=bounds[3],
                         resolution='h',projection='merc',lat_ts=clat)

    # attach new axes image to existing Basemap instance.
    m.ax = ax1
    # create light source object.

    ls = LightSource(azdeg = 90, altdeg = 20)
    # convert data to rgb array including shading from light source.
    # (must specify color map)
    rgb = ls.shade(topodata, cm.binary)
    im = m.imshow(rgb)
    #hold(True)
    print 'Hello'
    jet()
    im2 = m.imshow(mmidata,alpha=0.3)
    m.drawmapboundary(linewidth=2.0)
    #m.drawlsmask(ocean_color=[ 0.50980392,  0.60784314,  0.88235294])
    m.drawrivers(color=WATER)
    savefig('output.png')
    
if __name__ == '__main__':
    gridfile = sys.argv[1]
    topofile = sys.argv[2]
    main(gridfile,topofile)
