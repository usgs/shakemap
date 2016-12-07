import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.basemap import Basemap

# Temporary function for plotting maps.
# The idea is to separate plotting methods from modules like rupture.

def plot_rupture_wire3d(rupture, ax = None):
    """
    Method for making a simple representation of a Rupture instance. 
    This method draws the outline of each quadrilateral in 3D. 

    Args:
        rupture: A Rupture instance. 
        ax: A matplotlib axis (optional). 

    Returns:
        Matplotlib axis. 
    """

    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    else:
        if 'xlim3d' not in list(ax.properties().keys()):
            raise ShakeMapException(
                'Non-3d axes object passed to plot() method.')
    for quad in rupture.getQuadrilaterals():
        x = [p.longitude for p in quad]
        x.append(x[0])
        y = [p.latitude for p in quad]
        y.append(y[0])
        z = [-p.depth for p in quad]
        z.append(z[0])
        ax.plot(x, y, z, 'k')
    return ax

def map_rupture(rupture):
    """
    Method for making a simple representation of a Rupture instance.
    This method draws the surface projection of the rupture on a map.

    Args:
        rupture: A Rupture instance.

    """
    rlats = rupture.lats
    rlons = rupture.lons

    minbufx = 0.2
    minbufy = 0.2
    lat1 = np.nanmin(rupture.lats)
    lat2 = np.nanmax(rupture.lats)
    dlat = lat2 - lat1
    bufy = np.max([dlat/2, minbufy])
    lon1 = np.nanmin(rupture.lons)
    lon2 = np.nanmax(rupture.lons)
    dlon = lon2 - lon1
    bufx = np.max([dlon/2, minbufx])
    m = Basemap(llcrnrlat=lat1-bufy,urcrnrlat=lat2+bufy, llcrnrlon=lon1-bufy,urcrnrlon=lon2+bufy)
    m.arcgisimage(service='World_Shaded_Relief', xpixels = 500, verbose= True)
    x, y = m(rlons, rlats)
    m.plot(x, y, 'r')





