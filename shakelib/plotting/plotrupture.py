import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from mpl_toolkits.mplot3d import axes3d


MAP_WIDTH = 10.0  # map width in degrees

# Temporary function for plotting maps.
# The idea is to separate plotting methods from modules like rupture.


def plot_rupture_wire3d(rupture, ax=None):
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
        ax = axes3d.Axes3D(fig)
        # ax = fig.add_subplot(111, projection="3d")
    else:
        if "zaxis" not in list(ax.properties().keys()):
            raise TypeError("Non-3d axes object passed to plot() method.")
    for quad in rupture.getQuadrilaterals():
        x = [p.longitude for p in quad]
        x.append(x[0])
        y = [p.latitude for p in quad]
        y.append(y[0])
        z = [-p.depth for p in quad]
        z.append(z[0])
        ax.plot(x, y, z, "k")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_zlabel("Depth")

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

    # remove nans in these arrays
    rlats = rlats[~np.isnan(rlats)]
    rlons = rlons[~np.isnan(rlons)]

    minbufx = 0.2
    minbufy = 0.2
    lat1 = np.nanmin(rupture.lats)
    lat2 = np.nanmax(rupture.lats)
    dlat = lat2 - lat1
    bufy = np.max([dlat / 2, minbufy])
    lon1 = np.nanmin(rupture.lons)
    lon2 = np.nanmax(rupture.lons)
    dlon = lon2 - lon1
    bufx = np.max([dlon / 2, minbufx])

    ax = plt.axes(projection=ccrs.Robinson())
    bufx = MAP_WIDTH / 2.0
    bufy = MAP_WIDTH / 2.0
    ax.set_extent([lon1 - bufx, lon2 + bufx, lat1 - bufy, lat2 + bufy])
    ax.stock_img()
    ax.coastlines()
    ax.plot(rlons, rlats, "r", transform=ccrs.PlateCarree())
    return ax
