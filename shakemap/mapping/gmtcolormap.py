#!/usr/bin/env python

# stdlib imports
import os.path
import colorsys

from matplotlib.colors import LinearSegmentedColormap
import numpy as np


def read_cpt(cptfile):
    f = open(cptfile, 'rt')
    lines = f.readlines()
    f.close()

    x = []
    r = []
    g = []
    b = []
    colorModel = "RGB"
    for l in lines:
        ls = l.split()
        if not len(ls):
            continue
        if l[0] == "#":
            if ls[-1] == "HSV":
                colorModel = "HSV"
                continue
            else:
                continue
        if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
            pass
        else:
            x.append(float(ls[0]))
            r.append(float(ls[1]))
            g.append(float(ls[2]))
            b.append(float(ls[3]))
            xtemp = float(ls[4])
            rtemp = float(ls[5])
            gtemp = float(ls[6])
            btemp = float(ls[7])

    x.append(xtemp)
    r.append(rtemp)
    g.append(gtemp)
    b.append(btemp)

    nTable = len(r)
    x = np.array(x, np.float)
    r = np.array(r, np.float)
    g = np.array(g, np.float)
    b = np.array(b, np.float)
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    if colorModel == "HSV":
        for i in range(r.shape[0]):
            rr, gg, bb = colorsys.hsv_to_rgb(r[i] / 360., g[i], b[i])
            r[i] = rr
            g[i] = gg
            b[i] = bb
    if colorModel == "RGB":
        r = r / 255.
        g = g / 255.
        b = b / 255.
    xNorm = (x - x[0]) / (x[-1] - x[0])

    red = []
    blue = []
    green = []
    for i in range(len(x)):
        red.append([xNorm[i], r[i], r[i]])
        green.append([xNorm[i], g[i], g[i]])
        blue.append([xNorm[i], b[i], b[i]])
    colorDict = {"red": red, "green": green, "blue": blue}
    return (colorDict, x[0], x[-1])


class GMTColorMap(object):
    """
    Class to encapsulate linear segmented color maps with accompanying data min/max values,
    such as those found in GMT color palette table files, described here:
    http://gmt.soest.hawaii.edu/doc/5.1.0/GMT_Docs.html#b-gmt-file-formats
    """

    def __init__(self, cmap, vmin, vmax):
        """Create a GMTColorMap from an existing colormap and data min/max values.
        :param cmap:
          Matplotlib Colormap object.
        :param vmin:
          Minimum data value (associated with cmap(0.0))
        :param vmax:
          Minimum data value (associated with cmap(1.0))
        """
        self._vmin = vmin
        self._vmax = vmax
        self._cmap = cmap

    @classmethod
    def loadFromCPT(cls, cptfile, name=None):
        """Load a GMT-style color palette table (cpt) file.
        :param cptfile:
          Valid GMT cpt file.
        :param name:
          Desired Colormap name - if None, filename (without extension) will be used.
        :returns:
          GMTColorMap object.
        """
        if name is None:
            fpath, fname = os.path.split(cptfile)
            name = os.path.splitext(fname)[0]
        cdict, zmin, zmax = read_cpt(cptfile)
        cmap = LinearSegmentedColormap(name, cdict)
        return cls(cmap, zmin, zmax)

    def saveToCPT(self, filename, nsegments=None):
        """Save GMTColorMap object to a GMT cpt format.
        :param filename:
          Output file where colormap data should be written.
        :param nsegments:
          Number of desired segments (default is number in Colormap).
        """
        if nsegments is None:
            nsegments = self.cmap.N
        dv = (self.vmax - self.vmin) / nsegments
        zvalues = np.arange(self.vmin, self.vmax + dv, dv)
        f = open(filename, 'wt')
        for i in range(0, len(zvalues) - 1):
            z0 = zvalues[i]
            z1 = zvalues[i + 1]
            r0, g0, b0 = self.getRGBColor(z0)[0, 0:3]
            r1, g1, b1 = self.getRGBColor(z1)[0, 0:3]
            f.write('%.2f %i %i %i %.2f %i %i %i\n' %
                    (z0, r0, g0, b0, z1, r1, g1, b1))

        f.close()

    def getNorm(self, values):
        """Get the normalized (0-1) values associated with input physical values.
        :param values:
          Scalar or array-like set of physical values (height in meters, intensity in MMI, etc.)
        :returns:
          Normalized (0-1) value(s) corresponding to input physical values.
        """
        norm = np.array([(values - self.vmin) / (self.vmax - self.vmin)])
        return norm

    def getNormColor(self, values):
        """Get normalized (0-1) RGBA values associated with input physical values.
        :param values:
          Scalar or array-like set of physical values (height in meters, intensity in MMI, etc.)
        :returns:
          Array of normalized RGBA values corresponding to input physical values.
        """
        normvalues = self.getNorm(values)
        normcolors = self.cmap(normvalues)
        return normcolors

    def getRGBColor(self, values):
        """Get non-normalized (0-255) RGBA values associated with input physical values.
        :param values:
          Scalar or array-like set of physical values (height in meters, intensity in MMI, etc.)
        :returns:
          Array of non-normalized (0-255) RGBA values corresponding to input physical values.
        """
        rgba = self.getNormColor(values)
        rgba = np.round(rgba * 255).astype(np.int16)
        return rgba

    def getHexColor(self, values):
        rgba = self.getRGBColor(values)
        # make a list of strings
        hexcolors = []
        m, n = rgba.shape
        for i in range(0, m):
            r, g, b = rgba[i, 0:3]
            hexcolor = '#%02x%02x%02x' % (r, g, b)
            hexcolors.append(hexcolor)

        return hexcolors

    @property
    def vmin(self):
        """Get the minimum data value.
        :returns:
          vmin value.
        """
        return self._vmin

    @property
    def vmax(self):
        """Get the maximum data value.
        :returns:
          vmax value.
        """
        return self._vmax

    @property
    def cmap(self):
        """Get the matplotlib colormap.
        :returns:
          LinearSegmentedColormap derived from CPT file (or similar input).
        """
        return self._cmap
