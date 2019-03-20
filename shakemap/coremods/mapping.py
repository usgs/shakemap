# stdlib imports
import os.path
# from multiprocessing import Pool
import concurrent.futures as cf
import copy

# third party
from configobj import ConfigObj
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
import numpy as np
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from impactutils.colors.cpalette import ColorPalette
import rasterio.features

# local imports
# from mapio.gmt import GMTGrid
from mapio.reader import read
from impactutils.io.smcontainers import ShakeMapOutputContainer

from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error)
from .base import CoreModule, Contents
from shakemap.mapping.mapmaker import (draw_intensity, draw_contour)
from shakelib.utils.imt_string import oq_to_file


class MappingModule(CoreModule):
    """
    mapping -- Generate maps of the IMTs found in shake_result.hdf.
    """

    command_name = 'mapping'
    targets = [
        r'products/intensity\.jpg', r'products/intensity\.pdf',
        r'products/mmi_legend\.pdf',
        r'products/pga\.jpg', r'products/pga\.pdf',
        r'products/pgv\.jpg', r'products/pgv\.pdf',
        r'products/psa.*p.*\.jpg', r'products/psa.*p.*\.pdf']
    dependencies = [('products/shake_result.hdf', True)]
    configs = ['products.conf']

    def __init__(self, eventid):
        super(MappingModule, self).__init__(eventid)
        self.contents = Contents('Ground Motion Maps', 'maps', eventid)

    def execute(self):
        """
        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_result HDF file does not
                exist.
        """
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current', 'products')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)
        datafile = os.path.join(datadir, 'shake_result.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # Open the ShakeMapOutputContainer and extract the data
        container = ShakeMapOutputContainer.load(datafile)
        if container.getDataType() != 'grid':
            raise NotImplementedError('mapping module can only operate on '
                                      'gridded data, not sets of points')

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        config = ConfigObj(config_file, configspec=spec_file)
        results = config.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        # create contour files
        self.logger.debug('Mapping...')

        # get the filter size from the products.conf
        filter_size = config['products']['contour']['filter_size']

        # get the operator setting from config
        operator = config['products']['mapping']['operator']

        # get all of the pieces needed for the mapping functions
        layers = config['products']['mapping']['layers']
        oceanfile = layers['oceans']
        if 'topography' in layers and layers['topography'] != '':
            topofile = layers['topography']
        else:
            topofile = None

        # Get the number of parallel workers
        max_workers = config['products']['mapping']['max_workers']

        # Reading HDF5 files currently takes a long time, due to poor
        # programming in MapIO.  To save us some time until that issue is
        # resolved, we'll coarsely subset the topo grid once here and pass
        # it into both mapping functions
        # get the bounds of the map
        info = container.getMetadata()
        xmin = info['output']['map_information']['min']['longitude']
        xmax = info['output']['map_information']['max']['longitude']
        ymin = info['output']['map_information']['min']['latitude']
        ymax = info['output']['map_information']['max']['latitude']
        dy = float(info['output']['map_information']
                   ['grid_spacing']['latitude'])
        dx = float(info['output']['map_information']
                   ['grid_spacing']['longitude'])
        padx = 5*dx
        pady = 5*dy
        sxmin = float(xmin) - padx
        sxmax = float(xmax) + padx
        symin = float(ymin) - pady
        symax = float(ymax) + pady

        sampledict = GeoDict.createDictFromBox(sxmin, sxmax,
                                               symin, symax,
                                               dx, dy)
        if topofile:
            topogrid = read(topofile,
                            samplegeodict=sampledict,
                            resample=False)
        else:
            tdata = np.full([sampledict.ny, sampledict.nx], 0.0)
            topogrid = Grid2D(data=tdata, geodict=sampledict)

        model_config = container.getConfig()

        imtlist = container.getIMTs()

        alist = []
        for imtype in imtlist:
            component, imtype = imtype.split('/')
            comp = container.getComponents(imtype)[0]
            d = {'imtype': imtype,
                 'topogrid': topogrid,
                 'oceanfile': oceanfile,
                 'datadir': datadir,
                 'operator': operator,
                 'filter_size': filter_size,
                 'info': info,
                 'component': comp,
                 'imtdict': container.getIMTGrids(imtype, comp),
                 'ruptdict': copy.deepcopy(container.getRuptureDict()),
                 'stationdict': container.getStationDict(),
                 'config': model_config
                 }
            alist.append(d)
            if imtype == 'MMI':
                g = copy.deepcopy(d)
                g['imtype'] = 'thumbnail'
                alist.append(g)
                self.contents.addFile('intensityMap', 'Intensity Map',
                                      'Map of macroseismic intensity.',
                                      'intensity.jpg', 'image/jpeg')
                self.contents.addFile('intensityMap', 'Intensity Map',
                                      'Map of macroseismic intensity.',
                                      'intensity.pdf', 'application/pdf')
                self.contents.addFile('intensityThumbnail',
                                      'Intensity Thumbnail',
                                      'Thumbnail of intensity map.',
                                      'pin-thumbnail.png', 'image/png')
            else:
                fileimt = oq_to_file(imtype)
                self.contents.addFile(fileimt + 'Map',
                                      fileimt.upper() + ' Map',
                                      'Map of ' + imtype + '.',
                                      fileimt + '.jpg', 'image/jpeg')
                self.contents.addFile(fileimt + 'Map',
                                      fileimt.upper() + ' Map',
                                      'Map of ' + imtype + '.',
                                      fileimt + '.pdf', 'application/pdf')

        if max_workers > 0:
            with cf.ProcessPoolExecutor(max_workers=max_workers) as ex:
                results = ex.map(make_map, alist)
                list(results)
        else:
            for adict in alist:
                make_map(adict)

        container.close()


def make_map(adict):

    imtype = adict['imtype']
    if imtype == 'MMI':
        fig1, fig2 = draw_intensity(adict)
        # save to pdf/jpeg
        pdf_file = os.path.join(adict['datadir'], 'intensity.pdf')
        jpg_file = os.path.join(adict['datadir'], 'intensity.jpg')
        fig1.savefig(pdf_file, bbox_inches='tight')
        fig1.savefig(jpg_file, bbox_inches='tight')

        # save the legend file
        legend_file = os.path.join(adict['datadir'], 'mmi_legend.png')
        fig2.savefig(legend_file, bbox_inches='tight')
    elif imtype == 'thumbnail':
        make_pin_thumbnail(adict)
    else:
        fig1 = draw_contour(adict)
        fileimt = oq_to_file(imtype)
        pdf_file = os.path.join(adict['datadir'], '%s.pdf' % (fileimt))
        jpg_file = os.path.join(adict['datadir'], '%s.jpg' % (fileimt))
        fig1.savefig(pdf_file, bbox_inches='tight')
        fig1.savefig(jpg_file, bbox_inches='tight')


def make_pin_thumbnail(adict):
    """Make the artsy-thumbnail for the pin on the USGS webpages.
    """
    imtdict = adict['imtdict']
    grid = imtdict['mean']
    metadata = imtdict['mean_metadata']
    num_pixels = 300
    randx = np.random.rand(num_pixels)
    randy = np.random.rand(num_pixels)
    rx = (randx * metadata['nx']).astype(np.int)
    ry = (randy * metadata['ny']).astype(np.int)
    rvals = np.arange(num_pixels)

    x_grid = np.arange(400)
    y_grid = np.arange(400)

    mx_grid, my_grid = np.meshgrid(x_grid, y_grid)

    grid = griddata(np.hstack([randx.reshape((-1, 1)) * 400,
                               randy.reshape((-1, 1)) * 400]),
                    grid[ry, rx], (mx_grid, my_grid), method='nearest')
    grid = (grid * 10 + 0.5).astype(np.int).astype(np.float) / 10.0

    rgrid = griddata(np.hstack([randx.reshape((-1, 1)) * 400,
                                randy.reshape((-1, 1)) * 400]),
                     rvals, (mx_grid, my_grid), method='nearest')
    irgrid = rgrid.astype(np.int32)
    mypols = [p[0]['coordinates'] for p in rasterio.features.shapes(irgrid)]

    mmimap = ColorPalette.fromPreset('mmi')
    plt.figure(figsize=(2.75, 2.75), dpi=96, frameon=False)
    plt.axis('off')
    plt.tight_layout()
    plt.imshow(grid, cmap=mmimap.cmap, vmin=1.5, vmax=9.5)
    for pol in mypols:
        mycoords = list(zip(*pol[0]))
        plt.plot(mycoords[0], mycoords[1], color='#cccccc', linewidth=0.2)
    plt.savefig(
        os.path.join(adict['datadir'], "pin-thumbnail.png"),
        dpi=96,
        bbox_inches=matplotlib.transforms.Bbox(
            [[0.47, 0.39], [2.50, 2.50]]),
        pad_inches=0)
