# stdlib imports
import os.path
from collections import OrderedDict

# third party
from configobj import ConfigObj
from mapio.geodict import GeoDict
import numpy as np
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from impactutils.colors.cpalette import ColorPalette

# local imports
# from mapio.gmt import GMTGrid
from mapio.reader import read
from impactutils.io.smcontainers import ShakeMapOutputContainer

from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error)
from .base import CoreModule
from shakemap.mapping.mapmaker import (draw_intensity, draw_contour)


class MappingModule(CoreModule):
    """
    mapping -- Generate maps of the IMTs found in shake_result.hdf.
    """

    command_name = 'mapping'
    targets = [
        r'products/intensity\.jpg', r'products/intensity\.pdf',
        r'products/pga\.jpg', r'products/pga\.pdf',
        r'products/pgv\.jpg', r'products/pgv\.pdf',
        r'products/psa.*p.*\.jpg', r'products/psa.*p.*\.pdf']
    dependencies = [('products/shake_result.hdf', True)]
    configs = ['products.conf']

    # supply here a data structure with information about files that
    # can be created by this module.
    mapping_page = {
        'title': 'Ground Motion Maps',
        'slug': 'maps'
    }
    contents = OrderedDict.fromkeys(
        ['intensityMap',
         'intensityThumbnail',
         'pgaMap',
         'pgvMap',
         'psa[PERIOD]Map'])
    contents['intensityMap'] = {
        'title': 'Intensity Map',
        'caption': 'Map of macroseismic intensity.',
        'page': mapping_page,
        'formats': [{
            'filename': 'intensity.jpg',
            'type': 'image/jpeg'
        }, {
            'filename': 'intensity.pdf',
            'type': 'application/pdf'}
        ]
    }

    contents['intensityThumbnail'] = {
        'title': 'Intensity Thumbnail',
        'caption': 'Thumbnail of intensity map.',
        'page': mapping_page,
        'formats': [{
            'filename':
            'pin-thumbnail.png',
            'type': 'image/png'
        }]
    }

    contents['pgaMap'] = {
        'title': 'PGA Map',
        'caption': 'Map of peak ground acceleration (%g).',
        'page': mapping_page,
        'formats': [{
            'filename': 'pga.jpg',
            'type': 'image/jpeg'
        }, {
            'filename': 'pga.pdf',
            'type': 'image/jpeg'
        }]
    }
    contents['pgvMap'] = {
        'title': 'PGV Map',
        'caption': 'Map of peak ground velocity (cm/s).',
        'page': mapping_page,
        'formats': [{
            'filename': 'pgv.jpg',
            'type': 'image/jpeg'
        }, {
            'filename': 'pgv.pdf',
            'type': 'application/pdf'
        }]
    }
    psacap = 'Map of [FPERIOD] sec 5% damped pseudo-spectral acceleration(%g).'
    contents['psa[PERIOD]Map'] = {
        'title': 'PSA[PERIOD] Map',
        'page': mapping_page,
        'caption': psacap,
        'formats': [{
            'filename':
            'psa[0-9]p[0-9].jpg',
            'type': 'image/jpeg'
        }, {
            'filename':
            'psa[0-9]p[0-9].pdf',
            'type': 'application/pdf'
        }]
    }

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

        # get the contour filter_size setting from config
        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        config = ConfigObj(config_file, configspec=spec_file)
        results = config.validate(validator)
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        # get the filter size from the products.conf
        filter_size = config['products']['contour']['filter_size']

        # get the operator setting from config
        operator = config['products']['mapping']['operator']

        # get all of the pieces needed for the mapping functions
        layers = config['products']['mapping']['layers']
        oceanfile = layers['oceans']
        topofile = layers['topography']

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
        topogrid = read(topofile,
                        samplegeodict=sampledict,
                        resample=False)

        imtlist = container.getIMTs()
        for imtype in imtlist:
            component, imtype = imtype.split('/')
            if imtype == 'MMI':
                self.logger.debug('Drawing intensity map...')
                intensity_pdf, _ = draw_intensity(
                    container,
                    topogrid,
                    oceanfile,
                    datadir,
                    operator
                )
                self.logger.debug('Created intensity map %s' % intensity_pdf)
                self.make_pin_thumbnail(container, component, datadir)
            else:
                self.logger.debug('Drawing %s contour map...' % imtype)

                contour_pdf, contour_png = draw_contour(
                    container,
                    imtype, topogrid,
                    oceanfile,
                    datadir,
                    operator, filter_size
                )
                self.logger.debug('Created contour map %s' % contour_pdf)

        container.close()

    def make_pin_thumbnail(self, container, component, datadir):
        """Make the artsy-thumbnail for the pin on the USGS webpages.
        """
        imtdict = container.getIMTGrids("MMI", component)
        grid = imtdict['mean']
        metadata = imtdict['mean_metadata']
        rx = (np.random.rand(300) * metadata['nx']).astype(np.int)
        ry = (np.random.rand(300) * metadata['ny']).astype(np.int)
        rvals = np.arange(0, 30, 0.1)

        x_grid = np.linspace(0, metadata['nx'] - 1, 400)
        y_grid = np.linspace(0, metadata['ny'] - 1, 400)

        mx_grid, my_grid = np.meshgrid(x_grid, y_grid)

        grid = griddata(np.hstack([rx.reshape((-1, 1)), ry.reshape((-1, 1))]),
                        grid[ry, rx], (mx_grid, my_grid), method='nearest')
        grid = (grid * 10 + 0.5).astype(np.int).astype(np.float) / 10.0

        rgrid = griddata(np.hstack([rx.reshape((-1, 1)), ry.reshape((-1, 1))]),
                         rvals, (mx_grid, my_grid), method='nearest')

        mmimap = ColorPalette.fromPreset('mmi')
        plt.figure(figsize=(2.75, 2.75), dpi=96, frameon=False)
        plt.axis('off')
        plt.tight_layout()
        plt.imshow(grid, cmap=mmimap.cmap, vmin=1.5, vmax=9.5)
        plt.contour(rgrid, levels=rvals, colors='#cccccc', linewidths=0.02)
        plt.savefig(
            os.path.join(datadir, "pin-thumbnail.png"),
            dpi=96,
            bbox_inches=matplotlib.transforms.Bbox(
                [[0.47, 0.39], [2.50, 2.50]]),
            pad_inches=0)
