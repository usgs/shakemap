# stdlib imports
import argparse
import inspect
import os.path
# from multiprocessing import Pool
import concurrent.futures as cf
import copy
import json

# third party
from configobj import ConfigObj
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature
from cartopy.io.shapereader import Reader

# local imports
# from mapio.gmt import GMTGrid
from impactutils.io.smcontainers import ShakeMapOutputContainer
from impactutils.mapping.city import Cities


from shakemap.utils.config import (get_config_paths,
                                   get_configspec,
                                   get_custom_validator,
                                   config_error,
                                   check_extra_values,
                                   get_data_path)
from .base import CoreModule, Contents
from shakemap.mapping.mapmaker import draw_uncertainty_map
from shakelib.utils.imt_string import oq_to_file

WATERCOLOR = '#7AA1DA'


class UncertaintymapsModule(CoreModule):
    """
    uncertaintymaps -- Generate maps of the uncertainty of the IMTs
    found in shake_result.hdf.
    """

    command_name = 'uncertaintymaps'
    targets = [
        r'products/intensity_uncertainty\.jpg',
        r'products/intensity_uncertainty\.pdf',
        r'products/pga_uncertainty\.jpg', r'products/pga_uncertainty\.pdf',
        r'products/pgv_uncertainty\.jpg', r'products/pgv_uncertainty\.pdf',
        r'products/psa_uncertainty.*p.*\.jpg',
        r'products/psa_uncertainty.*p.*\.pdf']
    dependencies = [('products/shake_result.hdf', True)]
    configs = ['products.conf']

    display_magnitude = None

    def __init__(self, eventid):
        super(UncertaintymapsModule, self).__init__(eventid)
        self.contents = Contents('Ground Motion Uncertainty Maps',
                                 'uncertainty_maps', eventid)

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --display-magnitude flag
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-m', '--display-magnitude', type=float,
                            help='Override the magnitude displayed in '
                            'map labels.')
        #
        # This line should be in any modules that overrides this
        # one. It will collect up everything after the current
        # modules options in args.rem, which should be returned
        # by this function. Note: doing parser.parse_known_args()
        # will not work as it will suck up any later modules'
        # options that are the same as this one's.
        #
        parser.add_argument('rem', nargs=argparse.REMAINDER,
                            help=argparse.SUPPRESS)
        args = parser.parse_args(arglist)
        if args.display_magnitude:
            self.display_magnitude = args.display_magnitude
        return args.rem

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
            raise NotImplementedError('uncertaintymaps module can only '
                                      'operate on gridded data, not sets of '
                                      'points')

        # get the path to the products.conf file, load the config
        config_file = os.path.join(install_path, 'config', 'products.conf')
        spec_file = get_configspec('products')
        validator = get_custom_validator()
        config = ConfigObj(config_file, configspec=spec_file)
        results = config.validate(validator)
        check_extra_values(config, self.logger)
        if not isinstance(results, bool) or not results:
            config_error(config, results)

        # create contour files
        self.logger.debug('Uncertainty mapping...')

        # get the operator setting from config
        operator = config['products']['mapping']['operator']

        # get all of the pieces needed for the uncertainty mapping functions
        layers = config['products']['mapping']['layers']
        if 'countries' in layers and layers['countries'] != '':
            countries_file = layers['countries']
        else:
            countries_file = None
        if 'states_provs' in layers and layers['states_provs'] != '':
            states_provs_file = layers['states_provs']
        else:
            states_provs_file = None
        if 'oceans' in layers and layers['oceans'] != '':
            oceans_file = layers['oceans']
        else:
            oceans_file = None
        if 'lakes' in layers and layers['lakes'] != '':
            lakes_file = layers['lakes']
        else:
            lakes_file = None

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
        padx = 5 * dx
        pady = 5 * dy
        sxmin = float(xmin) - padx
        sxmax = float(xmax) + padx
        symin = float(ymin) - pady
        symax = float(ymax) + pady

        sampledict = GeoDict.createDictFromBox(sxmin, sxmax,
                                               symin, symax,
                                               dx, dy)
        tdata = np.full([sampledict.ny, sampledict.nx], 0.0)
        topogrid = Grid2D(data=tdata, geodict=sampledict)

        model_config = container.getConfig()

        imtlist = container.getIMTs()

        textfile = os.path.join(get_data_path(), 'mapping',
                                'map_strings.' +
                                config['products']['mapping']['language'])
        text_dict = get_text_strings(textfile)
        if config['products']['mapping']['fontfamily'] != '':
            matplotlib.rcParams['font.family'] = \
                config['products']['mapping']['fontfamily']
            matplotlib.rcParams['axes.unicode_minus'] = False

        allcities = Cities.fromDefault()
        states_provs = None
        countries = None
        oceans = None
        lakes = None
        faults = None
        roads = None
        if states_provs_file is not None:
            states_provs = ShapelyFeature(
                Reader(states_provs_file).geometries(),
                ccrs.PlateCarree(), facecolor='none')
        elif 'CALLED_FROM_PYTEST' not in os.environ:
            states_provs = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_1_states_provinces_lines',
                scale='10m',
                facecolor='none')
            # The feature constructor doesn't necessarily download the
            # data, but we want it to so that multiple threads don't
            # try to do it at once when they actually access the data.
            # So below we just call the geometries() method to trigger
            # the download if necessary.
            _ = states_provs.geometries()

        if countries_file is not None:
            countries = ShapelyFeature(
                Reader(countries_file).geometries(),
                ccrs.PlateCarree(), facecolor='none')
        elif 'CALLED_FROM_PYTEST' not in os.environ:
            countries = cfeature.NaturalEarthFeature(
                category='cultural',
                name='admin_0_countries',
                scale='10m',
                facecolor='none')
            _ = countries.geometries()

        if oceans_file is not None:
            oceans = ShapelyFeature(
                Reader(oceans_file).geometries(),
                ccrs.PlateCarree(), facecolor=WATERCOLOR)
        elif 'CALLED_FROM_PYTEST' not in os.environ:
            oceans = cfeature.NaturalEarthFeature(
                category='physical',
                name='ocean',
                scale='10m',
                facecolor=WATERCOLOR)
            _ = oceans.geometries()

        if lakes_file is not None:
            lakes = ShapelyFeature(
                Reader(lakes_file).geometries(),
                ccrs.PlateCarree(), facecolor=WATERCOLOR)
        elif 'CALLED_FROM_PYTEST' not in os.environ:
            lakes = cfeature.NaturalEarthFeature(
                category='physical',
                name='lakes',
                scale='10m',
                facecolor=WATERCOLOR)
            _ = lakes.geometries()

        alist = []
        llogo = config['products']['mapping'].get('license_logo') or None
        ltext = config['products']['mapping'].get('license_text') or None
        for imtype in imtlist:
            component, imtype = imtype.split('/')
            comp = container.getComponents(imtype)[0]
            d = {'imtype': imtype,
                 'topogrid': topogrid,
                 'allcities': allcities,
                 'states_provinces': states_provs,
                 'countries': countries,
                 'oceans': oceans,
                 'lakes': lakes,
                 'roads': roads,
                 'roadcolor': layers['roadcolor'],
                 'roadwidth': layers['roadwidth'],
                 'faults': faults,
                 'faultcolor': layers['faultcolor'],
                 'faultwidth': layers['faultwidth'],
                 'datadir': datadir,
                 'operator': operator,
                 'filter_size': 0,
                 'info': info,
                 'component': comp,
                 'imtdict': container.getIMTGrids(imtype, comp),
                 'ruptdict': copy.deepcopy(container.getRuptureDict()),
                 'stationdict': container.getStationDict(),
                 'config': model_config,
                 'tdict': text_dict,
                 'display_magnitude': self.display_magnitude,
                 'pdf_dpi': config['products']['mapping']['pdf_dpi'],
                 'img_dpi': config['products']['mapping']['img_dpi'],
                 'license_logo': llogo,
                 'license_text': ltext,
                 }
            alist.append(d)
            fileimt = oq_to_file(imtype) + '_uncertainty'
            self.contents.addFile(fileimt + 'UncertaintyMap',
                                  fileimt.upper() + ' Uncertainty Map',
                                  'Map of ' + imtype + ' uncertainty.',
                                  fileimt + '.jpg', 'image/jpeg')
            self.contents.addFile(fileimt + 'UncertaintyMap',
                                  fileimt.upper() + ' Uncertainty Map',
                                  'Map of ' + imtype + ' uncertainty.',
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

    fig1 = draw_uncertainty_map(adict)

    if imtype == 'MMI':
        # save to pdf/jpeg
        pdf_file = os.path.join(adict['datadir'], 'intensity_uncertainty.pdf')
        jpg_file = os.path.join(adict['datadir'], 'intensity_uncertainty.jpg')
    else:
        fileimt = oq_to_file(imtype)
        pdf_file = os.path.join(adict['datadir'],
                                '%s_uncertainty.pdf' % (fileimt))
        jpg_file = os.path.join(adict['datadir'],
                                '%s_uncertainty.jpg' % (fileimt))

    fig1.savefig(pdf_file, bbox_inches='tight', dpi=adict['pdf_dpi'])
    fig1.savefig(jpg_file, bbox_inches='tight', dpi=adict['img_dpi'])
    plt.close(fig1)


def get_text_strings(stringfile):
    """Read the file containing the translated text strings, remove the
    comments and parse as JSON.

    Args:
        stringfile (str): Path to the map_strings.xx file specified in the
            config. The file is assumend to be UTF-8.

    Returns:
        dict: A dictionary of strings for use in writing text to the maps.
    """
    if not os.path.isfile(stringfile):
        FileNotFoundError("File %s not found" % stringfile)
    f = open(stringfile, 'rt', encoding='utf-8-sig')
    jline = ''
    for line in f:
        if line.startswith('//'):
            continue
        jline += line
    f.close()
    return json.loads(jline)
