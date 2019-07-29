"""
Process a ShakeMap, based on the configuration and data found in
shake_data.hdf, and produce output in shake_result.hdf.
"""
import os
import argparse
import inspect
import os.path
import time as time
import copy
from time import gmtime, strftime
import shutil
from datetime import date
import json

import numpy as np
import numpy.ma as ma
from openquake.hazardlib import imt
import openquake.hazardlib.const as oqconst
import fiona
import cartopy.io.shapereader as shpreader
from shapely.geometry import shape

import concurrent.futures as cf

# local imports
from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D
from .base import CoreModule, Contents
from shakelib.rupture.point_rupture import PointRupture
from shakelib.sites import Sites
from shakelib.distance import (Distance,
                               get_distance,
                               get_distance_measures)
from shakelib.multigmpe import MultiGMPE
from shakelib.virtualipe import VirtualIPE
from shakelib.utils.utils import get_extent
from shakelib.utils.imt_string import oq_to_file
from shakelib.utils.containers import ShakeMapInputContainer
from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakelib.rupture import constants

from shakemap.utils.config import get_config_paths
from shakemap.utils.utils import get_object_from_config
from shakemap._version import get_versions
from shakemap.utils.generic_amp import get_generic_amp_factors
from shakemap.c.clib import (make_sigma_matrix,
                             geodetic_distance_fast,
                             make_sd_array)
# from shakemap.utils.exception import TerminateShakeMap

from shakelib.directivity.rowshandel2013 import Rowshandel2013

#
# TODO: Some constants; these should maybe be in a configuration
# or in a constants module or be GMICE-specific.
#
# default_mmi_stddev: the standard deviation of MMI values if it
#                     is not specified in the input
# min_mmi_convert: the minimum MMI to convert to PGM -- low
#                  intensities don't convert very accurately
# default_stddev_inter: This is a stand-in for tau when the gmpe set
#                       doesn't provide it. It is an educated guess based
#                       on the NGA-west, Akkar et al, and BC Hydro gmpes.
#                       It's not perfect, but probably isn't too far off.
#                       It is only used when the GMPE(s) don't provide a
#                       breakdown of the uncertainty terms. When used,
#                       this value is multiplied by the total standard
#                       deviation to get tau. The square of tau is then
#                       subtracted from the squared total stddev and the
#                       square root of the result is then used as the
#                       within-event stddev (phi).
#
SM_CONSTS = {
    'default_mmi_stddev': 0.3,
    'min_mmi_convert': 4.0,
    'default_stddev_inter': 0.55,
    'default_stddev_inter_mmi': 0.55
}


class DataFrame:
    def __init__(self):
        df = None  # noqa
        imts = None  # noqa
        sx = None  # noqa
        dx = None  # noqa


class ModelModule(CoreModule):
    """
    model -- Interpolate ground motions to a grid or list of locations.
    """

    command_name = 'model'
    targets = [r'products/shake_result\.hdf']
    dependencies = [('shake_data.hdf', True)]

    no_seismic = False
    no_macroseismic = False
    no_rupture = False
    use_simulations = False

    rock_vs30 = 760.0
    soil_vs30 = 180.0

    def __init__(self, eventid):
        super(ModelModule, self).__init__(eventid)
        self.contents = Contents(None, None, eventid)

    def parseArgs(self, arglist):
        """
        Set up the object to accept the --no_seismic, --no_macroseismic,
        and --no_rupture flags.
        """
        parser = argparse.ArgumentParser(
            prog=self.__class__.command_name,
            description=inspect.getdoc(self.__class__))
        parser.add_argument('-s', '--no_seismic', action='store_true',
                            help='Exclude instrumental seismic data from '
                            'the processing, ignoring any that may exist in '
                            'the input directory.')
        parser.add_argument('-m', '--no_macroseismic', action='store_true',
                            help='Exclude macroseismic data from the '
                            'processing, ignoring any that may exist in the '
                            'input directory.')
        parser.add_argument('-r', '--no_rupture', action='store_true',
                            help='Exclude a rupture model from the '
                            'processing, ignoring any that may exist in the '
                            'input directory.')
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
        if args.no_seismic:
            self.no_seismic = True
        if args.no_macroseismic:
            self.no_macroseismic = True
        if args.no_rupture:
            self.no_rupture = True
        return args.rem

    def execute(self):
        """
        Interpolate ground motions to a grid or list of locations.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_data HDF file does not exist.
        """
        self.logger.debug('Starting model...')
        # ---------------------------------------------------------------------
        # Make the input container and extract the config
        # ---------------------------------------------------------------------
        self._setInputContainer()
        self.config = self.ic.getConfig()

        self.sim_imt_paths = [x for x in self.ic.getArrays()
                              if 'simulations' in x]
        if len(self.sim_imt_paths):
            self.use_simulations = True

        # ---------------------------------------------------------------------
        # Clear away results from previous runs
        # ---------------------------------------------------------------------
        self._clearProducts()

        # ---------------------------------------------------------------------
        # Retrieve a bunch of config options and set them as attributes
        # ---------------------------------------------------------------------
        self._setConfigOptions()

        # ---------------------------------------------------------------------
        # Instantiate the gmpe, gmice, and ipe
        # Here we make a placeholder gmpe so that we can make the
        # rupture and distance contexts; later we'll make the
        # IMT-specific gmpes
        # ---------------------------------------------------------------------
        self.default_gmpe = MultiGMPE.from_config(self.config)

        self.gmice = get_object_from_config('gmice', 'modeling', self.config)

        if self.config['ipe_modules'][self.config['modeling']['ipe']][0] == \
                'VirtualIPE':
            pgv_imt = imt.from_string('PGV')
            ipe_gmpe = MultiGMPE.from_config(self.config, filter_imt=pgv_imt)
            self.ipe = VirtualIPE.fromFuncs(ipe_gmpe, self.gmice)
        else:
            ipe = get_object_from_config('ipe', 'modeling', self.config)
            if 'vs30' not in ipe.REQUIRES_SITES_PARAMETERS:
                ipe.REQUIRES_SITES_PARAMETERS.add('vs30')
            ipe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = \
                oqconst.IMC.GREATER_OF_TWO_HORIZONTAL
            self.ipe = MultiGMPE.from_list([ipe], [1.0])

        ipe_sd_types = self.ipe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if len(ipe_sd_types) == 1:
            self.ipe_total_sd_only = True
            self.ipe_stddev_types = [oqconst.StdDev.TOTAL]
        else:
            self.ipe_total_sd_only = False
            self.ipe_stddev_types = [oqconst.StdDev.TOTAL,
                                     oqconst.StdDev.INTER_EVENT,
                                     oqconst.StdDev.INTRA_EVENT]

        # ---------------------------------------------------------------------
        # Get the rupture object and rupture context
        # ---------------------------------------------------------------------
        self.rupture_obj = self.ic.getRuptureObject()
        # If the --no_rupture flag is used, switch to a PointRupture
        if self.no_rupture:
            self.rupture_obj = PointRupture(self.rupture_obj._origin)
        if 'mechanism' in self.config['modeling']:
            self.rupture_obj._origin.setMechanism(
                mech=self.config['modeling']['mechanism'])
        self.rx = self.rupture_obj.getRuptureContext([self.default_gmpe])
        # TODO: figure out how to not have to do this
        if self.rx.rake is None:
            self.rx.rake = 0

        #
        # Set up the coordinates for the attenuation curves
        #
        repi = np.logspace(-1, 3, 200)
        pt = self.rupture_obj._origin.getHypo()
        self.atten_coords = {
            'lons': np.full_like(repi, pt.x),
            'lats': np.array([pt.y + x / 111.0 for x in repi])
        }
        self.point_source = PointRupture(self.rupture_obj._origin)

        # ---------------------------------------------------------------------
        # The output locations: either a grid or a list of points
        # ---------------------------------------------------------------------
        self.logger.debug('Setting output params...')
        self._setOutputParams()

        landmask = self._getLandMask()
        # We used to do this, but we've decided not to. Leaving the code
        # in place in case we change our minds.
        # if landmask is not None and np.all(landmask):
        #     raise TerminateShakeMap("Mapped area is entirely water")

        # ---------------------------------------------------------------------
        # If the gmpe doesn't break down its stardard deviation into
        # within- and between-event terms, we need to handle things
        # somewhat differently.
        # ---------------------------------------------------------------------
        gmpe_sd_types = self.default_gmpe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if len(gmpe_sd_types) == 1:
            self.gmpe_total_sd_only = True
            self.gmpe_stddev_types = [oqconst.StdDev.TOTAL]
        else:
            self.gmpe_total_sd_only = False
            self.gmpe_stddev_types = [oqconst.StdDev.TOTAL,
                                      oqconst.StdDev.INTER_EVENT,
                                      oqconst.StdDev.INTRA_EVENT]

        # ---------------------------------------------------------------------
        # Are we going to include directivity?
        # ---------------------------------------------------------------------
        # Config option?
        dir_conf = self.config['modeling']['directivity']

        # Is the rupture not a point source?
        rup_check = not isinstance(self.rupture_obj, PointRupture)

        if dir_conf and rup_check:
            self.do_directivity = True
            # The following attribute will be used to store a list of tuples,
            # where each tuple will contain the 1) result of the directivity
            # model (for the periods defined by Rowshandel2013) and 2) the
            # associated distance context. The distance context is needed
            # within the _gmas function for figuring out which of the results
            # should be used when combining it with the GMPE result. We store
            # the pre-defined period first and interpolate later because there
            # is some optimization to doing it this way (some of the
            # calculation is period independent).
            self.dir_results = []
            # But if we want to save the results that were actually used for
            # each IMT, so we use a dictionary. This uses keys that are
            # the same ass self.outgrid.
            self.dir_output = {}
        else:
            self.do_directivity = False

        # ---------------------------------------------------------------------
        # Station data: Create DataFrame(s) with the input data:
        # df1 for instrumented data
        # df2 for non-instrumented data
        # ---------------------------------------------------------------------
        self.logger.debug('Setting data frames...')
        self._setDataFrames()

        # ---------------------------------------------------------------------
        # Add the predictions, etc. to the data frames
        # ---------------------------------------------------------------------
        self.logger.debug('Populating data frames...')
        self._populateDataFrames()

        # ---------------------------------------------------------------------
        # Try to make all the derived IMTs possible from MMI (if we have MMI)
        # ---------------------------------------------------------------------
        self._deriveIMTsFromMMI()
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        self._deriveMMIFromIMTs()

        self.logger.debug('Getting combined IMTs')

        # ---------------------------------------------------------------------
        # Get the combined set of input and output IMTs, their periods,
        # and an index dictionary, then make the cross-correlation function
        # ---------------------------------------------------------------------
        if self.use_simulations:
            #
            # Ignore what is in the configuration and make maps only for the
            # IMTs that are in the set of simulations (and MMI).
            #
            self.combined_imt_set = set(
                [x.split('/')[-1] for x in self.sim_imt_paths])
            self.sim_df = {}
            for imtstr in self.combined_imt_set:
                dset, _ = self.ic.getArray(['simulations'], imtstr)
                self.sim_df[imtstr] = dset
            self.combined_imt_set |= set(['MMI'])
        else:
            self.combined_imt_set = self.imt_out_set.copy()
            for ndf in self.dataframes:
                self.combined_imt_set |= getattr(self, ndf).imts

        self.imt_per, self.imt_per_ix = _get_period_arrays(
            self.combined_imt_set)
        self.ccf = get_object_from_config('ccf', 'modeling',
                                          self.config, self.imt_per)

        self.logger.debug('Doing bias')

        # ---------------------------------------------------------------------
        # Do the bias for all of the input and output IMTs. Hold on
        # to some of the products that will be used for the interpolation.
        # The "raw" values are the stddevs that have not been inflated by
        # the additional sigma (if any) of the point-source to finite
        # rupture approximation.
        # ---------------------------------------------------------------------
        self.nominal_bias = {}
        self.bias_num = {}  # The numerator term of the bias
        self.bias_den = {}  # The denominator of the bias (w/o the tau term)
        self.psd_raw = {}   # raw phi (intra-event stddev) of the output points
        self.psd = {}       # phi (intra-event stddev) of the output points
        self.tsd = {}       # tau (inter-event stddev) of the output points

        #
        # These are arrays (keyed by IMT) of the station data that will be
        # used to compute the bias and do the interpolation, they are filled
        # in the _fillDataArrays method
        #
        self.sta_period_ix = {}
        self.sta_lons_rad = {}
        self.sta_lats_rad = {}
        self.sta_resids = {}
        self.sta_phi = {}
        self.sta_sig_extra = {}
        self.sta_sig_total = {}
        self.sta_tau = {}
        self.sta_imtstr = {}
        self.sta_rrups = {}

        self._fillDataArrays()

        self._computeBias()

        self._computeDirectivityPredictionLocations()

        # ---------------------------------------------------------------------
        # Now do the MVN with the intra-event residuals
        # ---------------------------------------------------------------------
        self.outgrid = {}   # Holds the interpolated output arrays keyed by IMT
        self.outsd = {}     # Holds the standard deviation arrays keyed by IMT

        #
        # Places to put the results for the attenuation plots
        #
        self.atten_rock_mean = {}
        self.atten_soil_mean = {}
        self.atten_rock_sd = {}
        self.atten_soil_sd = {}

        self.logger.debug('Doing MVN...')
        if self.max_workers > 0:
            with cf.ThreadPoolExecutor(max_workers=self.max_workers) as ex:
                results = ex.map(self._computeMVN, self.imt_out_set)
                list(results)  # Check threads for possible exceptions, etc.
        else:
            for imt_str in self.imt_out_set:
                self._computeMVN(imt_str)

        # ---------------------------------------------------------------------
        # Output the data and metadata
        # ---------------------------------------------------------------------
        product_path = os.path.join(self.datadir, 'products')
        if not os.path.isdir(product_path):
            os.mkdir(product_path)
        oc = ShakeMapOutputContainer.create(os.path.join(
            product_path, 'shake_result.hdf'))

        # ---------------------------------------------------------------------
        # Might as well stick the whole config in the result
        # ---------------------------------------------------------------------
        oc.setConfig(self.config)

        #
        # We're going to need masked arrays of the output grids later, so
        # make them now.
        #
        moutgrid = self._getMaskedGrids(landmask)

        #
        # Get the info dictionary that will become info.json, and
        # store it in the output container
        #
        info = self._getInfo(moutgrid)
        oc.setMetadata(info)

        # ---------------------------------------------------------------------
        # Add the rupture JSON as a text string
        # ---------------------------------------------------------------------
        oc.setRuptureDict(self.rupture_obj._geojson)

        # ---------------------------------------------------------------------
        # Fill the station dictionary for stationlist.json and add it to
        # the output container
        # ---------------------------------------------------------------------
        sjdict = self._fillStationJSON()
        oc.setStationDict(sjdict)

        # ---------------------------------------------------------------------
        # Add the output grids or points to the output; include some
        # metadata.
        # ---------------------------------------------------------------------
        if self.do_grid:
            self._storeGriddedData(oc)
        else:
            self._storePointData(oc)

        self._storeAttenuationData(oc)

        oc.close()
        self.ic.close()

        self.contents.addFile('shakemapHDF',
                              'Comprehensive ShakeMap HDF Data File',
                              'HDF file containing all ShakeMap results.',
                              'shake_result.hdf', 'application/x-bag')
    # -------------------------------------------------------------------------
    # End execute()
    # -------------------------------------------------------------------------

    def _setInputContainer(self):
        """
        Open the input container and set
        the event's current data directory.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_data HDF file does not exist.
        """
        #
        # Find the shake_data.hdf file
        #
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        datafile = os.path.join(datadir, 'shake_data.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)
        self.datadir = datadir
        self.ic = ShakeMapInputContainer.load(datafile)

    def _clearProducts(self):
        """
        Function to delete an event's products directory if it exists.

        Returns:
            nothing
        """
        products_path = os.path.join(self.datadir, 'products')
        if os.path.isdir(products_path):
            shutil.rmtree(products_path, ignore_errors=True)
        pdl_path = os.path.join(self.datadir, 'pdl')
        if os.path.isdir(pdl_path):
            shutil.rmtree(pdl_path, ignore_errors=True)

    def _setConfigOptions(self):
        """
        Pull various useful configuration options out of the config
        dictionary.

        Returns:
            nothing
        """
        # ---------------------------------------------------------------------
        # Processing parameters
        # ---------------------------------------------------------------------
        self.max_workers = self.config['system']['max_workers']

        # ---------------------------------------------------------------------
        # Do we apply the generic amplification factors?
        # ---------------------------------------------------------------------
        self.apply_gafs = self.config['modeling']['apply_generic_amp_factors']

        # ---------------------------------------------------------------------
        # Bias parameters
        # ---------------------------------------------------------------------
        self.do_bias = self.config['modeling']['bias']['do_bias']
        self.bias_max_range = self.config['modeling']['bias']['max_range']
        self.bias_max_mag = self.config['modeling']['bias']['max_mag']
        self.bias_max_dsigma = \
            self.config['modeling']['bias']['max_delta_sigma']

        # ---------------------------------------------------------------------
        # Outlier parameters
        # ---------------------------------------------------------------------
        self.do_outliers = self.config['data']['outlier']['do_outliers']
        self.outlier_deviation_level = \
            self.config['data']['outlier']['max_deviation']
        self.outlier_max_mag = self.config['data']['outlier']['max_mag']
        self.outlier_valid_stations = \
            self.config['data']['outlier']['valid_stations']

        # ---------------------------------------------------------------------
        # These are the IMTs we want to make
        # ---------------------------------------------------------------------
        self.imt_out_set = set(self.config['interp']['imt_list'])

        # ---------------------------------------------------------------------
        # The x and y resolution of the output grid
        # ---------------------------------------------------------------------
        self.smdx = self.config['interp']['prediction_location']['xres']
        self.smdy = self.config['interp']['prediction_location']['yres']
        self.nmax = self.config['interp']['prediction_location']['nmax']

        # ---------------------------------------------------------------------
        # Get the Vs30 file name
        # ---------------------------------------------------------------------
        self.vs30default = self.config['data']['vs30default']
        self.vs30_file = self.config['data']['vs30file']
        if not self.vs30_file:
            self.vs30_file = None

    def _setOutputParams(self):
        """
        Set variables dealing with the output grid or points

        Returns:
            nothing
        """
        if self.use_simulations:
            self.do_grid = True
            imt_grp = self.sim_imt_paths[0]
            groups = imt_grp.split('/')
            myimt = groups[-1]
            del groups[-1]
            data, geodict = self.ic.getArray(groups, myimt)
            self.W = geodict['xmin']
            self.E = geodict['xmax']
            self.S = geodict['ymin']
            self.N = geodict['ymax']
            self.smdx = geodict['dx']
            self.smdy = geodict['dy']

            self.sites_obj_out = Sites.fromBounds(
                self.W, self.E, self.S,
                self.N, self.smdx, self.smdy,
                defaultVs30=self.vs30default,
                vs30File=self.vs30_file,
                padding=True, resample=True
            )
            self.smnx, self.smny = self.sites_obj_out.getNxNy()
            self.sx_out = self.sites_obj_out.getSitesContext()
            lons, lats = np.meshgrid(self.sx_out.lons,
                                     self.sx_out.lats)
            self.sx_out.lons = lons.copy()
            self.sx_out.lats = lats.copy()
            self.lons = lons.ravel()
            self.lats = lats.ravel()
            self.depths = np.zeros_like(lats)
            dist_obj_out = Distance.fromSites(self.default_gmpe,
                                              self.sites_obj_out,
                                              self.rupture_obj)
        elif (self.config['interp']['prediction_location']['file'] and
              self.config['interp']['prediction_location']['file'] != 'None'):
            #
            # FILE: Open the file and get the output points
            #
            self.do_grid = False
            in_sites = np.genfromtxt(
                self.config['interp']['prediction_location']['file'],
                autostrip=True, unpack=True,
                dtype=[np.double, np.double, np.double, '<U80'])
            if np.size(in_sites) == 0:
                self.logger.info('Points file is empty; nothing to do')
                return
            elif np.size(in_sites) == 1:
                lons, lats, vs30, idents = in_sites.item()
                self.idents = [idents]
            else:
                lons, lats, vs30, self.idents = zip(*in_sites)
            self.lons = np.array(lons).reshape(1, -1)
            self.lats = np.array(lats).reshape(1, -1)
            self.vs30 = np.array(vs30).reshape(1, -1)
            self.depths = np.zeros_like(self.lats)
            self.W = np.min(self.lons)
            self.E = np.max(self.lons)
            self.S = np.min(self.lats)
            self.N = np.max(self.lats)
            self.smnx = np.size(self.lons)
            self.smny = 1
            dist_obj_out = Distance(
                self.default_gmpe, self.lons, self.lats,
                self.depths, self.rupture_obj
            )

            self.sites_obj_out = Sites.fromBounds(
                self.W, self.E, self.S,
                self.N, self.smdx, self.smdy,
                defaultVs30=self.vs30default,
                vs30File=self.vs30_file,
                padding=True, resample=True
            )

            self.sx_out = self.sites_obj_out.getSitesContext(
                {'lats': self.lats,
                 'lons': self.lons}
            )
            # Replace the Vs30 from the grid (or default) with the Vs30
            # provided with the site list.
            if np.any(self.vs30 > 0):
                self.sx_out.vs30 = self.vs30
        else:
            #
            # GRID: Figure out the grid parameters and get output points
            #
            self.do_grid = True

            if self.config['interp']['prediction_location']['extent']:
                self.W, self.S, self.E, self.N = \
                    self.config['interp']['prediction_location']['extent']
            else:
                self.W, self.E, self.S, self.N = get_extent(
                    self.rupture_obj,
                    config=self.config
                )

            # Adjust resolution to be under nmax
            self._adjustResolution()

            self.sites_obj_out = Sites.fromBounds(
                self.W, self.E, self.S,
                self.N, self.smdx, self.smdy,
                defaultVs30=self.vs30default,
                vs30File=self.vs30_file,
                padding=True, resample=True
            )
            self.smnx, self.smny = self.sites_obj_out.getNxNy()
            self.sx_out = self.sites_obj_out.getSitesContext()
            lons, lats = np.meshgrid(self.sx_out.lons,
                                     self.sx_out.lats)
            self.sx_out.lons = lons.copy()
            self.sx_out.lats = lats.copy()
            self.lons = lons.ravel()
            self.lats = lats.ravel()
            self.depths = np.zeros_like(lats)
            dist_obj_out = Distance.fromSites(self.default_gmpe,
                                              self.sites_obj_out,
                                              self.rupture_obj)

        #
        # TODO: This will break if the IPE needs distance measures
        # that the GMPE doesn't; should make this a union of the
        # requirements of both
        #
        self.dx_out = dist_obj_out.getDistanceContext()
        #
        # Set up the sites and distance contexts for the attenuation curves
        #
        self.atten_sx_rock = self.sites_obj_out.getSitesContext(
            self.atten_coords, rock_vs30=self.rock_vs30)
        self.atten_sx_soil = self.sites_obj_out.getSitesContext(
            self.atten_coords, rock_vs30=self.soil_vs30)
        self.atten_dx = Distance(
            self.default_gmpe,
            self.atten_coords['lons'],
            self.atten_coords['lats'],
            np.zeros_like(self.atten_coords['lons']),
            rupture=self.point_source).getDistanceContext()

        self.lons_out_rad = np.radians(self.lons)
        self.lats_out_rad = np.radians(self.lats)
        self.flip_lons = False
        if self.W > 0 and self.E < 0:
            self.flip_lons = True
            self.lons_out_rad[self.lons_out_rad < 0] += 2 * np.pi

    def _setDataFrames(self):
        """
        Extract the StationList object from the input container and
        fill the DataFrame class and keep a list of dataframes.

            - df1 holds the instrumented data (PGA, PGV, SA)
            - df2 holds the non-instrumented data (MMI)
        """
        self.dataframes = []
        try:
            self.stations = self.ic.getStationList()
        except AttributeError:
            return
        if self.stations is None:
            return
        for dfid, val in (('df1', True), ('df2', False)):
            if dfid == 'df1' and self.no_seismic:
                continue
            if dfid == 'df2' and self.no_macroseismic:
                continue
            sdf, imts = self.stations.getStationDictionary(instrumented=val)
            if sdf is not None:
                df = DataFrame()
                df.df = sdf
                df.imts = imts
                setattr(self, dfid, df)
                self.dataframes.append(dfid)

        # Flag the stations in the bad stations list from the config
        if not hasattr(self, 'df1'):
            return
        evdt = date(self.rupture_obj._origin.time.year,
                    self.rupture_obj._origin.time.month,
                    self.rupture_obj._origin.time.day)
        nostart = date(1970, 1, 1)
        self.df1.df['flagged'] = np.full_like(self.df1.df['lon'], 0,
                                              dtype=np.bool)
        if 'bad_stations' not in self.config['data']:
            return
        for sid, dates in self.config['data']['bad_stations'].items():
            ondate, offdate = dates.split(':')
            year, month, day = map(int, ondate.split('-'))
            ondt = date(year, month, day)
            if offdate:
                year, month, day = map(int, offdate.split('-'))
                offdt = date(year, month, day)
            else:
                offdt = None
            bad = False
            if (ondt == nostart or ondt <= evdt) and \
                    (offdt is None or offdt >= evdt):
                bad = True
            if bad:
                self.df1.df['flagged'] |= (self.df1.df['id'] == sid)

    def _populateDataFrames(self):
        """
        Make the sites and distance contexts for each dataframe then
        compute the predictions for the IMTs in that dataframe.
        """
        for dfid in self.dataframes:
            dfn = getattr(self, dfid)
            df = dfn.df
            # -----------------------------------------------------------------
            # Get the sites and distance contexts
            # -----------------------------------------------------------------
            df['depth'] = np.zeros_like(df['lon'])
            lldict = {
                'lons': df['lon'],
                'lats': df['lat']
            }
            dfn.sx = self.sites_obj_out.getSitesContext(lldict)
            dfn.sx_rock = copy.deepcopy(dfn.sx)
            dfn.sx_rock.vs30 = np.full_like(dfn.sx.vs30, self.rock_vs30)
            dfn.sx_soil = copy.deepcopy(dfn.sx)
            dfn.sx_soil.vs30 = np.full_like(dfn.sx.vs30, self.soil_vs30)
            dist_obj = Distance(
                self.default_gmpe, df['lon'], df['lat'],
                df['depth'], self.rupture_obj
            )
            dfn.dx = dist_obj.getDistanceContext()

            # -----------------------------------------------------------------
            # Are we doing directivity?
            # -----------------------------------------------------------------
            if self.do_directivity is True:
                self.logger.info('Directivity for %s...' % dfid)
                time1 = time.time()
                dir_df = Rowshandel2013(
                    self.rupture_obj._origin, self.rupture_obj,
                    df['lat'].reshape((1, -1)),
                    df['lon'].reshape((1, -1)),
                    df['depth'].reshape((1, -1)),
                    dx=1.0,
                    T=Rowshandel2013.getPeriods(),
                    a_weight=0.5,
                    mtype=1
                )
                self.dir_results.append((dir_df, dfn.dx))
                directivity_time = time.time() - time1
                self.logger.debug(
                    'Directivity %s evaluation time: %f sec'
                    % (dfid, directivity_time))

            # -----------------------------------------------------------------
            # Do the predictions and other bookkeeping for each IMT
            # -----------------------------------------------------------------
            for imtstr in dfn.imts:
                oqimt = imt.from_string(imtstr)
                gmpe = None
                not_supported = False
                if imtstr != 'MMI':
                    try:
                        gmpe = MultiGMPE.from_config(self.config,
                                                     filter_imt=oqimt)
                    except KeyError:
                        self.logger.warn(
                            "Input IMT %s not supported by GMPE: ignoring" %
                            imtstr)
                        not_supported = True
                if not_supported:
                    pmean = np.full_like(df[imtstr], np.nan)
                    pmean_rock = np.full_like(df[imtstr], np.nan)
                    pmean_soil = np.full_like(df[imtstr], np.nan)
                    pstddev = [None] * 3
                    pstddev[0] = np.full_like(df[imtstr], np.nan)
                    pstddev[1] = np.full_like(df[imtstr], np.nan)
                    pstddev[2] = np.full_like(df[imtstr], np.nan)
                    pstddev_rock = [None] * 1
                    pstddev_soil = [None] * 1
                    pstddev_rock[0] = np.full_like(df[imtstr], np.nan)
                    pstddev_soil[0] = np.full_like(df[imtstr], np.nan)
                else:
                    pmean, pstddev = self._gmas(gmpe, dfn.sx, dfn.dx, oqimt,
                                                self.apply_gafs)
                    pmean_rock, pstddev_rock = self._gmas(gmpe, dfn.sx_rock,
                                                          dfn.dx, oqimt,
                                                          self.apply_gafs)
                    pmean_soil, pstddev_soil = self._gmas(gmpe, dfn.sx_soil,
                                                          dfn.dx, oqimt,
                                                          self.apply_gafs)
                df[imtstr + '_pred'] = pmean
                df[imtstr + '_pred_sigma'] = pstddev[0]
                df[imtstr + '_pred_rock'] = pmean_rock
                df[imtstr + '_pred_sigma_rock'] = pstddev_rock[0]
                df[imtstr + '_pred_soil'] = pmean_soil
                df[imtstr + '_pred_sigma_soil'] = pstddev_soil[0]

                if imtstr != 'MMI':
                    total_only = self.gmpe_total_sd_only
                    tau_guess = SM_CONSTS['default_stddev_inter']
                else:
                    total_only = self.ipe_total_sd_only
                    tau_guess = SM_CONSTS['default_stddev_inter_mmi']
                if total_only:
                    df[imtstr + '_pred_tau'] = tau_guess * pstddev[0]
                    df[imtstr + '_pred_phi'] = np.sqrt(
                        pstddev[0]**2 - df[imtstr + '_pred_tau']**2)
                else:
                    df[imtstr + '_pred_tau'] = pstddev[1]
                    df[imtstr + '_pred_phi'] = pstddev[2]
                #
                # Compute the total residual
                #
                df[imtstr + '_residual'] = \
                    df[imtstr] - df[imtstr + '_pred']
                # -------------------------------------------------------------
                # Do the outlier flagging if we have a fault, or we don't
                # have a fault but the event magnitude is under the limit
                # -------------------------------------------------------------
                if self.do_outliers and \
                        (not isinstance(self.rupture_obj, PointRupture) or
                         self.rx.mag <= self.outlier_max_mag):
                    #
                    # Make a boolean array of stations that have been
                    # manually rehabilitated by the operator
                    #
                    is_valid = np.full(np.shape(df['id']), False, dtype=bool)
                    for valid in self.outlier_valid_stations:
                        is_valid |= valid == df['id']
                    #
                    # turn off nan warnings for this statement
                    #
                    np.seterr(invalid='ignore')
                    flagged = (np.abs(df[imtstr + '_residual']) >
                               self.outlier_deviation_level *
                               df[imtstr + '_pred_sigma']) & (~is_valid)
                    np.seterr(invalid='warn')
                    #
                    # Add NaN values to the list of outliers
                    #
                    flagged |= np.isnan(df[imtstr + '_residual'])

                    self.logger.debug('IMT: %s, flagged: %d' %
                                      (imtstr, np.sum(flagged)))
                    df[imtstr + '_outliers'] = flagged
                else:
                    #
                    # Not doing outliers, but should still flag NaNs
                    #
                    flagged = np.isnan(df[imtstr + '_residual'])
                    df[imtstr + '_outliers'] = flagged
                #
                # If uncertainty hasn't been set for MMI, give it
                # the default value
                #
                if imtstr == 'MMI' and all(df['MMI_sd'] == 0):
                    df['MMI_sd'][:] = SM_CONSTS['default_mmi_stddev']
            #
            # Get the lons/lats in radians while we're at it
            #
            df['lon_rad'] = np.radians(df['lon'])
            df['lat_rad'] = np.radians(df['lat'])
            if self.flip_lons:
                df['lon_rad'][df['lon_rad'] < 0] += 2 * np.pi
            #
            # It will be handy later on to have the rupture distance
            # in the dataframes
            #
            dd = get_distance(['rrup'], df['lat'], df['lon'],
                              df['depth'], self.rupture_obj)
            df['rrup'] = dd['rrup']
            if dd['rrup_var'] is not None:
                df['rrup_var'] = dd['rrup_var']
            else:
                df['rrup_var'] = np.zeros_like(dd['rrup'])

    def _deriveIMTsFromMMI(self):
        """
        Compute all the IMTs possible from MMI
        TODO: This logic needs to be revisited. We should probably make what
        we have to to do the CMS to make the needed output IMTs, but
        for now, we're just going to use what we have and the ccf.
        """
        if 'df2' not in self.dataframes:
            return

        df2 = self.df2.df
        for gmice_imt in self.gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
            if imt.SA == gmice_imt:
                iterlist = self.gmice.DEFINED_FOR_SA_PERIODS
            else:
                iterlist = [None]
            for period in iterlist:
                if period:
                    oqimt = gmice_imt(period)
                else:
                    oqimt = gmice_imt()
                imtstr = str(oqimt)

                np.seterr(invalid='ignore')
                df2[imtstr], _ = self.gmice.getGMfromMI(df2['MMI'], oqimt,
                                                        dists=df2['rrup'],
                                                        mag=self.rx.mag)
                df2[imtstr][df2['MMI'] < SM_CONSTS['min_mmi_convert']] = \
                    np.nan
                np.seterr(invalid='warn')
                df2[imtstr + '_sd'] = \
                    np.full_like(df2['MMI'], self.gmice.getMI2GMsd()[oqimt])
                self.df2.imts.add(imtstr)
                #
                # Get the predictions and stddevs
                #
                not_supported = False
                try:
                    gmpe = MultiGMPE.from_config(self.config, filter_imt=oqimt)
                except KeyError:
                    self.logger.warn(
                        "Input IMT %s not supported by GMPE: ignoring" %
                        imtstr)
                    not_supported = True
                if not_supported:
                    pmean = np.full_like(df2['MMI'], np.nan)
                    pstddev = [None] * 3
                    pstddev[0] = np.full_like(df2['MMI'], np.nan)
                    pstddev[1] = np.full_like(df2['MMI'], np.nan)
                    pstddev[2] = np.full_like(df2['MMI'], np.nan)
                else:
                    pmean, pstddev = self._gmas(gmpe, self.df2.sx,
                                                self.df2.dx, oqimt,
                                                self.apply_gafs)
                    pmean_rock, pstddev_rock = self._gmas(gmpe,
                                                          self.df2.sx_rock,
                                                          self.df2.dx, oqimt,
                                                          self.apply_gafs)
                    pmean_soil, pstddev_soil = self._gmas(gmpe,
                                                          self.df2.sx_soil,
                                                          self.df2.dx, oqimt,
                                                          self.apply_gafs)
                df2[imtstr + '_pred'] = pmean
                df2[imtstr + '_pred_sigma'] = pstddev[0]
                df2[imtstr + '_pred_rock'] = pmean_rock
                df2[imtstr + '_pred_sigma_rock'] = pstddev_rock[0]
                df2[imtstr + '_pred_soil'] = pmean_soil
                df2[imtstr + '_pred_sigma_soil'] = pstddev_soil[0]
                if imtstr != 'MMI':
                    total_only = self.gmpe_total_sd_only
                    tau_guess = SM_CONSTS['default_stddev_inter']
                else:
                    total_only = self.ipe_total_sd_only
                    tau_guess = SM_CONSTS['default_stddev_inter_mmi']
                if total_only:
                    df2[imtstr + '_pred_tau'] = tau_guess * pstddev[0]
                    df2[imtstr + '_pred_phi'] = np.sqrt(
                        pstddev[0]**2 - df2[imtstr + '_pred_tau']**2)
                else:
                    df2[imtstr + '_pred_tau'] = pstddev[1]
                    df2[imtstr + '_pred_phi'] = pstddev[2]
                df2[imtstr + '_residual'] = df2[imtstr] - pmean
                df2[imtstr + '_outliers'] = np.isnan(df2[imtstr + '_residual'])

    def _deriveMMIFromIMTs(self):
        """
        Make derived MMI from each of the IMTs in the input (for
        which the GMICE is defined; then select the best MMI for
        each station based on a list of "preferred" IMTs; also
        calculate the predicted MMI and the residual.
        """
        if 'df1' not in self.dataframes:
            return
        df1 = self.df1.df
        gmice_imts = self.gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES
        gmice_pers = self.gmice.DEFINED_FOR_SA_PERIODS
        np.seterr(invalid='ignore')
        df1['MMI'] = self.gmice.getPreferredMI(df1, dists=df1['rrup'],
                                               mag=self.rx.mag)
        np.seterr(invalid='warn')
        df1['MMI_sd'] = self.gmice.getPreferredSD()
        if df1['MMI_sd'] is not None:
            df1['MMI_sd'] = np.full_like(df1['lon'], df1['MMI_sd'])
        for imtstr in self.df1.imts:
            oqimt = imt.from_string(imtstr)
            if not isinstance(oqimt, tuple(gmice_imts)):
                continue
            if isinstance(oqimt, imt.SA) and \
               oqimt.period not in gmice_pers:
                continue

            np.seterr(invalid='ignore')
            df1['derived_MMI_from_' + imtstr], _ = \
                self.gmice.getMIfromGM(df1[imtstr], oqimt,
                                       dists=df1['rrup'],
                                       mag=self.rx.mag)
            np.seterr(invalid='warn')
            df1['derived_MMI_from_' + imtstr + '_sd'] = \
                np.full_like(df1[imtstr], self.gmice.getGM2MIsd()[oqimt])

        preferred_imts = ['PGV', 'PGA', 'SA(1.0)', 'SA(0.3)', 'SA(3.0']
        if(df1['MMI'] is None):
            df1['MMI'] = np.full_like(df1['lon'], np.nan)
            df1['MMI_sd'] = np.full_like(df1['lon'], np.nan)
        for imtstr in preferred_imts:
            if 'derived_MMI_from_' + imtstr in df1:
                ixx = np.isnan(df1['MMI'])
                df1['MMI'][ixx] = df1['derived_MMI_from_' + imtstr][ixx]
                df1['MMI_sd'][ixx] = \
                    df1['derived_MMI_from_' + imtstr + '_sd'][ixx]
        self.df1.imts.add('MMI')
        #
        # Get the prediction and stddevs
        #
        gmpe = None
        pmean, pstddev = self._gmas(
            gmpe, self.df1.sx, self.df1.dx, imt.from_string('MMI'),
            self.apply_gafs)
        pmean_rock, pstddev_rock = self._gmas(
            gmpe, self.df1.sx_rock, self.df1.dx, imt.from_string('MMI'),
            self.apply_gafs)
        pmean_soil, pstddev_soil = self._gmas(
            gmpe, self.df1.sx_soil, self.df1.dx, imt.from_string('MMI'),
            self.apply_gafs)
        df1['MMI' + '_pred'] = pmean
        df1['MMI' + '_pred_sigma'] = pstddev[0]
        df1['MMI' + '_pred_rock'] = pmean_rock
        df1['MMI' + '_pred_sigma_rock'] = pstddev_rock[0]
        df1['MMI' + '_pred_soil'] = pmean_soil
        df1['MMI' + '_pred_sigma_soil'] = pstddev_soil[0]
        if self.ipe_total_sd_only:
            tau_guess = SM_CONSTS['default_stddev_inter_mmi']
            df1['MMI' + '_pred_tau'] = tau_guess * pstddev[0]
            df1['MMI' + '_pred_phi'] = np.sqrt(
                pstddev[0]**2 - df1['MMI' + '_pred_tau']**2)
        else:
            df1['MMI' + '_pred_tau'] = pstddev[1]
            df1['MMI' + '_pred_phi'] = pstddev[2]
        df1['MMI' + '_residual'] = df1['MMI'] - pmean
        df1['MMI' + '_outliers'] = np.isnan(df1['MMI' + '_residual'])

    def _fillDataArrays(self):
        """
        For each IMT get lists of the amplitudes that can contribute
        to the bias and the interpolation. Keep lists of, IMT, period
        index, lons, lats, residuals, tau, phi, additional uncertainty,
        and rupture distance.
        """
        #
        # Get the valid imts for each station
        #
        imtsets = {}
        sasets = {}
        for ndf in self.dataframes:
            imtsets[ndf], sasets[ndf] = _get_imt_lists(getattr(self, ndf))

        for imtstr in self.combined_imt_set:
            #
            # Fill the station arrays; here we use lists and append to
            # them because it is much faster than appending to a numpy
            # array; after the loop, the lists are converted to numpy
            # arrays:
            #
            imtlist = []    # The input IMT
            period_ix = []  # The index of the (pseudo-)period of the input IMT
            lons_rad = []   # longitude (in radians) of the input station
            lats_rad = []   # latitude (in radians) of the input station
            resids = []     # The residual of the input IMT
            tau = []        # The between-event stddev of the input IMT
            phi = []        # The within-event stddev of the input IMT
            sig_extra = []  # Additional stddev of the input IMT
            rrups = []      # The rupture distance of the input station
            for ndf in self.dataframes:
                sdf = getattr(self, ndf).df
                for i in range(np.size(sdf['lon'])):
                    #
                    # Each station can provide 0, 1, or 2 IMTs:
                    #
                    for imtin in _get_nearest_imts(imtstr, imtsets[ndf][i],
                                                   sasets[ndf][i]):
                        imtlist.append(imtin)
                        period_ix.append(self.imt_per_ix[imtin])
                        lons_rad.append(sdf['lon_rad'][i])
                        lats_rad.append(sdf['lat_rad'][i])
                        resids.append(sdf[imtin + '_residual'][i])
                        tau.append(sdf[imtin + '_pred_tau'][i])
                        phi.append(sdf[imtin + '_pred_phi'][i])
                        sig_extra.append(sdf[imtin + '_sd'][i])
                        rrups.append(sdf['rrup'][i])
            self.sta_imtstr[imtstr] = imtlist.copy()
            self.sta_period_ix[imtstr] = np.array(period_ix).reshape((-1, 1))
            self.sta_lons_rad[imtstr] = np.array(lons_rad).reshape((-1, 1))
            self.sta_lats_rad[imtstr] = np.array(lats_rad).reshape((-1, 1))
            if self.flip_lons:
                self.sta_lons_rad[imtstr][self.sta_lons_rad[imtstr] < 0] += \
                    2 * np.pi
            self.sta_resids[imtstr] = np.array(resids).reshape((-1, 1))
            self.sta_tau[imtstr] = np.array(tau).reshape((-1, 1))
            self.sta_phi[imtstr] = np.array(phi).reshape((-1, 1))
            self.sta_sig_extra[imtstr] = np.array(sig_extra).reshape((-1, 1))
            self.sta_sig_total[imtstr] = np.sqrt(
                self.sta_phi[imtstr]**2 + self.sta_sig_extra[imtstr]**2)
            self.sta_rrups[imtstr] = np.array(rrups)

    def _computeBias(self):
        """
        Compute a bias for all of the IMTs in the inputs and outputs
        """
        for imtstr in self.combined_imt_set:
            time1 = time.time()
            if np.size(self.sta_lons_rad) == 0:
                self.bias_num[imtstr] = 0.0
                self.bias_den[imtstr] = 0.0
                continue
            #
            # Get the index of the (pseudo-) period of the output IMT
            #
            outperiod_ix = self.imt_per_ix[imtstr]
            #
            # Get the distance-limited set of data for use in computing
            # the bias
            #
            dix = self.sta_rrups[imtstr] <= self.bias_max_range
            sta_phi_dl = self.sta_phi[imtstr][dix].reshape((-1, 1))
            sta_tau_dl = self.sta_tau[imtstr][dix].reshape((-1, 1))
            sta_sig_total_dl = self.sta_sig_total[imtstr][dix].reshape((-1, 1))
            sta_lons_rad_dl = self.sta_lons_rad[imtstr][dix].ravel()
            sta_lats_rad_dl = self.sta_lats_rad[imtstr][dix].ravel()
            sta_period_ix_dl = self.sta_period_ix[imtstr][dix].reshape((-1, 1))
            sta_resids_dl = self.sta_resids[imtstr][dix].reshape((-1, 1))
            if np.size(sta_lons_rad_dl) == 0:
                self.bias_num[imtstr] = 0.0
                self.bias_den[imtstr] = 0.0
                continue
            #
            # This builds the omega factors to apply to the covariance
            #
            corr_adj = sta_phi_dl / sta_sig_total_dl
            corr_adj22 = corr_adj * corr_adj.T
            np.fill_diagonal(corr_adj22, 1.0)
            #
            # Build the covariance matrix of the residuals
            #
            nsta = np.size(sta_lons_rad_dl)
            matrix22 = np.empty((nsta, nsta), dtype=np.double)
            geodetic_distance_fast(sta_lons_rad_dl,
                                   sta_lats_rad_dl,
                                   sta_lons_rad_dl,
                                   sta_lats_rad_dl, matrix22)
            ones = np.ones((1, nsta), dtype=np.long)
            t1_22 = sta_period_ix_dl * ones
            t2_22 = sta_period_ix_dl.T * ones.T
            self.ccf.getCorrelation(t1_22, t2_22, matrix22)
            sta_phi_dl_flat = sta_phi_dl.ravel()
            make_sigma_matrix(matrix22, corr_adj22,
                              sta_phi_dl_flat,
                              sta_phi_dl_flat)
            sigma22inv = np.linalg.pinv(matrix22)
            #
            # Compute the bias numerator and denominator pieces
            #
            if self.do_bias and (not isinstance(self.rupture_obj, PointRupture)
                                 or self.rx.mag <= self.bias_max_mag):
                #
                # Get the correlation between the inputs and outputs
                #
                out_ix_arr = np.full_like(sta_period_ix_dl, outperiod_ix)
                Z = np.zeros_like(out_ix_arr, dtype=float)
                self.ccf.getCorrelation(sta_period_ix_dl, out_ix_arr, Z)
                #
                # Scale the correlation factor (Z) by the correlation
                # adjustment due to observational uncertainty
                #
                Z *= corr_adj
                #
                # Compute the bias denominator and numerator terms
                #
                ztmp = Z.T.dot(sigma22inv)
                self.bias_num[imtstr] = ztmp.dot(Z * sta_resids_dl)[0][0]
                self.bias_den[imtstr] = ztmp.dot(Z)[0][0]
            else:
                self.bias_num[imtstr] = 0.0
                self.bias_den[imtstr] = 0.0

            nom_tau = np.mean(sta_tau_dl.ravel())
            nom_variance = 1.0 / ((1.0 / nom_tau**2) + self.bias_den[imtstr])
            self.nominal_bias[imtstr] = self.bias_num[imtstr] * nom_variance
            bias_time = time.time() - time1
            #
            # Write the nominal values of the bias and its stddev to log
            #
            self.logger.debug(
                '%s: nom bias %f nom stddev %f; %d stations (time=%f sec)'
                % (imtstr, self.nominal_bias[imtstr], np.sqrt(nom_variance),
                   np.size(sta_lons_rad_dl), bias_time))

    def _computeDirectivityPredictionLocations(self):
        """
        Figure out if we need the directivity factors, and if so, pre-calculate
        them. These will be used later in _computeMVN.
        """
        if self.do_directivity is True:
            self.logger.info('Directivity for prediction locations...')
            time1 = time.time()

            # Precompute directivity at all periods
            dir_out = Rowshandel2013(
                self.rupture_obj._origin, self.rupture_obj,
                self.lats.reshape((1, -1)),
                self.lons.reshape((1, -1)),
                np.zeros_like((len(self.lats), 1)),
                dx=1.0,
                T=Rowshandel2013.getPeriods(),
                a_weight=0.5,
                mtype=1
            )
            self.dir_results.append((dir_out, self.dx_out))
            # Precompute directivity for the attenuation curves
            dir_out = Rowshandel2013(
                self.rupture_obj._origin, self.rupture_obj,
                self.atten_coords['lats'].reshape((1, -1)),
                self.atten_coords['lons'].reshape((1, -1)),
                np.zeros_like((len(self.atten_coords['lats']), 1)),
                dx=1.0,
                T=Rowshandel2013.getPeriods(),
                a_weight=0.5,
                mtype=1
            )
            self.dir_results.append((dir_out, self.atten_dx))
            directivity_time = time.time() - time1
            self.logger.debug(
                'Directivity prediction evaluation time: %f sec'
                % directivity_time)
        else:
            self.directivity = None

    def _computeMVN(self, imtstr):
        """
        Do the MVN computations
        """
        self.logger.debug('computeMVN: doing IMT %s' % imtstr)
        time1 = time.time()
        #
        # Get the index of the (pesudo-) period of the output IMT
        #
        outperiod_ix = self.imt_per_ix[imtstr]

        #
        # Get the predictions at the output points
        #
        oqimt = imt.from_string(imtstr)
        if imtstr != 'MMI':
            gmpe = MultiGMPE.from_config(self.config, filter_imt=oqimt)
        else:
            gmpe = self.ipe

        pout_mean, pout_sd = self._gmas(
            gmpe, self.sx_out, self.dx_out, oqimt, self.apply_gafs)

        if self.use_simulations:
            if imtstr == 'MMI':
                pout_mean = self.gmice.getPreferredMI(self.sim_df,
                                                      dists=self.dx_out.rrup,
                                                      mag=self.rx.mag)
            else:
                pout_mean = self.sim_df[imtstr]

        #
        # While we have the gmpe for this IMT, we should make
        # the attenuation curves
        #
        x_mean, x_sd = self._gmas(gmpe,
                                  self.atten_sx_rock,
                                  self.atten_dx,
                                  oqimt, self.apply_gafs)
        self.atten_rock_mean[imtstr] = x_mean
        self.atten_rock_sd[imtstr] = x_sd[0]
        x_mean, x_sd = self._gmas(gmpe,
                                  self.atten_sx_soil,
                                  self.atten_dx,
                                  oqimt, self.apply_gafs)
        self.atten_soil_mean[imtstr] = x_mean
        self.atten_soil_sd[imtstr] = x_sd[0]

        #
        # If there are no data, just use the unbiased prediction
        # and the total stddev
        #
        if np.size(self.sta_lons_rad[imtstr]) == 0:
            self.outgrid[imtstr] = pout_mean
            self.outsd[imtstr] = pout_sd[0]
            return

        #
        # Get an array of the within-event standard deviations for the
        # output IMT at the output points
        #
        if imtstr != 'MMI':
            total_only = self.gmpe_total_sd_only
            tau_guess = SM_CONSTS['default_stddev_inter']
        else:
            total_only = self.ipe_total_sd_only
            tau_guess = SM_CONSTS['default_stddev_inter_mmi']
        if total_only:
            self.tsd[imtstr] = tau_guess * pout_sd[0]
            self.psd[imtstr] = np.sqrt(pout_sd[0]**2 - self.tsd[imtstr]**2)
            self.psd_raw[imtstr] = np.sqrt(pout_sd[1]**2 - self.tsd[imtstr]**2)
        else:
            self.psd[imtstr] = pout_sd[2]
            self.psd_raw[imtstr] = pout_sd[5]
            self.tsd[imtstr] = pout_sd[1]
        pout_sd2 = np.power(self.psd[imtstr], 2.0)

        #
        # Bias the predictions, and add the residual variance to
        # phi
        #
        out_bias_var = 1.0 / ((1.0 / self.tsd[imtstr]**2) +
                              self.bias_den[imtstr])
        out_bias = self.bias_num[imtstr] * out_bias_var
        pout_mean += out_bias.reshape(pout_mean.shape)
        self.psd[imtstr] = np.sqrt(self.psd[imtstr]**2 + out_bias_var)
        pout_sd2 += out_bias_var.reshape(pout_sd2.shape)
        self.psd_raw[imtstr] = np.sqrt(
            self.psd_raw[imtstr]**2 + out_bias_var.reshape(pout_sd2.shape))

        #
        # Unbias the station residuals and compute the
        # new phi that includes the variance of the bias
        #
        for i in range(np.size(self.sta_lons_rad[imtstr])):
            imtin = self.sta_imtstr[imtstr][i]
            in_bias_var = 1.0 / ((1.0 / self.sta_tau[imtstr][i, 0]**2) +
                                 self.bias_den[imtin])
            in_bias = self.bias_num[imtin] * in_bias_var
            self.sta_resids[imtstr][i, 0] -= in_bias
            self.sta_phi[imtstr][i, 0] = np.sqrt(
                self.sta_phi[imtstr][i, 0]**2 + in_bias_var)

        #
        # Update the omega factors to account for the bias and the
        # new value of phi
        #
        corr_adj = self.sta_phi[imtstr] / np.sqrt(
            self.sta_phi[imtstr]**2 + self.sta_sig_extra[imtstr]**2)
        corr_adj22 = corr_adj * corr_adj.T
        np.fill_diagonal(corr_adj22, 1.0)

        sta_lons_rad_flat = self.sta_lons_rad[imtstr].ravel()
        sta_lats_rad_flat = self.sta_lats_rad[imtstr].ravel()
        #
        # Re-build the covariance matrix of the residuals with
        # the full set of data
        #
        nsta = np.size(sta_lons_rad_flat)
        matrix22 = np.empty((nsta, nsta), dtype=np.double)
        geodetic_distance_fast(sta_lons_rad_flat,
                               sta_lats_rad_flat,
                               sta_lons_rad_flat,
                               sta_lats_rad_flat,
                               matrix22)
        ones = np.ones((1, nsta), dtype=np.long)
        t1_22 = self.sta_period_ix[imtstr] * ones
        t2_22 = self.sta_period_ix[imtstr].T * ones.T
        self.ccf.getCorrelation(t1_22, t2_22, matrix22)

        #
        # Rebuild sigma22_inv now that we have updated phi and
        # the correlation adjustment factors
        #
        sta_phi_flat = self.sta_phi[imtstr].ravel()
        make_sigma_matrix(matrix22, corr_adj22, sta_phi_flat, sta_phi_flat)
        sigma22inv = np.linalg.pinv(matrix22)

        #
        # Now do the MVN itself...
        #
        dtime = mtime = ddtime = ctime = stime = atime = 0

        ampgrid = np.zeros_like(pout_mean)
        sdgrid = np.zeros_like(pout_mean)
        corr_adj12 = corr_adj.T * np.ones((self.smnx, 1))
        # Stuff that doesn't change within the loop:
        lons_out_rad_flat = self.lons_out_rad.ravel()
        lats_out_rad_flat = self.lats_out_rad.ravel()
        d12_cols = self.smnx
        t2_12 = np.full((d12_cols, nsta), outperiod_ix, dtype=np.long)
        t1_12 = self.sta_period_ix[imtstr].T * np.ones((d12_cols, 1),
                                                       dtype=np.long)
        # sdsta is the standard deviation of the stations
        sdsta = self.sta_phi[imtstr].ravel()
        matrix12 = np.empty(t2_12.shape, dtype=np.double)
        rcmatrix = np.empty(t2_12.shape, dtype=np.double)
        for iy in range(self.smny):
            ss = iy * self.smnx
            se = (iy + 1) * self.smnx
            time4 = time.time()
            geodetic_distance_fast(
                sta_lons_rad_flat,
                sta_lats_rad_flat,
                lons_out_rad_flat[ss:se],
                lats_out_rad_flat[ss:se],
                matrix12)
            ddtime += time.time() - time4
            time4 = time.time()
            self.ccf.getCorrelation(t1_12, t2_12, matrix12)
            ctime += time.time() - time4
            time4 = time.time()
            # sdarr is the standard deviation of the output sites
            sdarr = self.psd[imtstr][iy, :].ravel()
            make_sigma_matrix(matrix12, corr_adj12, sdsta, sdarr)
            stime += time.time() - time4
            time4 = time.time()
            #
            # Sigma12 * Sigma22^-1 is known as the 'regression
            # coefficient' matrix (rcmatrix)
            #
            np.dot(matrix12, sigma22inv, out=rcmatrix)
            dtime += time.time() - time4
            time4 = time.time()

            #
            # This is the MVN solution for the conditional mean
            #
            adj_resid = corr_adj * self.sta_resids[imtstr]
            ampgrid[iy, :] = \
                pout_mean[iy, :] + rcmatrix.dot(adj_resid).reshape((-1,))
            atime += time.time() - time4
            time4 = time.time()
            #
            # We only want the diagonal elements of the conditional
            # covariance matrix, so there is no point in doing the
            # full solution with the dot product, e.g.:
            # sdgrid[ss:se] = pout_sd2[ss:se] -
            #       np.diag(rcmatrix.dot(sigma21))
            #
            make_sd_array(sdgrid, pout_sd2, iy, rcmatrix, matrix12)
            mtime += time.time() - time4

        self.outgrid[imtstr] = ampgrid
        self.outsd[imtstr] = sdgrid

        self.logger.debug('\ttime for %s distance=%f' % (imtstr, ddtime))
        self.logger.debug('\ttime for %s correlation=%f' % (imtstr, ctime))
        self.logger.debug('\ttime for %s sigma=%f' % (imtstr, stime))
        self.logger.debug('\ttime for %s rcmatrix=%f' % (imtstr, dtime))
        self.logger.debug('\ttime for %s amp calc=%f' % (imtstr, atime))
        self.logger.debug('\ttime for %s sd calc=%f' % (imtstr, mtime))
        self.logger.debug('total time for %s=%f' %
                          (imtstr, time.time() - time1))

    def _getLandMask(self):
        """
        Get the landmask for this map. Land will be False, water will
        be True (because of the way masked arrays work).
        """
        if not self.do_grid:
            return None
        gd = GeoDict.createDictFromBox(self.W, self.E, self.S, self.N,
                                       self.smdx, self.smdy)
        bbox = (gd.xmin, gd.ymin, gd.xmax, gd.ymax)
        if 'CALLED_FROM_PYTEST' in os.environ:
            return np.zeros((gd.ny, gd.nx), dtype=np.bool)

        oceans = shpreader.natural_earth(category='physical',
                                         name='ocean',
                                         resolution='10m')
        with fiona.open(oceans) as c:
            tshapes = list(c.items(bbox=bbox))
            shapes = []
            for tshp in tshapes:
                shapes.append(shape(tshp[1]['geometry']))
            if len(shapes):
                oceangrid = Grid2D.rasterizeFromGeometry(shapes, gd,
                                                         fillValue=0.0)
                return oceangrid.getData().astype(np.bool)
            else:
                return np.zeros((gd.ny, gd.nx), dtype=np.bool)

    def _getMaskedGrids(self, bmask):
        """
        For each grid in the output, generate a grid with the water areas
        masked out.
        """
        moutgrid = {}
        if not self.do_grid:
            for imtout in self.imt_out_set:
                moutgrid[imtout] = self.outgrid[imtout]
            return moutgrid
        for imtout in self.imt_out_set:
            moutgrid[imtout] = \
                ma.masked_array(self.outgrid[imtout], mask=bmask)
        return moutgrid

    def _getInfo(self, moutgrid):
        """
        Create an info dictionary that can be made into the info.json file.
        """
        #
        # Get the map grade
        #
        mean_rat, mygrade = _get_map_grade(
            self.do_grid, self.outsd, self.psd_raw, moutgrid)
        # ---------------------------------------------------------------------
        # This is the metadata for creating info.json
        # ---------------------------------------------------------------------
        st = 'strec'
        ip = 'input'
        ei = 'event_information'
        op = 'output'
        gm = 'ground_motions'
        mi = 'map_information'
        un = 'uncertainty'
        pp = 'processing'
        gmm = 'ground_motion_modules'
        ms = 'miscellaneous'
        sv = 'shakemap_versions'
        sr = 'site_response'
        info = self._info
        info[ip] = {}
        info[ip][ei] = {}
        info[ip][ei]['depth'] = str(self.rx.hypo_depth)
        info[ip][ei]['event_id'] = self._eventid

        # look for the presence of a strec_results file and read it in
        _, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        strecfile = os.path.join(datadir, 'strec_results.json')
        if os.path.isfile(strecfile):
            strec_results = json.load(open(strecfile, 'rt'))
            info[st] = strec_results

        # the following items are primarily useful for PDL
        origin = self.rupture_obj._origin
        info[ip][ei]['eventsource'] = origin.netid
        info[ip][ei]['netid'] = origin.netid
        # The netid could be a valid part of the eventsourcecode, so we have
        # to check here if it ***starts with*** the netid
        if origin.id.startswith(origin.netid):
            eventsourcecode = origin.id.replace(origin.netid, '', 1)
        else:
            eventsourcecode = origin.id
        info[ip][ei]['eventsourcecode'] = eventsourcecode
        info[ip][ei]['id'] = origin.id
        info[ip][ei]['productcode'] = origin.productcode
        info[ip][ei]['productsource'] = self.config['system']['source_network']
        info[ip][ei]['producttype'] = self.config['system']['product_type']

        info[ip][ei]['event_ref'] = getattr(origin, 'reference', None)
        info[ip][ei]['fault_ref'] = self.rupture_obj.getReference()
        if 'df2' in self.dataframes:
            info[ip][ei]['intensity_observations'] = \
                str(np.size(self.df2.df['lon']))
        else:
            info[ip][ei]['intensity_observations'] = '0'
        info[ip][ei]['latitude'] = str(self.rx.hypo_lat)
        info[ip][ei]['longitude'] = str(self.rx.hypo_lon)
        info[ip][ei]['location'] = origin.locstring
        info[ip][ei]['magnitude'] = str(self.rx.mag)
        info[ip][ei]['origin_time'] = origin.time.strftime(constants.TIMEFMT)
        if 'df1' in self.dataframes:
            info[ip][ei]['seismic_stations'] = \
                str(np.size(self.df1.df['lon']))
        else:
            info[ip][ei]['seismic_stations'] = '0'
        info[ip][ei]['src_mech'] = origin.mech
        # This AND locaction?
        info[ip][ei]['event_description'] = origin.locstring
        # This AND src_mech?
        # look at the origin information for indications that this
        # event is a scenario
        condition1 = hasattr(
            origin, 'event_type') and origin.event_type.lower() == 'scenario'
        condition2 = origin.id.endswith('_se')
        if condition1 or condition2:
            info[ip][ei]['event_type'] = 'SCENARIO'
        else:
            info[ip][ei]['event_type'] = 'ACTUAL'

        info[op] = {}
        info[op][gm] = {}
        for myimt in self.imt_out_set:
            info[op][gm][myimt] = {}
            if myimt == 'MMI':
                units = 'intensity'
            elif myimt == 'PGV':
                units = 'cms'
            else:
                units = 'g'
            info[op][gm][myimt]['units'] = units
            if myimt in self.nominal_bias:
                info[op][gm][myimt]['bias'] = \
                    _string_round(self.nominal_bias[myimt], 3)
            else:
                info[op][gm][myimt]['bias'] = None
            if myimt == 'MMI' or myimt == 'PGV':
                info[op][gm][myimt]['max_grid'] = \
                    _string_round(np.max(self.outgrid[myimt]), 3)
                info[op][gm][myimt]['max'] = \
                    _string_round(np.max(moutgrid[myimt]), 3)
            else:
                info[op][gm][myimt]['max_grid'] = \
                    _string_round(np.exp(np.max(self.outgrid[myimt])), 3)
                info[op][gm][myimt]['max'] = \
                    _string_round(np.exp(np.max(moutgrid[myimt])), 3)

        info[op][mi] = {}
        info[op][mi]['grid_points'] = {}
        info[op][mi]['grid_points']['longitude'] = str(self.smnx)
        info[op][mi]['grid_points']['latitude'] = str(self.smny)
        info[op][mi]['grid_points']['units'] = ''
        info[op][mi]['grid_spacing'] = {}
        info[op][mi]['grid_spacing']['longitude'] = _string_round(self.smdx, 7)
        info[op][mi]['grid_spacing']['latitude'] = _string_round(self.smdy, 7)
        info[op][mi]['grid_spacing']['units'] = 'degrees'
        info[op][mi]['grid_span'] = {}
        if self.E <= 0 and self.W >= 0:
            info[op][mi]['grid_span']['longitude'] = \
                _string_round(self.E + 360.0 - self.W, 3)
        else:
            info[op][mi]['grid_span']['longitude'] = \
                _string_round(self.E - self.W, 3)
        info[op][mi]['grid_span']['latitude'] = \
            _string_round(self.N - self.S, 3)
        info[op][mi]['grid_span']['units'] = 'degrees'
        info[op][mi]['min'] = {}
        info[op][mi]['min']['longitude'] = _string_round(self.W, 3)
        info[op][mi]['min']['latitude'] = _string_round(self.S, 3)
        info[op][mi]['min']['units'] = 'degrees'
        info[op][mi]['max'] = {}
        info[op][mi]['max']['longitude'] = _string_round(self.E, 3)
        info[op][mi]['max']['latitude'] = _string_round(self.N, 3)
        info[op][mi]['max']['units'] = 'degrees'
        info[op][un] = {}
        info[op][un]['grade'] = mygrade
        info[op][un]['mean_uncertainty_ratio'] = _string_round(mean_rat, 3)
        if 'df2' in self.dataframes:
            info[op][un]['total_flagged_mi'] = \
                str(np.sum(self.df2.df['MMI_outliers'] |
                           np.isnan(self.df2.df['MMI'])))
        else:
            info[op][un]['total_flagged_mi'] = '0'
        if 'df1' in self.dataframes:
            all_flagged = np.full(self.df1.df['lon'].shape, False,
                                  dtype=np.bool)
            for imtstr in self.df1.imts:
                if 'MMI' in imtstr:
                    continue
                all_flagged |= \
                    self.df1.df[imtstr + '_outliers'] | \
                    np.isnan(self.df1.df[imtstr])
            all_flagged |= self.df1.df['flagged']
            info[op][un]['total_flagged_pgm'] = str(np.sum(all_flagged))
        else:
            info[op][un]['total_flagged_pgm'] = '0'
        info[pp] = {}
        info[pp][gmm] = {}
        info[pp][gmm]['gmpe'] = {}
        info[pp][gmm]['gmpe']['module'] = str(self.config['modeling']['gmpe'])
        info[pp][gmm]['gmpe']['reference'] = ''
        info[pp][gmm]['ipe'] = {}
        info[pp][gmm]['ipe']['module'] = \
            str(self.config['ipe_modules'][self.config['modeling']['ipe']][0])
        info[pp][gmm]['ipe']['reference'] = ''
        info[pp][gmm]['gmice'] = {}
        info[pp][gmm]['gmice']['module'] = str(
            self.config['gmice_modules'][self.config['modeling']['gmice']][0])
        info[pp][gmm]['gmice']['reference'] = ''
        info[pp][gmm]['ccf'] = {}
        info[pp][gmm]['ccf']['module'] = \
            str(self.config['ccf_modules'][self.config['modeling']['ccf']][0])
        info[pp][gmm]['ccf']['reference'] = ''
        info[pp][gmm]['basin_correction'] = {}
        info[pp][gmm]['basin_correction']['module'] = 'None'
        info[pp][gmm]['basin_correction']['reference'] = ''
        info[pp][gmm]['directivity'] = {}
        info[pp][gmm]['directivity']['module'] = 'None'
        info[pp][gmm]['directivity']['reference'] = ''
        info[pp][ms] = {}
        info[pp][ms]['bias_max_dsigma'] = str(self.bias_max_dsigma)
        info[pp][ms]['bias_max_mag'] = str(self.bias_max_mag)
        info[pp][ms]['bias_max_range'] = str(self.bias_max_range)
        info[pp][ms]['median_dist'] = 'True'
        info[pp][ms]['do_outliers'] = self.do_outliers
        info[pp][ms]['outlier_deviation_level'] = str(
            self.outlier_deviation_level)
        info[pp][ms]['outlier_max_mag'] = str(self.outlier_max_mag)
        info[pp][sv] = {}
        info[pp][sv]['shakemap_revision'] = get_versions()['version']
        info[pp][sv]['shakemap_revision_id'] = \
            get_versions()['full-revisionid']
        info[pp][sv]['process_time'] = \
            strftime(constants.ALT_TIMEFMT, gmtime())
        info[pp][sv]['map_version'] = \
            self.ic.getVersionHistory()['history'][-1][2]
        info[pp][sv]['map_comment'] = \
            self.ic.getVersionHistory()['history'][-1][3]
        info[pp][sv]['map_data_history'] = \
            self.ic.getVersionHistory()['history']
        info[pp][sv]['map_status'] = self.config['system']['map_status']
        info[pp][sr] = {}
        info[pp][sr]['vs30default'] = str(self.vs30default)
        info[pp][sr]['site_correction'] = 'GMPE native'
        return info

    def _fillStationJSON(self):
        """
        Get the station JSON dictionary and then add a bunch of stuff to it.
        """
        if not hasattr(self, 'stations') or self.stations is None:
            return {'eventid': self._eventid,
                    'features': []}
        sjdict = {}
        # ---------------------------------------------------------------------
        # Compute a bias for all the IMTs in the data frames
        # ---------------------------------------------------------------------
        for ndf in self.dataframes:
            sdf = getattr(self, ndf).df
            imts = getattr(self, ndf).imts
            for myimt in imts:
                mybias = self.bias_num[myimt] / \
                    ((1.0 / sdf[myimt + '_pred_tau']**2) +
                     self.bias_den[myimt])
                sdf[myimt + '_bias'] = mybias.ravel()

        # ---------------------------------------------------------------------
        # Add the station data. The stationlist object has the original
        # data and produces a GeoJSON object (a dictionary, really), but
        # we need to add peak values and flagging that has been done here.
        # ---------------------------------------------------------------------
        #
        # First make a dictionary of distances
        #
        dist_dict = {
            'df1': {},
            'df2': {}
        }
        for ndf in self.dataframes:
            dx = getattr(self, ndf).dx
            for dm in get_distance_measures():
                dm_arr = getattr(dx, dm, None)
                if dm_arr is not None:
                    dist_dict[ndf][dm] = dm_arr
                else:
                    continue
                if dm in ('rrup', 'rjb'):
                    dm_var = getattr(dx, dm + '_var', None)
                    if dm_var is not None:
                        dist_dict[ndf][dm + '_var'] = dm_var
                    else:
                        dist_dict[ndf][dm + '_var'] = np.zeros_like(dm_arr)
        #
        # Get the index for each station ID
        #
        sjdict = self.stations.getGeoJson()
        sta_ix = {
            'df1': {},
            'df2': {}
        }
        for ndf in self.dataframes:
            sdf = getattr(self, ndf).df
            sta_ix[ndf] = dict(zip(sdf['id'], range(len(sdf['id']))))
        #
        # Now go through the GeoJSON and add various properties and
        # amps from the df_dict dictionaries
        #
        sjdict_copy = copy.copy(sjdict['features'])
        for station in sjdict_copy:
            if station['id'] in sta_ix['df1']:
                ndf = 'df1'
                station['properties']['station_type'] = 'seismic'
            elif station['id'] in sta_ix['df2']:
                ndf = 'df2'
                station['properties']['station_type'] = 'macroseismic'
            else:
                # We're probably using --no_seismic or --no_macroseismic
                if self.no_seismic or self.no_macroseismic:
                    sjdict['features'].remove(station)
                    continue
                else:
                    raise ValueError('Unknown station %s in stationlist' %
                                     (station['id']))
            dfx = getattr(self, ndf)
            sdf = dfx.df
            six = sta_ix[ndf][station['id']]
            #
            # Set the 'intensity', 'pga', and 'pga' peak properties
            #
            if 'MMI' in sdf and not sdf['MMI_outliers'][six] \
                    and not np.isnan(sdf['MMI'][six]):
                station['properties']['intensity'] = \
                    float("%.1f" % sdf['MMI'][six])
                station['properties']['intensity_stddev'] = \
                    sdf['MMI_sd'][six]
            else:
                station['properties']['intensity'] = 'null'
                station['properties']['intensity_stddev'] = 'null'

            if 'PGA' in sdf and not sdf['PGA_outliers'][six] \
                    and not np.isnan(sdf['PGA'][six]) \
                    and (ndf != 'df1' or not sdf['flagged'][six]):
                station['properties']['pga'] = \
                    _round_float(np.exp(sdf['PGA'][six]) * 100, 4)
            else:
                station['properties']['pga'] = 'null'

            if 'PGV' in sdf and not sdf['PGV_outliers'][six] \
                    and not np.isnan(sdf['PGV'][six]) \
                    and (ndf != 'df1' or not sdf['flagged'][six]):
                station['properties']['pgv'] = \
                    _round_float(np.exp(sdf['PGV'][six]), 4)
            else:
                station['properties']['pgv'] = 'null'
            #
            # Add vs30
            #
            station['properties']['vs30'] = _round_float(dfx.sx.vs30[six], 2)
            #
            # Add the predictions so we can plot residuals
            #
            station['properties']['predictions'] = []
            for key in sdf.keys():
                if not key.endswith('_pred'):
                    continue
                myamp = sdf[key][six]
                myamp_rock = sdf[key + '_rock'][six]
                myamp_soil = sdf[key + '_soil'][six]
                tau_str = 'ln_tau'
                phi_str = 'ln_phi'
                sigma_str = 'ln_sigma'
                sigma_str_rock = 'ln_sigma_rock'
                sigma_str_soil = 'ln_sigma_soil'
                bias_str = 'ln_bias'
                if key.startswith('PGV'):
                    value = np.exp(myamp)
                    value_rock = np.exp(myamp_rock)
                    value_soil = np.exp(myamp_soil)
                    units = 'cm/s'
                elif key.startswith('MMI'):
                    value = myamp
                    value_rock = myamp_rock
                    value_soil = myamp_soil
                    units = 'intensity'
                    tau_str = 'tau'
                    phi_str = 'phi'
                    sigma_str = 'sigma'
                    sigma_str_rock = 'sigma_rock'
                    sigma_str_soil = 'sigma_soil'
                    bias_str = 'bias'
                else:
                    value = np.exp(myamp) * 100
                    value_rock = np.exp(myamp_rock) * 100
                    value_soil = np.exp(myamp_soil) * 100
                    units = '%g'
                if self.gmpe_total_sd_only:
                    mytau = 0
                else:
                    mytau = sdf[key + '_tau'][six]
                myphi = sdf[key + '_phi'][six]
                mysigma = np.sqrt(mytau**2 + myphi**2)
                mysigma_rock = sdf[key + '_sigma_rock'][six]
                mysigma_soil = sdf[key + '_sigma_soil'][six]
                imt_name = key.lower().replace('_pred', '')
                mybias = sdf[imt_name.upper() + '_bias'][six]
                station['properties']['predictions'].append({
                    'name': imt_name,
                    'value': _round_float(value, 4),
                    'value_rock': _round_float(value_rock, 4),
                    'value_soil': _round_float(value_soil, 4),
                    'units': units,
                    tau_str: _round_float(mytau, 4),
                    phi_str: _round_float(myphi, 4),
                    sigma_str: _round_float(mysigma, 4),
                    sigma_str_rock: _round_float(mysigma_rock, 4),
                    sigma_str_soil: _round_float(mysigma_soil, 4),
                    bias_str: _round_float(mybias, 4),
                })
            #
            # For df1 stations, add the MMIs comverted from PGM
            #
            if ndf == 'df1':
                station['properties']['mmi_from_pgm'] = []
                for myimt in getattr(self, ndf).imts:
                    if myimt == 'MMI':
                        continue
                    dimtstr = 'derived_MMI_from_' + myimt
                    if dimtstr not in sdf:
                        continue
                    imt_name = myimt.lower()
                    myamp = sdf[dimtstr][six]
                    mysd = sdf[dimtstr + '_sd'][six]
                    if np.isnan(myamp):
                        myamp = 'null'
                        mysd = 'null'
                    station['properties']['mmi_from_pgm'].append({
                        'name': imt_name,
                        'value': _round_float(myamp, 2),
                        'sigma': _round_float(mysd, 2),
                    })

            #
            # For df2 stations, add the PGMs converted from MMI
            #
            if ndf == 'df2':
                station['properties']['pgm_from_mmi'] = []
                for myimt in getattr(self, ndf).imts:
                    if myimt == 'MMI':
                        continue
                    imt_name = myimt.lower()
                    myamp = sdf[myimt][six]
                    mysd = sdf[myimt + '_sd'][six]
                    if myimt == 'PGV':
                        value = np.exp(myamp)
                        units = 'cm/s'
                    else:
                        value = np.exp(myamp) * 100
                        units = '%g'
                    if np.isnan(value):
                        value = 'null'
                        mysd = 'null'
                    station['properties']['pgm_from_mmi'].append({
                        'name': imt_name,
                        'value': _round_float(value, 4),
                        'units': units,
                        'ln_sigma': _round_float(mysd, 4),
                    })
            #
            # Set the generic distance property (this is rrup)
            #
            station['properties']['distance'] = \
                _round_float(sdf['rrup'][six], 3)
            station['properties']['distance_stddev'] = \
                _round_float(np.sqrt(sdf['rrup_var'][six]), 3)
            #
            # Set the specific distances properties
            #
            station['properties']['distances'] = {}
            for dm, dm_arr in dist_dict[ndf].items():
                station['properties']['distances'][dm] = \
                    _round_float(dm_arr[six], 3)
            #
            # Set the outlier flags
            #
            mflag = '0'
            if ndf == 'df1' and sdf['flagged'][six]:
                mflag = 'ManuallyFlagged'
            for channel in station['properties']['channels']:
                for amp in channel['amplitudes']:
                    if amp['flag'] != '0':
                        amp['flag'] += ',' + mflag
                    else:
                        amp['flag'] = mflag
                    Name = amp['name'].upper()
                    if sdf[Name + '_outliers'][six]:
                        if amp['flag'] == '0':
                            amp['flag'] = 'Outlier'
                        else:
                            amp['flag'] += ',Outlier'
        sjdict['metadata'] = {'eventid': self._eventid}
        return sjdict

    def _storeGriddedData(self, oc):
        """
        Store gridded data in the output container.
        """
        metadata = {}
        metadata['xmin'] = self.W
        metadata['xmax'] = self.E
        metadata['ymin'] = self.S
        metadata['ymax'] = self.N
        metadata['nx'] = self.smnx
        metadata['ny'] = self.smny
        metadata['dx'] = self.smdx
        metadata['dy'] = self.smdy
        #
        # Put the Vs30 grid in the output container
        #
        _, units, digits = _get_layer_info('vs30')
        metadata['units'] = units
        metadata['digits'] = digits
        oc.setArray([], 'vs30', self.sx_out.vs30, metadata=metadata)
        #
        # Now do the distance grids
        #
        metadata['units'] = 'km'
        metadata['digits'] = 4
        for dm in get_distance_measures():
            dm_arr = getattr(self.dx_out, dm, None)
            if dm_arr is None:
                continue
            oc.setArray(['distances'], dm, dm_arr, metadata=metadata)
            if dm in ('rrup', 'rjb'):
                dm_var = getattr(self.dx_out, dm + '_var', None)
                if dm_var is None:
                    dm_var = np.zeros_like(dm_arr)
                dm_std = np.sqrt(dm_var)
                oc.setArray(['distances'], dm + '_std', dm_std,
                            metadata=metadata)

        #
        # Output the data and uncertainty grids
        #
        component = self.config['interp']['component']
        std_metadata = copy.copy(metadata)
        for key, value in self.outgrid.items():
            # set the data grid
            _, units, digits = _get_layer_info(key)
            metadata['units'] = units
            metadata['digits'] = digits

            # set the uncertainty grid
            std_layername, units, digits = _get_layer_info(key + '_sd')
            std_metadata['units'] = units
            std_metadata['digits'] = digits
            oc.setIMTGrids(key, value, metadata, self.outsd[key],
                           std_metadata, component)
        #
        # Directivity
        #
        del metadata['units']
        del metadata['digits']
        if self.do_directivity is True:
            for k, v in self.dir_output.items():
                imtstr, _, _ = _get_layer_info(k)
                oc.setArray(['directivity'], imtstr, v, metadata=metadata)

    def _storePointData(self, oc):
        """
        Store point data in the output container.
        """
        #
        # Store the Vs30
        #
        vs30_metadata = {
            'units': 'm/s',
            'digits': 4
        }
        oc.setArray([], 'vs30', self.sx_out.vs30.ravel(),
                    metadata=vs30_metadata)
        #
        # Store the distances
        #
        distance_metadata = {
            'units': 'km',
            'digits': 4
        }
        for dm in get_distance_measures():
            dm_arr = getattr(self.dx_out, dm, None)
            if dm_arr is None:
                continue
            oc.setArray(['distances'], dm, dm_arr.ravel(),
                        metadata=distance_metadata)
            if dm in ('rrup', 'rjb'):
                dm_var = getattr(self.dx_out, dm + '_var', None)
                if dm_var is None:
                    dm_var = np.zeros_like(dm_arr)
                oc.setArray(['distances'], dm + '_std',
                            np.sqrt(dm_var).ravel(),
                            metadata=distance_metadata)
        #
        # Store the IMTs
        #
        ascii_ids = np.array(
            [x.encode('ascii') for x in self.idents]).ravel()
        component = self.config['interp']['component']
        for key, value in self.outgrid.items():
            # set the data grid
            _, units, digits = _get_layer_info(key)
            mean_metadata = {
                'units': units,
                'digits': digits
            }
            # set the uncertainty grid
            std_layername, units, digits = _get_layer_info(key + '_sd')
            std_metadata = {
                'units': units,
                'digits': digits
            }
            oc.setIMTArrays(key,
                            self.sx_out.lons.ravel(),
                            self.sx_out.lats.ravel(),
                            ascii_ids,
                            value.ravel(), mean_metadata,
                            self.outsd[key].ravel(), std_metadata,
                            component)

    def _storeAttenuationData(self, oc):
        """
        Output arrays of rock and soil attenuation curves
        """

        for dist_type in ['repi', 'rhypo', 'rrup', 'rjb']:
            oc.setArray(['attenuation', 'distances'], dist_type,
                        getattr(self.atten_dx, dist_type, None))

        imtstrs = self.atten_rock_mean.keys()
        for imtstr in imtstrs:
            oc.setArray(['attenuation', 'rock', imtstr], 'mean',
                        self.atten_rock_mean[imtstr])
            oc.setArray(['attenuation', 'soil', imtstr], 'mean',
                        self.atten_soil_mean[imtstr])
            oc.setArray(['attenuation', 'rock', imtstr], 'std',
                        self.atten_rock_sd[imtstr])
            oc.setArray(['attenuation', 'soil', imtstr], 'std',
                        self.atten_soil_sd[imtstr])
        return

    #
    # Helper function to call get_mean_and_stddevs for the
    # appropriate object given the IMT and describe the
    # MultiGMPE structure.
    #
    def _gmas(self, gmpe, sx, dx, oqimt, apply_gafs):
        """
        This is a helper function to call get_mean_and_stddevs for the
        appropriate object given the IMT.

        Args:
            gmpe:
                A GMPE instance.
            sx:
                Sites context.
            dx:
                Distance context.
            oqimt:
                List of OpenQuake IMTs.
            apply_gafs (boolean):
                Whether or not to apply the generic
                amplification factors to the GMPE output.

        Returns:
            tuple: Tuple of two items:

                - Numpy array of means,
                - List of numpy array of standard deviations corresponding to
                  therequested stddev_types.

        """
        if 'MMI' in oqimt:
            pe = self.ipe
            sd_types = self.ipe_stddev_types
        else:
            pe = gmpe
            sd_types = self.gmpe_stddev_types

            # --------------------------------------------------------------------
            # Describe the MultiGMPE
            # --------------------------------------------------------------------
            if not hasattr(self, '_info'):
                self._info = {
                    'multigmpe': {}
                }
            self._info['multigmpe'][str(oqimt)] = gmpe.describe()

        mean, stddevs = pe.get_mean_and_stddevs(
            copy.deepcopy(sx), self.rx,
            copy.deepcopy(dx), oqimt,
            sd_types)

        # Include generic amp factors?
        if apply_gafs:
            gafs = get_generic_amp_factors(sx, str(oqimt))
            if gafs is not None:
                mean += gafs

        # Does directivity apply to this imt?
        row_pers = Rowshandel2013.getPeriods()

        if isinstance(oqimt, imt.PGA):
            imt_ok = False
        elif isinstance(oqimt, imt.PGV) or isinstance(oqimt, imt.MMI):
            tper = 1.0
            imt_ok = True
        elif isinstance(oqimt, imt.SA):
            tper = oqimt.period
            min_per = np.min(row_pers)
            max_per = np.max(row_pers)
            if (tper >= min_per) and (tper <= max_per):
                imt_ok = True
            else:
                imt_ok = False
        else:
            imt_ok = False

        # Did we calculate directivity?
        calc_dir = self.do_directivity

        if calc_dir and imt_ok:
            # Use distance context to figure out which directivity result
            # we need to use.
            all_fd = None
            for dirdf, tmpdx in self.dir_results:
                if dx == tmpdx:
                    all_fd = dirdf.getFd()
                    break
            if all_fd is None:
                raise RuntimeError(
                    "Failed to detect dataframe for directivity calculation.")

            # Does oqimt match any of those periods?
            if tper in row_pers:
                fd = all_fd[row_pers.index(tper)]
            else:
                # Log(period) interpolation.
                apers = np.array(row_pers)
                per_below = np.max(apers[apers < tper])
                per_above = np.min(apers[apers > tper])
                fd_below = all_fd[row_pers.index(per_below)]
                fd_above = all_fd[row_pers.index(per_above)]
                x1 = np.log(per_below)
                x2 = np.log(per_above)
                fd = fd_below + (np.log(tper) - x1) * \
                    (fd_above - fd_below)/(x2 - x1)
            # Reshape to match the mean
            fd = fd.reshape(mean.shape)
            # Store the interpolated grid
            imtstr = str(oqimt)
            self.dir_output[imtstr] = fd
            if isinstance(oqimt, imt.MMI):
                mean *= np.exp(fd)
            else:
                mean += fd

        return mean, stddevs

    def _adjustResolution(self):
        """
        This is a helper function to adjust the resolution to be under
        the maximum value specified in the config.
        """
        # We want to only use resolutions that are multiples of 1 minute or
        # an integer division of 1 minute.
        one_minute = 1/60
        multiples = np.arange(1, 11)
        divisions = 1/multiples
        factors = np.sort(np.unique(np.concatenate((divisions, multiples))))
        ok_res = one_minute * factors
        latspan = self.N - self.S

        # Deal with possible 180 longitude disontinuity
        if self.E > self.W:
            lonspan = self.E - self.W
        else:
            xmax = self.E + 360
            lonspan = xmax - self.W
        nx = np.floor(lonspan / self.smdx) + 1
        ny = np.floor(latspan / self.smdy) + 1
        ngrid = nx * ny
        nmax = self.nmax
        if ngrid > nmax:
            self.logger.info('Extent and resolution of shakemap results in '
                             'too many grid points. Adjusting resolution...')
            self.logger.info('Longitude span: %f' % lonspan)
            self.logger.info('Latitude span: %f' % latspan)
            self.logger.info('Current dx: %f' % self.smdx)
            self.logger.info('Current dy: %f' % self.smdy)
            self.logger.info('Current number of grid points: %i' % ngrid)
            self.logger.info('Max grid points allowed: %i' % nmax)
            target_res = \
                (-(latspan + lonspan) -
                 np.sqrt(latspan**2 + lonspan**2 +
                         2 * latspan * lonspan *
                         (2 * nmax - 1))) / (2 * (1 - nmax))

            if np.any(ok_res > target_res):
                sel_res = np.min(ok_res[ok_res > target_res])
            else:
                sel_res = np.max(ok_res)
            self.smdx = sel_res
            self.smdy = sel_res
            self.logger.info('Updated dx: %f' % self.smdx)
            self.logger.info('Updatd dy: %f' % self.smdy)
            nx = np.floor(lonspan / self.smdx) + 1
            ny = np.floor(latspan / self.smdy) + 1
            self.logger.info('Updated number of grid points: %i'
                             % (nx * ny))


def _round_float(val, digits):
    if ma.is_masked(val) or val == '--' or val == 'null' or np.isnan(val):
        return None
    return float(('%.' + str(digits) + 'f') % val)


def _string_round(val, digits):
    if ma.is_masked(val) or val == '--' or val == 'null' or np.isnan(val):
        return None
    return str(_round_float(val, digits))


def _get_period_arrays(*args):
    """
    Return 1) a sorted array of the periods represented by the IMT list(s)
    in the input, and 2) a dictionary of the IMTs and their indices.

    Args:
        *args (list): One or more lists of IMTs.

    Returns:
        array, dict: Numpy array of the (sorted) periods represented by the
        IMTs in the input list(s), and a dictionary of the IMTs and their
        indices into the period array.
    """
    imt_per = set()
    imt_per_ix = {}
    for imt_list in args:
        if imt_list is None:
            continue
        for imtstr in imt_list:
            if imtstr == 'PGA':
                period = 0.01
            elif imtstr in ('PGV', 'MMI'):
                period = 1.0
            else:
                period = float(imtstr.replace('SA(', '').replace(')', ''))
            imt_per.add(period)
            imt_per_ix[imtstr] = period
    imt_per = sorted(imt_per)
    for imtstr, period in imt_per_ix.items():
        imt_per_ix[imtstr] = imt_per.index(period)
    return np.array(imt_per), imt_per_ix


def _get_period_from_imt(imtstr):
    """
    Return a float representing the period of the SA IMT in the input.

    Args:
        imtstr (str): A string representing an SA IMT.

    Returns:
        float: The period of the SA IMT as a float.
    """
    return float(imtstr.replace('SA(', '').replace(')', ''))


def _get_nearest_imts(imtstr, imtset, saset):
    """
    Return the input IMT, or it's closest surrogarte (or bracket) found
    in imtset.

    Args:
        imtstr (str): An (OQ-style) IMT string.
        imtset (list): A list of IMTs to search for imtstr or its closest
            surrogate (or bracket).
        saset (list): The SA IMTs found in imtset.

    Returns:
        tuple: The IMT, it's closest surrogate, or a bracket of SAs with
        periods on either side of the IMT's period, from the IMTs in intset.
    """
    if imtstr in imtset:
        return (imtstr, )
    #
    # If we're here, then we know that IMT isn't in the inputs. Try
    # some alternatives.
    #
    if imtstr == 'PGA':
        #
        # Use the highest frequency in the inputs, otherwise use PGV
        #
        if len(saset):
            return (sorted(saset, key=_get_period_from_imt)[0], )
        elif 'PGV' in imtset:
            return ('PGV', )
        else:
            return ()
    elif imtstr == 'PGV':
        #
        # Use 1.0 sec SA (or its bracket) if it's there, otherwise
        # use PGA
        #
        if 'SA(1.0)' in saset:
            return ('SA(1.0)', )
        sa_tuple = _get_sa_bracket('SA(1.0)', saset)
        if sa_tuple != ():
            return sa_tuple
        if 'PGA' in imtset:
            return ('PGA', )
        else:
            return ()
    elif imtstr == 'MMI':
        #
        # Use PGV if it's there, otherwise use 1.0 sec SA (or its
        # bracket)
        #
        if 'PGV' in imtset:
            return ('PGV', )
        if 'SA(1.0)' in saset:
            return ('SA(1.0)', )
        return _get_sa_bracket('SA(1.0)', saset)
    elif imtstr.startswith('SA('):
        #
        # We know the actual IMT isn't here, so get the bracket
        #
        return _get_sa_bracket(imtstr, saset)
    else:
        raise ValueError('Unknown IMT %s in get_imt_bracket' % imtstr)


def _get_sa_bracket(myimt, saset):
    """
    For a given SA IMT, look through the input SAs and return a tuple of
    a) a pair of IMT strings representing the periods bracketing the given
    period; or c) the single IMT representing the first or last period in
    the input list if the given period is off the end of the list.

    Args:
        myper (float): The period to search for in the input lists.
        saset (list): A list of SA IMTs.

    Returns:
        tuple: One or two strings representing the IMTs closest to or
        bracketing the input IMT.

    """
    if not len(saset):
        return ()
    #
    # Stick the target IMT into a copy of the list of SAs, then sort
    # the list by period.
    #
    ss = saset.copy()
    ss.append(myimt)
    tmplist = sorted(ss, key=_get_period_from_imt)
    nimt = len(tmplist)
    #
    # Get the index of the target IMT in the sorted list
    #
    myix = tmplist.index(myimt)
    #
    # If the target IMT is off the end of the list, return the
    # appropriate endpoint; else return the pair of IMTs that
    # bracket the target.
    #
    if myix == 0:
        return (tmplist[1], )
    elif myix == nimt - 1:
        return (tmplist[-2], )
    else:
        return (tmplist[myix - 1], tmplist[myix + 1])


def _get_imt_lists(df):
    """
    Given a data frame, return a list of lists of valid IMTS for
    each station in the dataframe; also return a list of the valid
    SA IMTs for each station.

    Args:
        df (DataFrame): A DataFrame.

    Returns:
        list, list: Two lists of lists: each list contains lists
        corresponding to the stations in the data frame: the first
        list contains all of the valid IMTs for that station, the
        second list contains just the valid SA IMTs for the station.
    """
    imtlist = []
    salist = []
    nlist = np.size(df.df['lon'])
    for ix in range(nlist):
        valid_imts = []
        sa_imts = []
        if 'flagged' not in df.df or not df.df['flagged'][ix]:
            for this_imt in df.imts:
                if not np.isnan(df.df[this_imt + '_residual'][ix]) and \
                        not df.df[this_imt + '_outliers'][ix]:
                    valid_imts.append(this_imt)
                    if this_imt.startswith('SA('):
                        sa_imts.append(this_imt)
        imtlist.append(valid_imts)
        salist.append(sa_imts)
    return imtlist, salist


def _get_map_grade(do_grid, outsd, psd, moutgrid):
    """
    Computes a 'grade' for the map. Essentially looks at the ratio of
    the computed PGA uncertainty to the predicted PGA uncertainty for
    the area inside the MMI 6 contour. If the maximum MMI is less than
    6, or the map is not a grid, the grade and mean ratio are set to '--'.

    Args:
        do_grid (bool): Is the map a grid (True) or a list of points
            (False)?

        outsd (dict): A dictionary of computed uncertainty arrays.

        psd (dict): A dictionary of predicted uncertainty arrays.

        moutgrid (dict): A dictionary of landmasked output ground
            motion arrays.

    Returns:
        tuple: The mean uncertainty ratio and the letter grade.
    """
    mean_rat = '--'
    mygrade = '--'
    if not do_grid or 'PGA' not in outsd or 'PGA' not in psd \
            or 'MMI' not in moutgrid:
        return mean_rat, mygrade
    sd_rat = outsd['PGA'] / psd['PGA']
    mmimasked = ma.masked_less(moutgrid['MMI'], 6.0)
    mpgasd_rat = ma.masked_array(sd_rat, mask=mmimasked.mask)
    if not np.all(mpgasd_rat.mask):
        gvals = [0.96, 0.98, 1.05, 1.25]
        grades = ['A', 'B', 'C', 'D', 'F']
        mean_rat = mpgasd_rat.mean()
        for ix, val in enumerate(gvals):
            if mean_rat <= val:
                mygrade = grades[ix]
                break
        if mygrade == '--':
            mygrade = 'F'
    return mean_rat, mygrade


# we need a way to get units information about intensity measure types
# and translate between openquake naming convention and ShakeMap grid naming
# convention.
def _get_layer_info(layer):
    """
    We need a way to get units information about intensity measure types
    and translate between OpenQuake naming convention and ShakeMap grid naming
    convention.

    Args:
        layer (str): ShakeMap grid name.

    Returns:
        tuple: Tuple including:

            - OpenQuake naming convention,
            - units,
            - significant digits.

    """
    layer_out = layer
    layer_units = ''
    layer_digits = 4  # number of significant digits

    if layer.endswith('_sd'):
        layer_out = oq_to_file(layer.replace('_sd', ''))
        layer_out = layer_out + '_sd'
    else:
        layer_out = oq_to_file(layer)
    if layer.startswith('SA'):
        layer_units = 'ln(g)'
    elif layer.startswith('PGA'):
        layer_units = 'ln(g)'
    elif layer.startswith('PGV'):
        layer_units = 'ln(cm/s)'
    elif layer.startswith('MMI'):
        layer_units = 'intensity'
        layer_digits = 2
    elif layer.startswith('vs30'):
        layer_units = 'm/s'
    else:
        raise ValueError('Unknown layer type: %s' % layer)

    return (layer_out, layer_units, layer_digits)
