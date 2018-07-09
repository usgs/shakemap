"""
Process a ShakeMap, based on the configuration and data found in
shake_data.hdf, and produce output in shake_result.hdf.
"""
import warnings
import os.path
import time as time
import copy
from time import gmtime, strftime
import shutil
from collections import OrderedDict

import numpy as np
import numpy.ma as ma
import numexpr as ne
from mpl_toolkits.basemap import maskoceans
from openquake.hazardlib import imt
import openquake.hazardlib.const as oqconst

import concurrent.futures as cf

# local imports
from .base import CoreModule
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
from shakelib.utils.containers import ShakeMapOutputContainer
from shakelib.utils.distance import geodetic_distance_fast

from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D

from shakemap.utils.config import get_config_paths
from shakemap.utils.utils import get_object_from_config
from shakemap._version import get_versions
from shakemap.utils.generic_amp import get_generic_amp_factors

#
# TODO: Some constants; these should maybe be in a configuration
# or in a constants module or be GMICE-specific.
#
# default_mmi_stddev: the standard deviation of MMI values if it
#                     is not specified in the input
# min_mmi_convert: the minimum MMI to convert to PGM -- low
#                  intensities don't convert very accurately
# default_stddev_inter: This is a stand-in for tau when the gmpe set
#                       doesn't provide it. It is an educated guess
#                       based on the NGA-east work and BC Hydro gmpe.
#                       It's not perfect, but probably isn't too far off.
#                       It is only used when the GMPEs don't provide a
#                       breakdown of the uncertainty terms.
#
SM_CONSTS = {'default_mmi_stddev': 0.3,
             'min_mmi_convert': 4.0,
             'default_stddev_inter': 0.35}

TIMEFMT = '%Y-%m-%dT%H:%M:%SZ'


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

    # supply here a data structure with information about files that
    # can be created by this module.
    contents = OrderedDict.fromkeys(['shakemapHDF'])
    contents['shakemapHDF'] = {
        'title': 'Comprehensive ShakeMap HDF Data File',
        'caption': 'HDF file containing all ShakeMap results.',
        'formats': [{'filename': 'shake_result.hdf',
                     'type': 'application/x-bag'}]
    }

    def execute(self):
        """
        Interpolate ground motions to a grid or list of locations.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_data HDF file does not exist.
        """
        self.logger.debug('Starting model...')
        # ------------------------------------------------------------------
        # Make the input container and extract the config
        # ------------------------------------------------------------------
        self._setInputContainer()
        self.config = self.ic.getConfig()
        # ------------------------------------------------------------------
        # Clear away results from previous runs
        # ------------------------------------------------------------------
        self._clearProducts()
        # ------------------------------------------------------------------
        # Retrieve a bunch of config options and set them as attributes
        # ------------------------------------------------------------------
        self._setConfigOptions()
        # ------------------------------------------------------------------
        # Instantiate the gmpe, gmice, and ipe
        # Here we make a placeholder gmpe so that we can make the
        # rupture and distance contexts; later we'll make the
        # IMT-specific gmpes
        # ------------------------------------------------------------------
        self.default_gmpe = MultiGMPE.from_config(self.config)

        self.gmice = get_object_from_config('gmice', 'modeling', self.config)

        if self.config['ipe_modules'][self.config['modeling']['ipe']][0] == \
                'VirtualIPE':
            pgv_imt = imt.from_string('PGV')
            ipe_gmpe = MultiGMPE.from_config(self.config, filter_imt=pgv_imt)
            self.ipe = VirtualIPE.fromFuncs(ipe_gmpe, self.gmice)
        else:
            self.ipe = get_object_from_config('ipe', 'modeling', self.config)
        # ------------------------------------------------------------------
        # Get the rupture object and rupture context
        # ------------------------------------------------------------------
        self.rupture_obj = self.ic.getRuptureObject()
        if 'mechanism' in self.config['modeling']:
            self.rupture_obj._origin.setMechanism(
                mech=self.config['modeling']['mechanism'])
        self.rx = self.rupture_obj.getRuptureContext([self.default_gmpe])
        # TODO: figure out how to not have to do this
        if self.rx.rake is None:
            self.rx.rake = 0
        # ------------------------------------------------------------------
        # The output locations: either a grid or a list of points
        # ------------------------------------------------------------------
        self.logger.debug('Setting output params...')
        self._setOutputParams()
        # ------------------------------------------------------------------
        # If the gmpe doesn't break down its stardard deviation into
        # within- and between-event terms, we need to handle things
        # somewhat differently.
        # ------------------------------------------------------------------
        gmpe_sd_types = self.default_gmpe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if len(gmpe_sd_types) == 1:
            self.total_sd_only = True
            self.stddev_types = [oqconst.StdDev.TOTAL]
        else:
            self.total_sd_only = False
            self.stddev_types = [oqconst.StdDev.TOTAL,
                                 oqconst.StdDev.INTER_EVENT,
                                 oqconst.StdDev.INTRA_EVENT]
        # ------------------------------------------------------------------
        # Station data: Create DataFrame(s) with the input data:
        # df1 for instrumented data
        # df2 for non-instrumented data
        # ------------------------------------------------------------------
        self.logger.debug('Setting data frames...')
        self._setDataFrames()
        # ------------------------------------------------------------------
        # Add the predictions, etc. to the data frames
        # ------------------------------------------------------------------
        self.logger.debug('Populating data frames...')
        self._populateDataFrames()
        # ------------------------------------------------------------------
        # Try to make all the derived IMTs possible from MMI (if we have MMI)
        # ------------------------------------------------------------------
        self._deriveIMTsFromMMI()
        # ------------------------------------------------------------------
        # ------------------------------------------------------------------
        self._deriveMMIFromIMTs()

        self.logger.debug('Getting combined IMTs')
        # ------------------------------------------------------------------
        # Get the combined set of input and output IMTs, their periods,
        # and an index dictionary, then make the cross-correlation function
        # ------------------------------------------------------------------
        self.combined_imt_set = self.imt_out_set.copy()
        for ndf in self.dataframes:
            self.combined_imt_set |= getattr(self, ndf).imts

        self.imt_per, self.imt_per_ix = _get_period_arrays(
            self.combined_imt_set)
        self.ccf = get_object_from_config('ccf', 'modeling',
                                          self.config, self.imt_per)

        self.logger.debug('Doing bias')
        # ------------------------------------------------------------------
        # Do the bias for all of the input and output IMTs. Hold on
        # to some of the products that will be used for the interpolation.
        # ------------------------------------------------------------------
        self.nominal_bias = {}
        self.bias_num = {}   # The numerator term of the bias
        self.bias_den = {}   # The denominator of the bias (w/o the tau term)
        self.psd = {}        # phi of the output points
        self.tsd = {}        # tau of the output points

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

        # ------------------------------------------------------------------
        # Now do the MVN with the intra-event residuals
        # ------------------------------------------------------------------
        self.outgrid = {}   # Holds the interpolated output arrays keyed by IMT
        self.outsd = {}     # Holds the standard deviation arrays keyed by IMT

        #
        # Places to put the results for the regression plots
        #
        self.rockgrid = {}
        self.soilgrid = {}
        self.rocksd = {}
        self.soilsd = {}

        self.logger.debug('Doing MVN...')
        with cf.ThreadPoolExecutor(max_workers=self.max_workers) as ex:
            ex.map(self._computeMVN, self.imt_out_set)
#        self._computeMVN()

        # ------------------------------------------------------------------
        # Output the data and metadata
        # ------------------------------------------------------------------
        product_path = os.path.join(self.datadir, 'products')
        if not os.path.isdir(product_path):
            os.mkdir(product_path)
        oc = ShakeMapOutputContainer.create(os.path.join(
            product_path, 'shake_result.hdf'))
        # ------------------------------------------------------------------
        # Might as well stick the whole config in the result
        # ------------------------------------------------------------------
        oc.setConfig(self.config)

        #
        # We're going to need masked arrays of the output grids later, so
        # make them now.
        #
        moutgrid = self._getMaskedGrids()

        #
        # Get the info dictionary that will become info.json, and
        # store it in the output container
        #
        info = self._getInfo(moutgrid)
        oc.setDictionary('info.json', info)

        # ------------------------------------------------------------------
        # Add the rupture JSON as a text string
        # ------------------------------------------------------------------
        oc.setRupture(self.rupture_obj)

        # ------------------------------------------------------------------
        # Fill the station dictionary for stationlist.json and add it to
        # the output container
        # ------------------------------------------------------------------
        sjdict = self._fillStationJSON()
        oc.setStationDict(sjdict)

        # ------------------------------------------------------------------
        # Add the output grids or points to the output; include some
        # metadata.
        # ------------------------------------------------------------------
        if self.do_grid:
            self._storeGriddedData(oc)
        else:
            self._storePointData(oc)

        if self.do_grid:
            self._storeRegressionData(oc)

        oc.close()
        self.ic.close()
    # ------------------------------------------------------------------
    # End execute()
    # ------------------------------------------------------------------

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
        # ------------------------------------------------------------------
        # Processing parameters
        # ------------------------------------------------------------------
        self.max_workers = self.config['system']['max_workers']
        # ------------------------------------------------------------------
        # Do we apply the generic amplification factors?
        # ------------------------------------------------------------------
        self.apply_gafs = self.config['modeling']['apply_generic_amp_factors']
        # ------------------------------------------------------------------
        # Bias parameters
        # ------------------------------------------------------------------
        self.do_bias = self.config['modeling']['bias']['do_bias']
        self.bias_max_range = self.config['modeling']['bias']['max_range']
        self.bias_max_mag = self.config['modeling']['bias']['max_mag']
        self.bias_max_dsigma = \
            self.config['modeling']['bias']['max_delta_sigma']
        # ------------------------------------------------------------------
        # Outlier parameters
        # ------------------------------------------------------------------
        self.do_outliers = self.config['data']['outlier']['do_outliers']
        self.outlier_deviation_level = \
            self.config['data']['outlier']['max_deviation']
        self.outlier_max_mag = self.config['data']['outlier']['max_mag']
        # ------------------------------------------------------------------
        # These are the IMTs we want to make
        # ------------------------------------------------------------------
        self.imt_out_set = set(self.config['interp']['imt_list'])
        # ------------------------------------------------------------------
        # The x and y resolution of the output grid
        # ------------------------------------------------------------------
        self.smdx = self.config['interp']['prediction_location']['xres']
        self.smdy = self.config['interp']['prediction_location']['yres']
        # ------------------------------------------------------------------
        # Get the Vs30 file name
        # ------------------------------------------------------------------
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
        if self.config['interp']['prediction_location']['file'] and \
           self.config['interp']['prediction_location']['file'] != 'None':
            #
            # FILE: Open the file and get the output points
            #
            self.do_grid = False
            in_sites = np.genfromtxt(
                self.config['interp']['prediction_location']['file'],
                autostrip=True, unpack=True,
                dtype=[np.float, np.float, np.float, '<U80'])
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
            dist_obj_out = Distance(self.default_gmpe, self.lons, self.lats,
                                    self.depths, self.rupture_obj)

            self.sites_obj_out = Sites.fromBounds(self.W, self.E, self.S,
                                                  self.N, self.smdx, self.smdy,
                                                  defaultVs30=self.vs30default,
                                                  vs30File=self.vs30_file,
                                                  padding=True, resample=True)

            self.sx_out = self.sites_obj_out.getSitesContext(
                {'lats': self.lats,
                 'lons': self.lons})
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
                self.W, self.E, self.S, self.N = get_extent(self.rupture_obj,
                                                            config=self.config)

            self.sites_obj_out = Sites.fromBounds(self.W, self.E, self.S,
                                                  self.N, self.smdx, self.smdy,
                                                  defaultVs30=self.vs30default,
                                                  vs30File=self.vs30_file,
                                                  padding=True, resample=True)
            self.smnx, self.smny = self.sites_obj_out.getNxNy()
            self.sx_out = self.sites_obj_out.getSitesContext()
            #
            # Grids on rock and soil for the regression plots
            #
            self.sx_rock = self.sites_obj_out.getSitesContext(rock_vs30=760)
            self.sx_soil = self.sites_obj_out.getSitesContext(rock_vs30=180)
            lons, lats = np.meshgrid(self.sx_out.lons,
                                     self.sx_out.lats)
            self.sx_out.lons = np.flipud(lons.copy())
            self.sx_out.lats = np.flipud(lats.copy())
            self.lons = np.flipud(lons).flatten()
            self.lats = np.flipud(lats).flatten()
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

        self.lons_out_rad = np.radians(self.lons)
        self.lats_out_rad = np.radians(self.lats)

    def _setDataFrames(self):
        """
        Extract the StationList object from the input container and
        fill the DataFrame class and keep a list of dataframes.

            - df1 holds the instrumented data (PGA, PGV, SA)
            - df2 holds the non-instrumented data (MMI)
        """
        self.dataframes = []
        self.stations = self.ic.getStationList()
        if self.stations is None:
            return
        for dfid, val in (('df1', True), ('df2', False)):
            sdf, imts = self.stations.getStationDictionary(instrumented=val)
            if sdf is not None:
                df = DataFrame()
                df.df = sdf
                df.imts = imts
                setattr(self, dfid, df)
                self.dataframes.append(dfid)

    def _populateDataFrames(self):
        """
        Make the sites and distance contexts for each dataframe then
        compute the predictions for the IMTs in that dataframe.
        """
        for dfid in self.dataframes:
            dfn = getattr(self, dfid)
            df = dfn.df
            #
            # Get the sites and distance contexts
            #
            df['depth'] = np.zeros_like(df['lon'])
            lldict = {'lons': df['lon'], 'lats': df['lat']}
            dfn.sx = self.sites_obj_out.getSitesContext(lldict)
            dist_obj = Distance(self.default_gmpe, df['lon'], df['lat'],
                                df['depth'], self.rupture_obj)
            dfn.dx = dist_obj.getDistanceContext()
            #
            # Do the predictions and other bookkeeping for each IMT
            #
            for imtstr in dfn.imts:
                oqimt = imt.from_string(imtstr)
                gmpe = None
                if imtstr != 'MMI':
                    gmpe = MultiGMPE.from_config(self.config, filter_imt=oqimt)
                pmean, pstddev = self._gmas(
                    gmpe, dfn.sx, dfn.dx, oqimt, self.apply_gafs)
                df[imtstr + '_pred'] = pmean
                df[imtstr + '_pred_sigma'] = pstddev[0]
                if self.total_sd_only:
                    tau_guess = SM_CONSTS['default_stddev_inter']
                    df[imtstr + '_pred_tau'] = np.full_like(
                        df[imtstr + '_pred'], tau_guess)
                    df[imtstr + '_pred_phi'] = np.sqrt(
                        pstddev[0]**2 - tau_guess**2)
                else:
                    df[imtstr + '_pred_tau'] = pstddev[1]
                    df[imtstr + '_pred_phi'] = pstddev[2]
                #
                # Compute the total residual
                #
                df[imtstr + '_residual'] = \
                    df[imtstr] - df[imtstr + '_pred']
                # ----------------------------------------------------------
                # Do the outlier flagging if we have a fault, or we don't
                # have a fault but the event magnitude is under the limit
                # ----------------------------------------------------------
                if self.do_outliers and \
                        (not isinstance(self.rupture_obj, PointRupture) or
                         self.rx.mag <= self.outlier_max_mag):
                    #
                    # turn off nan warnings for this statement
                    #
                    np.seterr(invalid='ignore')
                    flagged = \
                        np.abs(df[imtstr + '_residual']) > \
                        self.outlier_deviation_level * \
                        df[imtstr + '_pred_sigma']
                    np.seterr(invalid='warn')

                    self.logger.debug('IMT: %s, flagged: %d' %
                                      (imtstr, np.sum(flagged)))
                    df[imtstr + '_outliers'] = flagged
                else:
                    df[imtstr + '_outliers'] = np.full(df[imtstr].shape,
                                                       False,
                                                       dtype=np.bool)
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
            #
            # It will be handy later on to have the rupture distance
            # in the dataframes
            #
            dd = get_distance(['rrup'], df['lat'], df['lon'],
                              df['depth'], self.rupture_obj)
            df['rrup'] = dd['rrup']

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
                oqimt = gmice_imt(period)
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
                gmpe = MultiGMPE.from_config(self.config, filter_imt=oqimt)
                pmean, pstddev = self._gmas(
                    gmpe, self.df2.sx, self.df2.dx, oqimt, self.apply_gafs)
                df2[imtstr + '_pred'] = pmean
                df2[imtstr + '_pred_sigma'] = pstddev[0]
                if self.total_sd_only:
                    tau_guess = SM_CONSTS['default_stddev_inter']
                    df2[imtstr + '_pred_tau'] = np.full_like(
                        df2[imtstr + '_pred'], tau_guess)
                    df2[imtstr + '_pred_phi'] = np.sqrt(pstddev[0]**2 -
                                                        tau_guess**2)
                else:
                    df2[imtstr + '_pred_tau'] = pstddev[1]
                    df2[imtstr + '_pred_phi'] = pstddev[2]
                df2[imtstr + '_residual'] = df2[imtstr] - pmean
                df2[imtstr + '_outliers'] = np.full(
                    pmean.shape, False, dtype=np.bool)

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
        df1['MMI' + '_pred'] = pmean
        df1['MMI' + '_pred_sigma'] = pstddev[0]
        if self.total_sd_only:
            tau_guess = SM_CONSTS['default_stddev_inter']
            df1['MMI' + '_pred_tau'] = np.full_like(
                df1['MMI' + '_pred'], tau_guess)
            df1['MMI' + '_pred_phi'] = np.sqrt(
                pstddev[0]**2 - tau_guess**2)
        else:
            df1['MMI' + '_pred_tau'] = pstddev[1]
            df1['MMI' + '_pred_phi'] = pstddev[2]
        df1['MMI' + '_residual'] = df1['MMI'] - pmean
        df1['MMI' + '_outliers'] = np.full(pmean.shape, False,
                                           dtype=np.bool)

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
            sta_lons_rad_dl = self.sta_lons_rad[imtstr][dix].reshape((-1, 1))
            sta_lats_rad_dl = self.sta_lats_rad[imtstr][dix].reshape((-1, 1))
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
            dist22 = geodetic_distance_fast(sta_lons_rad_dl,
                                            sta_lats_rad_dl,
                                            sta_lons_rad_dl.T,
                                            sta_lats_rad_dl.T)
            d22_rows, d22_cols = np.shape(dist22)  # should be square
            t1_22 = np.tile(sta_period_ix_dl, (1, d22_cols))
            t2_22 = np.tile(sta_period_ix_dl.T, (d22_rows, 1))
            corr22 = self.ccf.getCorrelation(t1_22, t2_22, dist22)
            sigma22 = corr22 * corr_adj22 * (sta_phi_dl * sta_phi_dl.T)
            sigma22inv = np.linalg.pinv(sigma22)
            #
            # Compute the bias numerator and denominator pieces
            #
            if self.do_bias and (not isinstance(self.rupture_obj, PointRupture)
                                 or self.rx.mag <= self.bias_max_mag):
                #
                # Get the correlation between the inputs and outputs
                #
                out_ix_arr = np.full_like(sta_period_ix_dl, outperiod_ix)
                dist_arr = np.zeros_like(out_ix_arr, dtype=float)
                Z = self.ccf.getCorrelation(sta_period_ix_dl, out_ix_arr,
                                            dist_arr)
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

            nom_tau = np.mean(sta_tau_dl.flatten())
            nom_variance = 1.0 / ((1.0 / nom_tau**2) + self.bias_den[imtstr])
            self.nominal_bias[imtstr] = self.bias_num[imtstr] * nom_variance
            bias_time = time.time() - time1
            #
            # Print the nominal values of the bias and its stddev
            #
            self.logger.debug(
                '%s: nom bias %f nom stddev %f; %d stations (time=%f sec)'
                % (imtstr, self.nominal_bias[imtstr], np.sqrt(nom_variance),
                   np.size(sta_lons_rad_dl), bias_time))

    def _computeMVN(self, imtstr):
        """
        Do the MVN computations
        """
        time1 = time.time()
        #
        # Get the index of the (pesudo-) period of the output IMT
        #
        outperiod_ix = self.imt_per_ix[imtstr]
        #
        # Get the predictions at the output points
        #
        oqimt = imt.from_string(imtstr)
        gmpe = None
        if imtstr != 'MMI':
            gmpe = MultiGMPE.from_config(self.config, filter_imt=oqimt)
        pout_mean, pout_sd = self._gmas(
            gmpe, self.sx_out, self.dx_out, oqimt, self.apply_gafs)
        if self.do_grid:
            #
            # Fill the grids for the regression plots
            #
            x_mean, x_sd = self._gmas(
                gmpe, self.sx_rock, self.dx_out, oqimt, False)
            self.rockgrid[imtstr] = x_mean
            self.rocksd[imtstr] = x_sd[0]
            x_mean, x_sd = self._gmas(
                gmpe, self.sx_soil, self.dx_out, oqimt, False)
            self.soilgrid[imtstr] = x_mean
            self.soilsd[imtstr] = x_sd[0]

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
        if self.total_sd_only:
            self.psd[imtstr] = np.sqrt(
                pout_sd[0]**2 - SM_CONSTS['default_stddev_inter']**2)
            self.tsd[imtstr] = np.full_like(
                self.psd[imtstr], SM_CONSTS['default_stddev_inter'])
        else:
            self.psd[imtstr] = pout_sd[2]
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
        #
        # Re-build the covariance matrix of the residuals with
        # the full set of data
        #
        dist22 = geodetic_distance_fast(self.sta_lons_rad[imtstr],
                                        self.sta_lats_rad[imtstr],
                                        self.sta_lons_rad[imtstr].T,
                                        self.sta_lats_rad[imtstr].T)
        d22_rows, d22_cols = np.shape(dist22)  # should be square
        t1_22 = np.tile(self.sta_period_ix[imtstr], (1, d22_cols))
        t2_22 = np.tile(self.sta_period_ix[imtstr].T, (d22_rows, 1))
        corr22 = self.ccf.getCorrelation(t1_22, t2_22, dist22)
        #
        # Rebuild sigma22_inv now that we have updated phi and
        # the correlation adjustment factors
        #
        sigma22 = corr22 * corr_adj22 * \
            (self.sta_phi[imtstr] * self.sta_phi[imtstr].T)
        sigma22inv = np.linalg.pinv(sigma22)
        #
        # Now do the MVN itself...
        #
        dtime = mtime = ddtime = ctime = stime = atime = 0

        ampgrid = np.zeros_like(pout_mean)
        sdgrid = np.zeros_like(pout_mean)
        corr_adj12 = corr_adj * np.ones((1, self.smnx))  # noqa
        for iy in range(self.smny):
            ss = iy * self.smnx
            se = (iy + 1) * self.smnx
            time4 = time.time()
            dist12 = geodetic_distance_fast(
                self.lons_out_rad[ss:se].reshape(1, -1),
                self.lats_out_rad[ss:se].reshape(1, -1),
                self.sta_lons_rad[imtstr],
                self.sta_lats_rad[imtstr])
            t2_12 = np.full(dist12.shape, outperiod_ix, dtype=np.int)
            _, d12_cols = np.shape(dist12)
            t1_12 = np.tile(self.sta_period_ix[imtstr], (1, d12_cols))
            ddtime += time.time() - time4
            time4 = time.time()
            corr12 = self.ccf.getCorrelation(t1_12, t2_12, dist12)  # noqa
            ctime += time.time() - time4
            time4 = time.time()
            # sdarr is the standard deviation of the output sites
            sdarr = self.psd[imtstr][iy, :].reshape((1, -1))  # noqa
            # sdsta is the standard deviation of the stations
            sdsta = self.sta_phi[imtstr]  # noqa
            sigma12 = ne.evaluate(
                "corr12 * corr_adj12 * (sdsta * sdarr)").T
            stime += time.time() - time4
            time4 = time.time()
            #
            # Sigma12 * Sigma22^-1 is known as the 'regression
            # coefficient' matrix (rcmatrix)
            #
            rcmatrix = sigma12.dot(sigma22inv)
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
            sdgrid[iy, :] = \
                pout_sd2[iy, :] - np.sum(rcmatrix * sigma12, axis=1)
            mtime += time.time() - time4

        self.outgrid[imtstr] = ampgrid
        sdgrid[sdgrid < 0] = 0
        self.outsd[imtstr] = np.sqrt(sdgrid)

        self.logger.debug('\ttime for %s distance=%f' % (imtstr, ddtime))
        self.logger.debug('\ttime for %s correlation=%f' % (imtstr, ctime))
        self.logger.debug('\ttime for %s sigma=%f' % (imtstr, stime))
        self.logger.debug('\ttime for %s rcmatrix=%f' % (imtstr, dtime))
        self.logger.debug('\ttime for %s amp calc=%f' % (imtstr, atime))
        self.logger.debug('\ttime for %s sd calc=%f' % (imtstr, mtime))
        self.logger.debug('total time for %s=%f' %
                          (imtstr, time.time() - time1))

    def _getMaskedGrids(self):
        """
        For each grid in the output, generate a grid with the water areas
        masked out.
        """
        moutgrid = {}
        # We only need to make the mask for one IMT, then we can
        # apply it to all of the other grids.
        refimt = list(self.imt_out_set)[0]
        #
        # Don't know what to do about this warning; hopefully someone
        # will fix maskoceans()
        #
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    category=np.VisibleDeprecationWarning)
            moutgrid[refimt] = \
                maskoceans(self.sx_out.lons, self.sx_out.lats,
                           self.outgrid[refimt], inlands=False, grid=1.25)
        for imtout in self.imt_out_set:
            if imtout == refimt:
                continue
            moutgrid[imtout] = \
                ma.masked_array(self.outgrid[imtout],
                                mask=copy.copy(moutgrid[refimt].mask))
        return moutgrid

    def _getInfo(self, moutgrid):
        """
        Create an info dictionary that can be made into the info.json file.
        """
        #
        # Get the map grade
        #
        mean_rat, mygrade = _get_map_grade(self.do_grid, self.outsd, self.psd,
                                           moutgrid)
        # ------------------------------------------------------------------
        # This is the metadata for creating info.json
        # ------------------------------------------------------------------
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

        # the following items are primarily useful for PDL
        info[ip][ei]['eventsource'] = self.rupture_obj._origin.netid
        info[ip][ei]['netid'] = self.rupture_obj._origin.netid
        info[ip][ei]['eventsourcecode'] = \
            self.rupture_obj._origin.id
        info[ip][ei]['id'] = \
            self.rupture_obj._origin.id
        info[ip][ei]['productcode'] = self.rupture_obj._origin.productcode
        info[ip][ei]['productsource'] = self.config['system']['source_network']
        info[ip][ei]['producttype'] = self.config['system']['product_type']

        info[ip][ei]['fault_ref'] = self.rupture_obj.getReference()
        if 'df2' in self.dataframes:
            info[ip][ei]['intensity_observations'] = \
                str(np.size(self.df2.df['lon']))
        else:
            info[ip][ei]['intensity_observations'] = '0'
        info[ip][ei]['latitude'] = str(self.rx.hypo_lat)
        info[ip][ei]['longitude'] = str(self.rx.hypo_lon)
        info[ip][ei]['location'] = self.rupture_obj._origin.locstring
        info[ip][ei]['magnitude'] = str(self.rx.mag)
        info[ip][ei]['origin_time'] = \
            self.rupture_obj._origin.time.strftime(TIMEFMT)
        if 'df1' in self.dataframes:
            info[ip][ei]['seismic_stations'] = \
                str(np.size(self.df1.df['lon']))
        else:
            info[ip][ei]['seismic_stations'] = '0'
        info[ip][ei]['src_mech'] = self.rupture_obj._origin.mech
        # This AND locaction?
        info[ip][ei]['event_description'] = self.rupture_obj._origin.locstring
        # This AND src_mech?
        info[ip][ei]['event_type'] = self.rupture_obj._origin.mech
        info[op] = {}
        info[op][gm] = {}
        for myimt in self.imt_out_set:
            info[op][gm][myimt] = {}
            if myimt == 'MMI':
                units = 'intensity'
            elif myimt == 'PGV':
                units = 'cms'
            else:
                units = 'ln(g)'
            info[op][gm][myimt]['units'] = units
            if myimt in self.nominal_bias:
                info[op][gm][myimt]['bias'] = str(self.nominal_bias[myimt])
            else:
                info[op][gm][myimt]['bias'] = '-'
            info[op][gm][myimt]['max_grid'] = str(np.max(self.outgrid[myimt]))
            info[op][gm][myimt]['max'] = str(np.max(moutgrid[myimt]))
        info[op][mi] = {}
        info[op][mi]['grid_points'] = {}
        info[op][mi]['grid_points']['longitude'] = str(self.smnx)
        info[op][mi]['grid_points']['latitude'] = str(self.smny)
        info[op][mi]['grid_points']['units'] = ''
        info[op][mi]['grid_spacing'] = {}
        info[op][mi]['grid_spacing']['longitude'] = str(self.smdx)
        info[op][mi]['grid_spacing']['latitude'] = str(self.smdy)
        info[op][mi]['grid_spacing']['units'] = 'degrees'
        info[op][mi]['grid_span'] = {}
        info[op][mi]['grid_span']['longitude'] = str(self.E - self.W)
        info[op][mi]['grid_span']['latitude'] = str(self.N - self.S)
        info[op][mi]['grid_span']['units'] = 'degrees'
        info[op][mi]['min'] = {}
        info[op][mi]['min']['longitude'] = str(self.W)
        info[op][mi]['min']['latitude'] = str(self.S)
        info[op][mi]['min']['units'] = 'degrees'
        info[op][mi]['max'] = {}
        info[op][mi]['max']['longitude'] = str(self.E)
        info[op][mi]['max']['latitude'] = str(self.N)
        info[op][mi]['max']['units'] = 'degrees'
        info[op][un] = {}
        info[op][un]['grade'] = mygrade
        info[op][un]['mean_uncertainty_ratio'] = mean_rat
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
            strftime(TIMEFMT, gmtime())
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
        sjdict = {}
        if self.stations is None:
            return sjdict
        # ------------------------------------------------------------------
        # Compute a bias for all the IMTs in the data frames
        # ------------------------------------------------------------------
        for ndf in self.dataframes:
            sdf = getattr(self, ndf).df
            imts = getattr(self, ndf).imts
            for myimt in imts:
                mybias = self.bias_num[myimt] / \
                    ((1.0 / sdf[myimt + '_pred_tau']**2) +
                     self.bias_den[myimt])
                sdf[myimt + '_bias'] = mybias.flatten()
        # ------------------------------------------------------------------
        # Add the station data. The stationlist object has the original
        # data and produces a GeoJSON object (a dictionary, really), but
        # we need to add peak values and flagging that has been done here.
        # ------------------------------------------------------------------
        #
        # First make a dictionary of distances
        #
        dist_dict = {'df1': {}, 'df2': {}}
        for ndf in self.dataframes:
            dx = getattr(self, ndf).dx
            for dm in get_distance_measures():
                dm_arr = getattr(dx, dm, None)
                if dm_arr is None:
                    continue
                else:
                    dist_dict[ndf][dm] = dm_arr
        #
        # Get the index for each station ID
        #
        sjdict = self.stations.getGeoJson()
        sta_ix = {'df1': {}, 'df2': {}}
        for ndf in self.dataframes:
            sdf = getattr(self, ndf).df
            sta_ix[ndf] = dict(zip(sdf['id'], range(len(sdf['id']))))
        #
        # Now go through the GeoJSON and add various properties and
        # amps from the df_dict dictionaries
        #
        for station in sjdict['features']:
            if station['id'] in sta_ix['df1']:
                ndf = 'df1'
                station['properties']['station_type'] = 'seismic'
            elif station['id'] in sta_ix['df2']:
                ndf = 'df2'
                station['properties']['station_type'] = 'macroseismic'
            else:
                raise ValueError('Unknown station %s in stationlist' %
                                 (station['id']))
            sdf = getattr(self, ndf).df
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
                    and not np.isnan(sdf['PGA'][six]):
                station['properties']['pga'] = \
                    _round_float(np.exp(sdf['PGA'][six]) * 100, 4)
            else:
                station['properties']['pga'] = 'null'

            if 'PGV' in sdf and not sdf['PGV_outliers'][six] \
                    and not np.isnan(sdf['PGV'][six]):
                station['properties']['pgv'] = \
                    _round_float(np.exp(sdf['PGV'][six]), 4)
            else:
                station['properties']['pgv'] = 'null'
            #
            # Add the predictions so we can plot residuals
            #
            station['properties']['predictions'] = []
            for key in sdf.keys():
                if not key.endswith('_pred'):
                    continue
                myamp = sdf[key][six]
                tau_str = 'ln_tau'
                phi_str = 'ln_phi'
                sigma_str = 'ln_sigma'
                bias_str = 'ln_bias'
                if key.startswith('PGV'):
                    value = np.exp(myamp)
                    units = 'cm/s'
                elif key.startswith('MMI'):
                    value = myamp
                    units = 'intensity'
                    tau_str = 'tau'
                    phi_str = 'phi'
                    sigma_str = 'sigma'
                    bias_str = 'bias'
                else:
                    value = np.exp(myamp) * 100
                    units = '%g'
                if self.total_sd_only:
                    mytau = 0
                else:
                    mytau = sdf[key + '_tau'][six]
                myphi = sdf[key + '_phi'][six]
                mysigma = np.sqrt(mytau**2 + myphi**2)
                imt_name = key.lower().replace('_pred', '')
                mybias = sdf[imt_name.upper() + '_bias'][six]
                station['properties']['predictions'].append({
                    'name': imt_name,
                    'value': _round_float(value, 4),
                    'units': units,
                    tau_str: _round_float(mytau, 4),
                    phi_str: _round_float(myphi, 4),
                    sigma_str: _round_float(mysigma, 4),
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
            for channel in station['properties']['channels']:
                for amp in channel['amplitudes']:
                    Name = amp['name'].upper()
                    if sdf[Name + '_outliers'][six]:
                        if amp['flag'] == '0':
                            amp['flag'] = 'T'
                        else:
                            amp['flag'] += 'T'
                    if not amp['flag']:
                        amp['flag'] = '0'
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
        gdict = GeoDict(metadata)
        #
        # Put the Vs30 grid in the output container
        #
        _, units, digits = _get_layer_info('vs30')
        vs30 = Grid2D(self.sx_out.vs30, gdict)
        vs30_metadata = {}
        vs30_metadata['units'] = units
        vs30_metadata['digits'] = digits
        oc.setGrid('vs30', vs30, metadata=vs30_metadata)
        #
        # Now do the distance grids
        #
        dist_metadata = {}
        dist_metadata['units'] = 'km'
        dist_metadata['digits'] = 4
        for dm in get_distance_measures():
            dm_arr = getattr(self.dx_out, dm, None)
            if dm_arr is None:
                continue
            dm_arr_2d = Grid2D(dm_arr.copy(), gdict)
            oc.setGrid('distance_' + dm, dm_arr_2d, metadata=dist_metadata)
        #
        # Output the data and uncertainty grids
        #
        component = self.config['interp']['component']
        for key, value in self.outgrid.items():
            # set the data grid
            mean_grid = Grid2D(value, gdict)
            mean_layername, units, digits = _get_layer_info(key)
            mean_metadata = {'units': units,
                             'digits': digits}

            # set the uncertainty grid
            std_layername, units, digits = _get_layer_info(key + '_sd')
            std_metadata = {'units': units,
                            'digits': digits}
            std_grid = Grid2D(self.outsd[key], gdict.copy())
            oc.setIMTGrids(key,
                           mean_grid, mean_metadata,
                           std_grid, std_metadata,
                           component)

    def _storePointData(self, oc):
        """
        Store point data in the output container.
        """
        #
        # Store the Vs30
        #
        vs30_metadata = {'units': 'm/s',
                         'digits': 4}
        oc.setArray('vs30', self.sx_out.vs30.flatten(),
                    metadata=vs30_metadata)
        #
        # Store the distances
        #
        distance_metadata = {'units': 'km',
                             'digits': 4}
        for dm in get_distance_measures():
            dm_arr = getattr(self.dx_out, dm, None)
            if dm_arr is None:
                continue
            oc.setArray('distance_' + dm, dm_arr.copy().flatten(),
                        metadata=distance_metadata)
        #
        # Store the IMTs
        #
        ascii_ids = np.array(
            [x.encode('ascii') for x in self.idents]).flatten()
        component = self.config['interp']['component']
        for key, value in self.outgrid.items():
            # set the data grid
            mean_layername, units, digits = _get_layer_info(key)
            mean_metadata = {'units': units,
                             'digits': digits}
            # set the uncertainty grid
            std_layername, units, digits = _get_layer_info(key + '_sd')
            std_metadata = {'units': units,
                            'digits': digits}
            oc.setIMTArrays(key, self.sx_out.lons.flatten(),
                            self.sx_out.lats.flatten(),
                            ascii_ids, value.flatten(), mean_metadata,
                            self.outsd[key].flatten(), std_metadata, component)

    def _storeRegressionData(self, oc):
        """
        Average the regression grids in distance bins and output the arrays
        """
        rrup = self.dx_out.rrup
        dx = self.smdx * 111.0
        dmax = np.max(rrup)
        dists = np.arange(0, dmax + dx, dx)
        ndists = np.size(dists)
        rockmean = {}
        soilmean = {}
        rocksd = {}
        soilsd = {}
        mean_dists = np.zeros(ndists - 1)
        for imtstr in self.rockgrid.keys():
            rockmean[imtstr] = np.zeros(ndists - 1)
            soilmean[imtstr] = np.zeros(ndists - 1)
            rocksd[imtstr] = np.zeros(ndists - 1)
            soilsd[imtstr] = np.zeros(ndists - 1)

        bad_ix = []
        for ix in range(ndists - 1):
            ixx = (rrup >= dists[ix]) & (rrup < dists[ix+1])
            mean_dists[ix] = (dists[ix] + dists[ix+1]) / 2.0
            if not np.any(ixx):
                bad_ix.append(ix)
                continue
            for imtstr in self.rockgrid.keys():
                rockmean[imtstr][ix] = np.mean(self.rockgrid[imtstr][ixx])
                soilmean[imtstr][ix] = np.mean(self.soilgrid[imtstr][ixx])
                rocksd[imtstr][ix] = np.mean(self.rocksd[imtstr][ixx])
                soilsd[imtstr][ix] = np.mean(self.soilsd[imtstr][ixx])

        mean_dists = np.delete(mean_dists, bad_ix)
        oc.setArray('regression_distances', mean_dists)
        for imtstr in self.rockgrid.keys():
            rockmean[imtstr] = np.delete(rockmean[imtstr], bad_ix)
            soilmean[imtstr] = np.delete(soilmean[imtstr], bad_ix)
            rocksd[imtstr] = np.delete(rocksd[imtstr], bad_ix)
            soilsd[imtstr] = np.delete(soilsd[imtstr], bad_ix)
            oc.setArray('regression_' + imtstr + '_rock_mean',
                        rockmean[imtstr])
            oc.setArray('regression_' + imtstr + '_soil_mean',
                        soilmean[imtstr])
            oc.setArray('regression_' + imtstr + '_rock_sd',
                        rocksd[imtstr])
            oc.setArray('regression_' + imtstr + '_soil_sd',
                        soilsd[imtstr])

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
        else:
            pe = gmpe

            # -------------------------------------------------------------------------
            # Describe the MultiGMPE
            # -------------------------------------------------------------------------
            if not hasattr(self, '_info'):
                self._info = {
                    'multigmpe': {}
                }
            self._info['multigmpe'][str(oqimt)] = gmpe.describe()

        mean, stddevs = pe.get_mean_and_stddevs(
            copy.deepcopy(sx), self.rx,
            copy.deepcopy(dx), oqimt,
            self.stddev_types)

        if apply_gafs:
            gafs = get_generic_amp_factors(sx, str(oqimt))
            if gafs is not None:
                mean += gafs

        return mean, stddevs


def _round_float(val, digits):
    if val == 'null':
        return val
    return float(('%.' + str(digits) + 'f') % val)


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
        for this_imt in df.imts:
            if not np.isnan(df.df[this_imt][ix]) and \
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
    6, or the map is not a grid, the grade and mean ratio are set to '-'.

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
    mean_rat = '-'
    mygrade = '-'
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
        if mygrade == '-':
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
