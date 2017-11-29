"""
Process a ShakeMap, based on the configuration and data found in
shake_data.hdf, and produce output in shake_result.hdf.
"""
import warnings
import os.path
import time as time
import copy
from time import gmtime, strftime
import json

import numpy as np
import numpy.ma as ma
import numexpr as ne
from mpl_toolkits.basemap import maskoceans
from openquake.hazardlib import imt
import openquake.hazardlib.const as oqconst

# local imports
from .base import CoreModule
from shakelib.rupture.point_rupture import PointRupture
from shakelib.sites import Sites
from shakelib.distance import Distance, get_distance
from shakelib.multigmpe import MultiGMPE
from shakelib.virtualipe import VirtualIPE
from shakelib.utils.utils import get_extent
from shakelib.utils.imt_string import oq_to_file
from shakelib.utils.containers import InputContainer, OutputContainer
from shakelib.utils.distance import geodetic_distance_fast

from mapio.geodict import GeoDict
from mapio.grid2d import Grid2D

from shakemap.utils.config import get_config_paths
from shakemap.utils.utils import get_object_from_config
from shakemap._version import get_versions

#
# TODO: Some constants; these should maybe be in a configuration
# or in a constants module or be GMICE-specific.
#
SM_CONSTS = {'default_mmi_std': 0.3,
             'min_mmi_convert': 4.0}


class ModelModule(CoreModule):
    """
    **model** -- Interpolate ground motions to a grid or list of locations.
    """
    command_name = 'model'

    def execute(self):
        """Interpolate ground motions to a grid or list of locations.

        Raises:
            NotADirectoryError: When the event data directory does not exist.
            FileNotFoundError: When the the shake_data HDF file does not exist.
        """
        #
        # Find the shake_data file
        #
        self.logger.info('Inside model')
        install_path, data_path = get_config_paths()
        datadir = os.path.join(data_path, self._eventid, 'current')
        if not os.path.isdir(datadir):
            raise NotADirectoryError('%s is not a valid directory.' % datadir)

        datafile = os.path.join(datadir, 'shake_data.hdf')
        if not os.path.isfile(datafile):
            raise FileNotFoundError('%s does not exist.' % datafile)

        # ------------------------------------------------------------------
        # Make the input container and extract the config
        # ------------------------------------------------------------------
        ic = InputContainer.load(datafile)
        config = ic.getConfig()
        # ------------------------------------------------------------------
        # Instantiate the gmpe, gmice, ipe, and ccf
        # Here we make a placeholder gmpe so that we can make the
        # rupture and distance contexts; later we'll make the
        # IMT-specific gmpes
        # ------------------------------------------------------------------
        default_gmpe = MultiGMPE.from_config(config)

        gmice = get_object_from_config('gmice', 'modeling', config)

        if config['ipe_modules'][config['modeling']['ipe']][0] == 'VirtualIPE':
            pgv_imt = imt.from_string('PGV')
            ipe_gmpe = MultiGMPE.from_config(config, filter_imt=pgv_imt)
            ipe = VirtualIPE.fromFuncs(ipe_gmpe, gmice)
        else:
            ipe = get_object_from_config('ipe', 'modeling', config)
        # ------------------------------------------------------------------
        # Bias parameters
        # ------------------------------------------------------------------
        do_bias         = config['modeling']['bias']['do_bias']
        bias_max_range  = config['modeling']['bias']['max_range']
        bias_max_mag    = config['modeling']['bias']['max_mag']
        bias_max_dsigma = config['modeling']['bias']['max_delta_sigma']
        # ------------------------------------------------------------------
        # Outlier parameters
        # ------------------------------------------------------------------
        outlier_deviation_level = config['data']['outlier']['max_deviation']
        outlier_max_mag         = config['data']['outlier']['max_mag']
        # ------------------------------------------------------------------
        # These are the IMTs we want to make
        # ------------------------------------------------------------------
        imt_out_set_str = set(config['interp']['imt_list'])
        # ------------------------------------------------------------------
        # Get the rupture object and rupture context
        # ------------------------------------------------------------------
        rupture_obj = ic.getRupture()
        if 'mechanism' in config['modeling']:
            rupture_obj._origin.setMechanism(
                    mech=config['modeling']['mechanism'])
        rx = rupture_obj.getRuptureContext([default_gmpe])
        if rx.rake is None:
            rx.rake = 0
        # ------------------------------------------------------------------
        # Get the Vs30 file name
        # ------------------------------------------------------------------
        vs30default = config['data']['vs30default']
        vs30_file = config['data']['vs30file']
        if not vs30_file:
            vs30_file = None
        # ------------------------------------------------------------------
        # The output locations: either a grid or a list of points
        # ------------------------------------------------------------------
        smdx = config['interp']['prediction_location']['xres']
        smdy = config['interp']['prediction_location']['yres']
        # TODO: functionize these blocks?
        if config['interp']['prediction_location']['file'] and \
           config['interp']['prediction_location']['file'] != 'None':
            #
            # FILE: Open the file and get the output points
            #
            do_grid = False
            in_sites = np.genfromtxt(
                            config['interp']['prediction_location']['file'],
                            autostrip=True, unpack=True,
                            dtype=[np.float, np.float, np.float, '<U80'])
            lons, lats, vs30, idents = zip(*in_sites)
            lons = np.array(lons).reshape(1, -1)
            lats = np.array(lats).reshape(1, -1)
            vs30 = np.array(vs30).reshape(1, -1)
            depths = np.zeros_like(lats)
            W = np.min(lons)
            E = np.max(lons)
            S = np.min(lats)
            N = np.max(lats)
            smnx = np.size(lons)
            smny = 1
            dist_obj_out = Distance(default_gmpe, lons, lats, depths,
                                    rupture_obj)

            sites_obj_out = Sites.fromBounds(W, E, S, N, smdx, smdy,
                                             defaultVs30=vs30default,
                                             vs30File=vs30_file)

            sx_out_soil = sites_obj_out.getSitesContext({'lats': lats,
                                                         'lons': lons})
        else:
            #
            # GRID: Figure out the grid parameters and get output points
            #
            do_grid = True
            smdx = config['interp']['prediction_location']['xres']
            smdy = config['interp']['prediction_location']['yres']

            if config['interp']['prediction_location']['extent']:
                W, S, E, N = config['interp']['prediction_location']['extent']
            else:
                W, E, S, N = get_extent(rupture_obj)

            sites_obj_out = Sites.fromBounds(W, E, S, N, smdx, smdy,
                                             defaultVs30=vs30default,
                                             vs30File=vs30_file)
            smnx, smny = sites_obj_out.getNxNy()
            sx_out_soil = sites_obj_out.getSitesContext()
            lons, lats = np.meshgrid(sx_out_soil.lons, sx_out_soil.lats)
            lons = np.flipud(lons).flatten()
            lats = np.flipud(lats).flatten()
            depths = np.zeros_like(lats)
            dist_obj_out = Distance.fromSites(default_gmpe, sites_obj_out,
                                              rupture_obj)

        dx_out = dist_obj_out.getDistanceContext()

        lons_out_rad = np.radians(lons)
        lats_out_rad = np.radians(lats)
        # ------------------------------------------------------------------
        # Station data
        # ------------------------------------------------------------------
        stations = ic.getStationList()
        gmpe_sd_types = default_gmpe.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if len(gmpe_sd_types) == 1:
            total_sd_only = True
            stddev_types = [oqconst.StdDev.TOTAL]
        else:
            total_sd_only = False
            stddev_types = [oqconst.StdDev.TOTAL, oqconst.StdDev.INTER_EVENT,
                            oqconst.StdDev.INTRA_EVENT]
        #
        # df1 holds the instrumented data (PGA, PGV, SA)
        # df2 holds the non-instrumented data (MMI)
        #
        df_dict = {'df1': None, 'df2': None}
        imt_in_str_dict = {'df1': None, 'df2': None}
        imt_in_dict = {'df1': None, 'df2': None}
        imt_in_str_set = set()
        sx_dict = {'df1': None, 'df2': None}
        dx_dict = {'df1': None, 'df2': None}
        if stations is not None:
            df_dict['df1'] = stations.getStationDictionary(instrumented=True)
            df_dict['df2'] = stations.getStationDictionary(instrumented=False)
            #
            # Make the sites and distance contexts for each dictionary then
            # compute the predictions for the IMTs in that dictionary.
            #
            for ndf, df in df_dict.items():
                if not df:
                    imt_in_str_dict[ndf] = set()
                    imt_in_dict[ndf] = set()
                    continue
                #
                # Make lists of the input IMTs
                #
                imt_in_str_dict[ndf] = set(
                        [x for x in df.keys() if x in ('PGA', 'PGV', 'MMI') or
                         x.startswith('SA(')])
                imt_in_dict[ndf] = set(
                        [imt.from_string(x) for x in imt_in_str_dict[ndf]])
                imt_in_str_set |= imt_in_str_dict[ndf]
                #
                # Get the sites and distance contexts
                #
                df['depth'] = np.zeros_like(df['lon'])
                lldict = {'lons': df['lon'], 'lats': df['lat']}
                sx_dict[ndf] = sites_obj_out.getSitesContext(lldict)
                dist_obj = Distance(default_gmpe, df['lon'], df['lat'],
                                    df['depth'], rupture_obj)
                dx_dict[ndf] = dist_obj.getDistanceContext()
                #
                # Do the predictions and other bookkeeping for each IMT
                #
                for imtstr in imt_in_str_dict[ndf]:
                    oqimt = imt.from_string(imtstr)
                    gmpe = None
                    if imtstr != 'MMI':
                        gmpe = MultiGMPE.from_config(config, filter_imt=oqimt)
                    pmean, pstddev = gmas(ipe, gmpe, sx_dict[ndf], rx,
                                          dx_dict[ndf], oqimt, stddev_types)
                    df[imtstr + '_pred'] = pmean
                    df[imtstr + '_pred_sigma'] = pstddev[0]
                    if total_sd_only is True:
                        tau_guess = get_default_std_inter()
                        df[imtstr + '_pred_tau'] = np.full_like(
                                df[imtstr + '_pred'], tau_guess)
                        df[imtstr + '_pred_phi'] = np.sqrt(pstddev[0]**2 -
                                                           tau_guess**2)
                    else:
                        df[imtstr + '_pred_tau'] = pstddev[1]
                        df[imtstr + '_pred_phi'] = pstddev[2]
                    #
                    # Compute the total residual
                    #
                    df[imtstr + '_residual'] = \
                        df[imtstr] - df[imtstr + '_pred']
                    # ----------------------------------------------------------
                    # Do the outlier flagging if we don't have a fault and
                    # the event magnitude is over the limit
                    # ----------------------------------------------------------
                    if not isinstance(rupture_obj, PointRupture) or \
                       rx.mag <= outlier_max_mag:
                        #
                        # turn off nan warnings for this statement
                        #
                        np.seterr(invalid='ignore')
                        flagged = \
                            np.abs(df[imtstr + '_residual']) > \
                            outlier_deviation_level * df[imtstr + '_pred_sigma']
                        np.seterr(invalid='warn')

                        self.logger.info('IMT: %s, flagged: %d' %
                                    (imtstr, np.sum(flagged)))
                        df[imtstr + '_outliers'] = flagged
                    else:
                        df[imtstr + '_outliers'] = np.full(df[imtstr].shape,
                                                           False,
                                                           dtype=np.bool)
                    #
                    # Make the uncertainty arrays for any IMTs that don't
                    # have them.
                    #
                    if (imtstr + '_sd') not in df:
                        if imtstr == 'MMI':
                            sdval = SM_CONSTS['default_mmi_std']
                        else:
                            sdval = 0.0
                        df[imtstr + '_sd'] = np.full_like(df['lon'], sdval)
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
                                  df['depth'], rupture_obj)
                df['rrup'] = dd['rrup']
        df1 = df_dict['df1']
        df2 = df_dict['df2']

    # %%
        # ------------------------------------------------------------------
        # Compute all the IMTs possible from MMI
        # This logic needs to be revisited. We should probably make what
        # we have to to do the CMS to make the needed output IMTs, but
        # for now, we're just going to use what we have and the ccf.
        # ------------------------------------------------------------------
        if df2:
            for gmice_imt in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
                if imt.SA == gmice_imt:
                    iterlist = gmice.DEFINED_FOR_SA_PERIODS
                else:
                    iterlist = [None]
                for period in iterlist:
                    oqimt = gmice_imt(period)
                    imtstr = str(oqimt)

                    np.seterr(invalid='ignore')
                    df2[imtstr], _ = gmice.getGMfromMI(df2['MMI'], oqimt,
                                                       dists=df2['rrup'],
                                                       mag=rx.mag)
                    df2[imtstr][df2['MMI'] < SM_CONSTS['min_mmi_convert']] = \
                        np.nan
                    np.seterr(invalid='warn')
                    df2[imtstr + '_sd'] = \
                        np.full_like(df2['MMI'], gmice.getMI2GMsd()[oqimt])
                    imt_in_str_dict['df2'].add(imtstr)
                    #
                    # Get the predictions and stddevs
                    #
                    gmpe = MultiGMPE.from_config(config, filter_imt=oqimt)
                    pmean, pstddev = gmas(ipe, gmpe, sx_dict['df2'], rx,
                                          dx_dict['df2'], oqimt, stddev_types)
                    df2[imtstr + '_pred'] = pmean
                    df2[imtstr + '_pred_sigma'] = pstddev[0]
                    if total_sd_only is True:
                        tau_guess = get_default_std_inter()
                        df2[imtstr + '_pred_tau'] = np.full_like(
                                df2[imtstr + '_pred'], tau_guess)
                        df2[imtstr + '_pred_phi'] = np.sqrt(pstddev[0]**2 -
                                                            tau_guess**2)
                    else:
                        df2[imtstr + '_pred_tau'] = pstddev[1]
                        df2[imtstr + '_pred_phi'] = pstddev[2]
                    df2[imtstr + '_residual'] = df2[str(oqimt)] - pmean
                    df2[imtstr + '_outliers'] = np.full(pmean.shape, False,
                                                        dtype=np.bool)

        #
        # Now make derived MMI from the best available PGM; This is ugly and
        # it would be nice to have a more deterministic way of doing it
        #
        if df1:
            imtstr = None
            if 'PGV' in df1 \
                    and imt.PGV in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
                imtstr = 'PGV'
            elif 'PGA' in df1 \
                    and imt.PGA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES:
                imtstr = 'PGA'
            elif 'SA(1.0)' in df1 \
                    and imt.SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES \
                    and 1.0 in gmice.DEFINED_FOR_SA_PERIODS:
                imtstr = 'SA(1.0)'
            elif 'SA(0.3)' in df1 \
                    and imt.SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES \
                    and 0.3 in gmice.DEFINED_FOR_SA_PERIODS:
                imtstr = 'SA(0.3)'
            elif 'SA(3.0)' in df1 \
                    and imt.SA in gmice.DEFINED_FOR_INTENSITY_MEASURE_TYPES \
                    and 3.0 in gmice.DEFINED_FOR_SA_PERIODS:
                imtstr = 'SA(3.0)'

            if imtstr is not None:
                oqimt = imt.from_string(imtstr)
                np.seterr(invalid='ignore')
                df1['MMI'], _ = gmice.getMIfromGM(df1[imtstr], oqimt,
                                                  dists=df1['rrup'],
                                                  mag=rx.mag)
                np.seterr(invalid='warn')
                df1['MMI_sd'] = np.full_like(df1[imtstr],
                                             gmice.getGM2MIsd()[oqimt])
                imt_in_str_dict['df1'].add('MMI')
                #
                # Get the prediction and stddevs
                #
                gmpe = None
                pmean, pstddev = gmas(ipe, gmpe, sx_dict['df1'], rx,
                                      dx_dict['df1'], imt.from_string('MMI'),
                                      stddev_types)
                df1['MMI' + '_pred'] = pmean
                df1['MMI' + '_pred_sigma'] = pstddev[0]
                if total_sd_only is True:
                    tau_guess = get_default_std_inter()
                    df1['MMI' + '_pred_tau'] = np.full_like(
                            df1['MMI' + '_pred'], tau_guess)
                    df1['MMI' + '_pred_phi'] = np.sqrt(pstddev[0]**2 -
                                                       tau_guess**2)
                else:
                    df1['MMI' + '_pred_tau'] = pstddev[1]
                    df1['MMI' + '_pred_phi'] = pstddev[2]
                df1['MMI' + '_residual'] = df1['MMI'] - pmean
                df1['MMI' + '_outliers'] = np.full(pmean.shape, False,
                                                   dtype=np.bool)

        imt_per = get_period_array(imt_in_str_dict['df1'],
                                   imt_in_str_dict['df2'],
                                   imt_out_set_str)
        ccf = get_object_from_config('ccf', 'modeling', config, imt_per)

        imt_per_ix = {str(per): ix for ix, per in enumerate(imt_per)}
    # %%
        # ------------------------------------------------------------------
        # Do the MVN
        # ------------------------------------------------------------------
        #
        # First do the bias for all of the input and output IMTs. Hold on
        # to some of the products that will be used for the interpolation.
        #
        nominal_bias = {}
        bias_num = {}   # The numerator term of the bias
        bias_den = {}   # The denominator term of the bias
        outgrid = {}
        outsd = {}
        psd = {}
        tsd = {}

        sta_period_ix = {}
        sta_lons_rad = {}
        sta_lats_rad = {}
        sta_resids = {}
        sta_phi = {}
        sta_sig_total = {}
        sta_tau = {}
        sta_imtstr = {}
        corr_adj = {}
        corr22 = {}

        #
        # Compute a bias for all of the IMTs in the inputs and outputs
        #
        combined_imt_str = imt_out_set_str
        if imt_in_str_dict['df1']:
            combined_imt_str |= imt_in_str_dict['df1']
        if imt_in_str_dict['df2']:
            combined_imt_str |= imt_in_str_dict['df2']
        for imtstr in combined_imt_str:
            time1 = time.time()
            #
            # Get the index of the (pseudo-) period of the output IMT
            #
            outperiod_ix = get_period_index_from_imt_str(imtstr, imt_per_ix)
            #
            # Fill the station arrays
            #
            period_ix = []
            lons_rad = []
            lats_rad = []
            resids = []
            phi = []
            sig_total = []
            tau = []
            imtlist = []
            for ndf, sdf in df_dict.items():
                if not sdf:
                    continue
                for i in range(np.size(sdf['lon'])):
                    for imtin in get_sta_imts(imtstr, sdf, i,
                                              imt_in_str_dict[ndf]):
                        imtlist.append(imtin)
                        inperiod_ix = get_period_index_from_imt_str(
                                imtin, imt_per_ix)
                        period_ix.append(inperiod_ix)
                        lons_rad.append(sdf['lon_rad'][i])
                        lats_rad.append(sdf['lat_rad'][i])
                        resids.append(sdf[imtin + '_residual'][i])
                        tau.append(sdf[imtin + '_pred_tau'][i])
                        _phi = sdf[imtin + '_pred_phi'][i]
                        phi.append(_phi)
                        sig_total.append(
                                np.sqrt(_phi**2 + sdf[imtin + '_sd'][i]**2))
            sta_imtstr[imtstr] = np.array(imtlist)
            sta_period_ix[imtstr] = np.array(period_ix).reshape((-1, 1))
            sta_lons_rad[imtstr] = np.array(lons_rad).reshape((-1, 1))
            sta_lats_rad[imtstr] = np.array(lats_rad).reshape((-1, 1))
            sta_resids[imtstr] = np.array(resids).reshape((-1, 1))
            sta_tau[imtstr] = np.array(tau).reshape((-1, 1))
            sta_phi[imtstr] = np.array(phi).reshape((-1, 1))
            sta_sig_total[imtstr] = np.array(sig_total).reshape((-1, 1))
            if len(lons_rad) == 0:
                bias_num[imtstr] = 0.0
                bias_den[imtstr] = 0.0
                continue
            #
            # This builds the omega factors to apply to the covariance
            #
            corr_adj[imtstr] = (sta_phi[imtstr] / sta_sig_total[imtstr])
            corr_adj22 = corr_adj[imtstr] * corr_adj[imtstr].T
            np.fill_diagonal(corr_adj22, 1.0)
            #
            # Build the covariance matrix of the residuals
            #
            dist22 = geodetic_distance_fast(sta_lons_rad[imtstr],
                                            sta_lats_rad[imtstr],
                                            sta_lons_rad[imtstr].T,
                                            sta_lats_rad[imtstr].T)
            d22_rows, d22_cols = np.shape(dist22)  # should be square
            t1_22 = np.tile(sta_period_ix[imtstr], (1, d22_cols))
            t2_22 = np.tile(sta_period_ix[imtstr].T, (d22_rows, 1))
            corr22[imtstr] = \
                ccf.getCorrelation(t1_22, t2_22, dist22) * corr_adj22
            sigma22 = corr22[imtstr] * (sta_phi[imtstr] * sta_phi[imtstr].T)
            sigma22inv = np.linalg.pinv(sigma22)
            #
            # Compute the bias numerator and denominator pieces
            #
            if do_bias and (not isinstance(rupture_obj, PointRupture)
                            or rx.mag <= bias_max_mag):
                #
                # Get the correlation between the inputs and outputs
                #
                out_ix_arr = np.full_like(sta_period_ix[imtstr], outperiod_ix)
                dist_arr = np.zeros_like(out_ix_arr, dtype=float)
                P = ccf.getCorrelation(sta_period_ix[imtstr], out_ix_arr,
                                       dist_arr)
                bias_den[imtstr] = P.T.dot(sigma22inv.dot(P))
                bias_num[imtstr] = \
                    P.T.dot(sigma22inv.dot(P * sta_resids[imtstr]))
            else:
                bias_num[imtstr] = 0.0
                bias_den[imtstr] = 0.0

            nom_tau = np.mean(sta_tau[imtstr].flatten())
            nom_variance = 1.0 / ((1.0 / nom_tau**2) + bias_den[imtstr])
            nominal_bias[imtstr] = bias_num[imtstr] * nom_variance
            bias_time = time.time() - time1
            #
            # Print the nominal values of the bias and its stddev
            #
            self.logger.info('%s: nom bias %f nom stddev %f; %d stations (time=%f sec)' %
                 (imtstr, nominal_bias[imtstr], np.sqrt(nom_variance),
                 np.size(sta_lons_rad[imtstr]), bias_time))
        #
        # End bias
        #
        #
        # Now do the MVN with the intra-event residuals
        #
        for imtstr in imt_out_set_str:
            time1 = time.time()
            #
            # Get the index of the (pesudo-) period of the output IMT
            #
            outperiod_ix = get_period_index_from_imt_str(imtstr, imt_per_ix)
            #
            # Get the predictions at the output points
            #
            oqimt = imt.from_string(imtstr)
            gmpe = None
            if imtstr != 'MMI':
                gmpe = MultiGMPE.from_config(config, filter_imt=oqimt)
            pout_mean, pout_sd = gmas(ipe, gmpe, sx_out_soil, rx, dx_out,
                                      oqimt, stddev_types)
            #
            # If there are no data, just use the unbiased prediction
            # and the total stddev
            #
            if np.size(sta_lons_rad[imtstr]) == 0:
                outgrid[imtstr] = pout_mean
                outsd[imtstr] = pout_sd[0]
                continue
            #
            # Get an array of the within-event standard deviations for the
            # output IMT at the output points
            #
            if total_sd_only is True:
                psd[imtstr] = np.sqrt(pout_sd[0]**2 -
                                      get_default_std_inter()**2)
                tsd[imtstr] = np.full_like(psd[imtstr],
                                           get_default_std_inter())
            else:
                psd[imtstr] = pout_sd[2]
                tsd[imtstr] = pout_sd[1]
            pout_sd2 = np.power(psd[imtstr], 2.0)
            #
            # Bias the predictions, and add the residual variance to
            # phi
            #
            out_bias_var = 1.0 / ((1.0 / tsd[imtstr]**2) + bias_den[imtstr])
            out_bias = bias_num[imtstr] * out_bias_var
            pout_mean += out_bias
            psd[imtstr] = np.sqrt(psd[imtstr]**2 + out_bias_var)
            #
            # Unbias the residual array
            #
            for i in range(np.size(sta_lons_rad[imtstr])):
                imtin = sta_imtstr[imtstr][i]
                in_bias_var = 1.0 / ((1.0 / sta_tau[imtstr][i, 0]**2) +
                                     bias_den[imtin])
                in_bias = bias_num[imtin] * in_bias_var
                sta_resids[imtstr][i, 0] -= in_bias
                sta_phi[imtstr][i, 0] = np.sqrt(sta_phi[imtstr][i, 0]**2 +
                                                in_bias_var)
            #
            # Now do the MVN itself...
            #
            dtime = mtime = ddtime = ctime = stime = atime = 0
            #
            # Rebuild sigma22_inv now that we have updated phi
            #
            sigma22 = corr22[imtstr] * (sta_phi[imtstr] * sta_phi[imtstr].T)
            sigma22inv = np.linalg.pinv(sigma22)

            ampgrid = np.zeros_like(pout_mean)
            sdgrid = np.zeros_like(pout_mean)
            corr_adj12 = corr_adj[imtstr] * np.ones((1, smnx))
            for iy in range(smny):
                ss = iy * smnx
                se = (iy + 1) * smnx
                time4 = time.time()
                dist12 = geodetic_distance_fast(
                        lons_out_rad[ss:se].reshape(1, -1),
                        lats_out_rad[ss:se].reshape(1, -1),
                        sta_lons_rad[imtstr],
                        sta_lats_rad[imtstr])
                t2_12 = np.full(dist12.shape, outperiod_ix, dtype=np.int)
                d12_rows, d12_cols = np.shape(dist12)
                t1_12 = np.tile(sta_period_ix[imtstr], (1, d12_cols))
                ddtime += time.time() - time4
                time4 = time.time()
                corr12 = ccf.getCorrelation(t1_12, t2_12, dist12)
                ctime += time.time() - time4
                time4 = time.time()
                sdarr = psd[imtstr][iy, :].reshape((1, -1))
                ss = sta_phi[imtstr]
                sigma12 = ne.evaluate("corr12 * corr_adj12 * (ss * sdarr)").T
                stime += time.time() - time4
                time4 = time.time()
                rcmatrix = sigma12.dot(sigma22inv)
                dtime += time.time() - time4
                time4 = time.time()
                ampgrid[iy, :] = \
                    pout_mean[iy, :] + \
                    (corr_adj12.T *
                     rcmatrix).dot(sta_resids[imtstr]).reshape((-1,))
                atime += time.time() - time4
                time4 = time.time()
#        sdgrid[ss:se] = pout_sd2[ss:se] - np.diag(rcmatrix.dot(sigma12))
                sdgrid[iy, :] = \
                    pout_sd2[iy, :] - np.sum(rcmatrix * sigma12, axis=1)
                mtime += time.time() - time4

            outgrid[imtstr] = ampgrid
            sdgrid[sdgrid < 0] = 0
            outsd[imtstr] = np.sqrt(sdgrid)

            self.logger.debug('\ttime for %s distance=%f' % (imtstr, ddtime))
            self.logger.debug('\ttime for %s correlation=%f' % (imtstr, ctime))
            self.logger.debug('\ttime for %s sigma=%f' % (imtstr, stime))
            self.logger.debug('\ttime for %s rcmatrix=%f' % (imtstr, dtime))
            self.logger.debug('\ttime for %s amp calc=%f' % (imtstr, atime))
            self.logger.debug('\ttime for %s sd calc=%f' % (imtstr, mtime))
            self.logger.info('total time for %s=%f' % (imtstr, time.time() - time1))

    # %%
        # ------------------------------------------------------------------
        # Output the data and metadata
        # ------------------------------------------------------------------
        product_path = os.path.join(datadir, 'products')
        if not os.path.isdir(product_path):
            os.mkdir(product_path)
        oc = OutputContainer.create(os.path.join(product_path,
                                                 'shake_result.hdf'))
        # ------------------------------------------------------------------
        # Might as well stick the whole config in the result
        # ------------------------------------------------------------------
        oc.setDictionary('config', config)

        #
        # We're going to need masked arrays of the output grids later, so
        # make them now. We only need to make the mask once, then we can
        # apply it to all of the other grids.
        #
        moutgrid = {}
        refimt = list(imt_out_set_str)[0]
        #
        # Don't know what to do about this warning; hopefully someone
        # will fix maskoceans()
        #
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",
                                    category=np.VisibleDeprecationWarning)
            moutgrid[refimt] = \
                maskoceans(lons.reshape((smny, smnx)),
                           lats.reshape((smny, smnx)),
                           outgrid[refimt], inlands=False, grid=1.25)
        for imtout in imt_out_set_str:
            if imtout == refimt:
                continue
            moutgrid[imtout] = \
                ma.masked_array(outgrid[imtout],
                                mask=copy.copy(moutgrid[refimt].mask))
        #
        # Get the map grade
        #
        mean_rat, mygrade = get_map_grade(do_grid, outsd, psd, moutgrid)
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
        info = {}
        info[ip] = {}
        info[ip][ei] = {}
        info[ip][ei]['depth'] = str(rx.hypo_depth)
        info[ip][ei]['event_id'] = self._eventid
        info[ip][ei]['fault_ref'] = rupture_obj.getReference()
        info[ip][ei]['intensity_observations'] = \
            str(np.size(df2['lon'])) if df2 else '0'
        info[ip][ei]['latitude'] = str(rx.hypo_lat)
        info[ip][ei]['longitude'] = str(rx.hypo_lon)
        info[ip][ei]['location'] = rupture_obj._origin.locstring
        info[ip][ei]['magnitude'] = str(rx.mag)
        info[ip][ei]['origin_time'] = \
            rupture_obj._origin.time.strftime('%Y-%m-%d %H:%M:%S')
        info[ip][ei]['seismic_stations'] = \
            str(np.size(df1['lon'])) if df1 else '0'
        info[ip][ei]['src_mech'] = rupture_obj._origin.mech
        info[ip][ei]['event_description'] = rupture_obj._origin.locstring  ### This AND locaction?
        info[ip][ei]['event_type'] = rupture_obj._origin.mech              ### This AND src_mech?
        info[op] = {}
        info[op][gm] = {}
        for myimt in imt_out_set_str:
            info[op][gm][myimt] = {}
            if myimt == 'MMI':
                units = 'intensity'
            elif myimt == 'PGV':
                units = 'cms'
            else:
                units = 'ln(g)'
            info[op][gm][myimt]['units'] = units
            info[op][gm][myimt]['bias'] = \
                str(nominal_bias[myimt]) if myimt in nominal_bias else '-'
            info[op][gm][myimt]['max_grid'] = str(np.max(outgrid[myimt]))
            info[op][gm][myimt]['max'] = str(np.max(moutgrid[myimt]))
        info[op][mi] = {}
        info[op][mi]['grid_points'] = {}
        info[op][mi]['grid_points']['longitude'] = str(smnx)
        info[op][mi]['grid_points']['latitude'] = str(smny)
        info[op][mi]['grid_points']['units'] = ''
        info[op][mi]['grid_spacing'] = {}
        info[op][mi]['grid_spacing']['longitude'] = str(smdx)
        info[op][mi]['grid_spacing']['latitude'] = str(smdy)
        info[op][mi]['grid_spacing']['units'] = 'degrees'
        info[op][mi]['grid_span'] = {}
        info[op][mi]['grid_span']['longitude'] = str(E - W)
        info[op][mi]['grid_span']['latitude'] = str(N - S)
        info[op][mi]['grid_span']['units'] = 'degrees'
        info[op][mi]['min'] = {}
        info[op][mi]['min']['longitude'] = str(W)
        info[op][mi]['min']['latitude'] = str(S)
        info[op][mi]['min']['units'] = 'degrees'
        info[op][mi]['max'] = {}
        info[op][mi]['max']['longitude'] = str(E)
        info[op][mi]['max']['latitude'] = str(N)
        info[op][mi]['max']['units'] = 'degrees'
        info[op][un] = {}
        info[op][un]['grade'] = mygrade
        info[op][un]['mean_uncertainty_ratio'] = mean_rat
        info[op][un]['total_flagged_mi'] = \
            str(np.sum(df2['MMI_outliers'] |
                       np.isnan(df2['MMI']))) if df2 else '0'
        if df1:
            all_flagged = np.full(df1['lon'].shape, False, dtype=np.bool)
            for imtstr in imt_in_str_dict['df1']:
                if 'MMI' in imtstr:
                    continue
                all_flagged |= \
                    df1[imtstr + '_outliers'] | np.isnan(df1[imtstr])
            info[op][un]['total_flagged_pgm'] = str(np.sum(all_flagged))
        else:
            info[op][un]['total_flagged_pgm'] = '0'
        info[pp] = {}
        info[pp][gmm] = {}
        info[pp][gmm]['gmpe'] = {}
        info[pp][gmm]['gmpe']['module'] = str(config['modeling']['gmpe'])
        info[pp][gmm]['gmpe']['reference'] = ''
        info[pp][gmm]['ipe'] = {}
        info[pp][gmm]['ipe']['module'] = \
            str(config['ipe_modules'][config['modeling']['ipe']][0])
        info[pp][gmm]['ipe']['reference'] = ''
        info[pp][gmm]['gmice'] = {}
        info[pp][gmm]['gmice']['module'] = \
            str(config['gmice_modules'][config['modeling']['gmice']][0])
        info[pp][gmm]['gmice']['reference'] = ''
        info[pp][gmm]['ccf'] = {}
        info[pp][gmm]['ccf']['module'] = \
            str(config['ccf_modules'][config['modeling']['ccf']][0])
        info[pp][gmm]['ccf']['reference'] = ''
        info[pp][gmm]['basin_correction'] = {}
        info[pp][gmm]['basin_correction']['module'] = 'None'
        info[pp][gmm]['basin_correction']['reference'] = ''
        info[pp][gmm]['directivity'] = {}
        info[pp][gmm]['directivity']['module'] = 'None'
        info[pp][gmm]['directivity']['reference'] = ''
        info[pp][ms] = {}
        info[pp][ms]['bias_max_dsigma'] = str(bias_max_dsigma)
        info[pp][ms]['bias_max_mag'] = str(bias_max_mag)
        info[pp][ms]['bias_max_range'] = str(bias_max_range)
        info[pp][ms]['median_dist'] = 'True'
        info[pp][ms]['outlier_deviation_level'] = str(outlier_deviation_level)
        info[pp][ms]['outlier_max_mag'] = str(outlier_max_mag)
        info[pp][sv] = {}
        info[pp][sv]['shakemap_revision'] = get_versions()['version']
        info[pp][sv]['shakemap_revision_id'] = \
            get_versions()['full-revisionid']
        info[pp][sv]['process_time'] = \
            strftime("%Y-%m-%d %H:%M:%S%Z", gmtime())
        info[pp][sv]['map_version'] = ic.getVersionHistory()['history'][-1][2]
        info[pp][sv]['map_data_history'] = ic.getVersionHistory()['history']
        info[pp][sv]['map_status'] = config['system']['map_status']
        info[pp][sr] = {}
        info[pp][sr]['vs30default'] = str(vs30default)
        info[pp][sr]['site_correction'] = 'GMPE native'

        oc.setString('info.json', json.dumps(info))

        # ------------------------------------------------------------------
        # Add the rupture JSON as a text string
        # ------------------------------------------------------------------
        oc.setString('rupture.json', json.dumps(rupture_obj._geojson))

        # ------------------------------------------------------------------
        # Add the station data. The stationlist object has the original
        # data and produces a GeoJSON object (a dictionary, really), but
        # we need to add peak values and flagging that has been done here.
        # ------------------------------------------------------------------
        if stations:
            sjdict = stations.getGeoJson()
            sta_ix = {}
            for ndf, sdf in df_dict.items():
                if not sdf:
                    sta_ix[ndf] = {}
                else:
                    sta_ix[ndf] = dict(zip(sdf['id'], range(len(sdf['id']))))
            for station in sjdict['features']:
                if station['id'] in sta_ix['df1']:
                    sdf = df1
                    six = sta_ix['df1'][station['id']]
                elif station['id'] in sta_ix['df2']:
                    sdf = df2
                    six = sta_ix['df2'][station['id']]
                else:
                    raise ValueError('Unknown station %s in stationlist' %
                                     (station['id']))
                if 'MMI' in sdf and not sdf['MMI_outliers'][six]:
                    station['properties']['intensity'] = \
                        '%.1f' % sdf['MMI'][six]
                else:
                    station['properties']['intensity'] = 'null'

                if 'PGA' in sdf and not sdf['PGA_outliers'][six]:
                    station['properties']['pga'] = \
                        '%.4f' % (np.exp(sdf['PGA'][six]) * 100)
                else:
                    station['properties']['pga'] = 'null'
                if 'PGV' in sdf and not sdf['PGV_outliers'][six]:
                    station['properties']['pgv'] = '%.4f' % \
                                                (np.exp(sdf['PGV'][six]))
                else:
                    station['properties']['pgv'] = 'null'
                #
                # Add the predictions so we can plot residuals
                #
                station['properties']['predictions'] = {}
                for key in sdf.keys():
                    if not key.endswith('_pred'):
                        continue
                    myamp = sdf[key][six]
                    if key.startswith('PGV'):
                        value = np.exp(myamp)
                        units = 'cm/s'
                    elif key.startswith('MMI'):
                        value = myamp
                        units = 'intensity'
                    else:
                        value = np.exp(myamp) * 100
                        units = '%g'
                    station['properties']['predictions'][key.lower()] = {
                            'value': value,
                            'units': units
                    }
                station['properties']['distance'] = '%.2f' % sdf['rrup'][six]
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
        else:
            sjdict = {}

        oc.setString('stationlist.json', json.dumps(sjdict))

        # ------------------------------------------------------------------
        # Add the output grids or points to the output; include some
        # metadata.
        # ------------------------------------------------------------------
        metadata = {}
        if do_grid:
            metadata['xmin'] = W
            metadata['xmax'] = E
            metadata['ymin'] = S
            metadata['ymax'] = N
            metadata['nx'] = smnx
            metadata['ny'] = smny
            metadata['dx'] = smdx
            metadata['dy'] = smdy
            gdict = GeoDict(metadata)
            layername, units, digits = get_layer_info('vs30')
            vs30 = Grid2D(sx_out_soil.vs30, gdict)
            vs30_metadata = {}
            vs30_metadata['units'] = units
            vs30_metadata['digits'] = digits
            oc.setGrid('vs30', vs30, metadata=vs30_metadata)

            for key, value in outgrid.items():
                component = config['interp']['component']
                # set the data grid
                mean_grid = Grid2D(value, gdict)
                mean_layername, units, digits = get_layer_info(key)
                mean_metadata = {'units': units,
                                 'digits': digits}

                # set the uncertainty grid
                std_layername, units, digits = get_layer_info(key+'_sd')
                std_metadata = {'units': units,
                                'digits': digits}
                std_grid = Grid2D(outsd[key], gdict.copy())
                oc.setIMT(key,
                          mean_grid, mean_metadata,
                          std_grid, std_metadata,
                          component)

        else:
            metadata['type'] = 'points'
            metadata['lons'] = lons.flatten()
            metadata['lats'] = lats.flatten()
            metadata['facility_ids'] = [x.encode('ascii') for x in idents]
            oc.setArray('vs30', vs30.flatten(), metadata=metadata)

        oc.close()

        self.logger.info('done')
        # ------------------------------------------------------------------
        # End interp()
        # ------------------------------------------------------------------


def get_period_index_from_imt_str(imtstr, imt_per_ix):
        """
        Get the index for the period of the specified IMT.

        Args:
            imtstr (str): The (OQ-style) IMT string.
            imt_per_ix (dict): Dictionary relating periods to their
                indices.

        Returns:
            int: The index corresponding to the period of the IMT.
        """
        if imtstr == 'PGA':
            return imt_per_ix['0.01']
        elif imtstr in ('PGV', 'MMI'):
            return imt_per_ix['1.0']
        else:
            return imt_per_ix[imtstr.replace('SA(', '').replace(')', '')]


def get_period_array(*args):
    """
    Return an array of the periods represented by the IMT list(s) in
    the input.

    Args:
        *args (list): One or more lists of IMTs.

    Returns:
        array: Numpy array of the (sorted) periods represented by the
        IMTs in the input list(s).
    """
    imt_per = set()
    for imt_list in args:
        if imt_list is None:
            continue
        for imtstr in imt_list:
            if imtstr == 'PGA':
                imt_per.add(0.01)
            elif imtstr == 'PGV' or imtstr == 'MMI':
                imt_per.add(1.0)
            else:
                imt_per.add(float(imtstr.replace('SA(', '').replace(')', '')))
    return np.array(sorted(imt_per))


def get_imts(imtstr, imtset):
    """
    Return the input imt or its closest surrogarte (or bracket) found
    in imtset.

    Args:
        imtstr (str): An (OQ-style) IMT string.
        imtsset (set): A set of IMTs to search for imtstr or its closest
            surrogate (or bracket).

    Returns:
        tuple: The IMT, it's closest surrogate, or a bracket of periods on
        either side of the IMT's period, from the IMTs in intset.
    """

    if imtstr in imtset:
        return (imtstr, )

    salist = [x for x in imtset if x.startswith('SA(')]
    periodlist = [float(x.replace('SA(', '').replace(')', '')) for x in salist]
    periodlist = sorted(periodlist)
    periodlist_str = [str(x) for x in periodlist]

    #
    # If we're here, then we know that IMT isn't in the inputs. Try
    # some alternatives.
    #
    if imtstr == 'PGA':
        #
        # Use the highest frequency in the inputs, otherwise use PGV
        #
        if len(salist):
            return ('SA(' + periodlist_str[0] + ')', )
        elif 'PGV' in imtset:
            return ('PGV', )
        else:
            return ()
    elif imtstr == 'PGV':
        #
        # Use 1.0 sec SA (or its bracket) if it's there, otherwise
        # use PGA
        #
        sa_tuple = get_sa_bracket(1.0, periodlist, periodlist_str)
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
        return get_sa_bracket(1.0, periodlist, periodlist_str)
    elif imtstr.startswith('SA('):
        myper = float(imtstr.replace('SA(', '').replace(')', ''))
        return get_sa_bracket(myper, periodlist, periodlist_str)
    else:
        raise ValueError('Unknown IMT %s in get_imt_bracket' % imtstr)


def get_sa_bracket(myper, plist, plist_str):

    """
    For a given period, look through the input periods and return a tuple of
    a) the single IMT string representing the period itself if it is found
    in the input; b) a pair of IMT strings representing the periods
    bracketing the given period; or c) the single IMT representing the first
    or last period in the input list if the given period is off the end of
    the list.

    Args:
        myper (float): The period to search for in the input lists.
        plist (array): A sorted list of periods as floats.
        plist_str (array): The periods in 'plist' (above) as strings.

    Returns:
        tuple: One or two strings representing the IMTs closest to or
        bracketing the input period.

    """
    if not len(plist):
        return ()
    try:
        return ('SA(' + plist_str[plist.index(myper)] + ')', )
    except ValueError:
        pass
    for i, v in enumerate(plist):
        if v > myper:
            break
    if i == 0 or v < myper:
        return ('SA(' + plist_str[i] + ')', )
    else:
        return ('SA(' + plist_str[i-1] + ')', 'SA(' + plist_str[i] + ')')


def get_sta_imts(imtstr, sdf, ix, imtset):
    """
    Get the IMT, its closest surrogate, or its bracket for a stataion
    in a particular station dataframe, accounting for missing or
    flagged data.

    Args:
        imtstr (str): The desired IMT as an OQ-style string.

        sdf (dict): The dataframe containing the station.

        ix (int): The index of the station in the dataframe.

        imtset (set): The list of IMTs (as OQ-style strings) in the
            dataframe.

    Returns:
        tuple: The IMT, its closest surrogate, or the pair of IMTs
        bracketing it in period, gathered from the valid data for the
        station.
    """
    myimts = set()
    for this_imt in imtset:
        if not np.isnan(sdf[this_imt][ix]) and \
           not sdf[this_imt + '_outliers'][ix]:
            myimts.add(this_imt)
    return get_imts(imtstr, myimts)


def get_map_grade(do_grid, outsd, psd, moutgrid):
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
    if not do_grid or 'PGA' not in outsd or 'PGA' not in psd:
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
def get_layer_info(layer):
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
        layer_out = oq_to_file(layer.replace('_sd',''))
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


#
# Helper function to call get_mean_and_stddevs for the
# appropriate object given the IMT
#
def gmas(ipe, gmpe, sx, rx, dx, oqimt, stddev_types):
    """
    This is a helper function to call get_mean_and_stddevs for the
    appropriate object given the IMT.

    Args:
        ipe: An IPE instance.
        gmpe: A GMPE instance.
        sx: Sites context.
        rx: Rupture context.
        dx: Distance context.
        oqimt: List of OpenQuake IMTs.
        stddev_types: list of OpenQuake standard deviation types.

    Returns:
        tuple: Tuple of two items:

            - Numpy array of means,
            - List of numpy array of standard deviations corresponding to the
              requested stddev_types.

    """
    if 'MMI' in oqimt:
        pe = ipe
    else:
        pe = gmpe
    return pe.get_mean_and_stddevs(sx, rx, dx, oqimt, stddev_types)


#
# This is a stand-in for tau when the gmpe set doesn't provide it
# It is an educated guess based on the NGA-east work and BC Hydro
# gmpe. It's not perfect, but probably isn't too far off. It is
# only used when the GMPEs don't provide a breakdown of the
# uncertainty terms.
#
def get_default_std_inter():
    """
    This is a stand-in for tau when the gmpe set doesn't provide it
    It is an educated guess based on the NGA-east work and BC Hydro
    gmpe. It's not perfect, but probably isn't too far off. It is
    only used when the GMPEs don't provide a breakdown of the
    uncertainty terms.
    """
    return 0.35
