"""
Module to simplify importing the OQ implementation of the
NGA-East GMPE suite.
"""

import os
import re

import pandas as pd
import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib import const
from openquake.hazardlib import imt as IMT
from openquake.hazardlib.gsim.usgs_ceus_2019 import NGAEastUSGSGMPE


class NGAEast(GMPE):
    """
    Returns NGA East GMPE that combines all of the individual NGAEastUSGSGMPE
    GMPEs.
    """
    DEFINED_FOR_TECTONIC_REGION_TYPE = const.TRT.STABLE_CONTINENTAL
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = const.IMC.RotD50
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = set([
        IMT.PGA,
        # IMT.PGV,
        IMT.SA
    ])
    # DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
    #     const.StdDev.TOTAL,
    # ])

    DEFINED_FOR_STANDARD_DEVIATION_TYPES = set([
        const.StdDev.TOTAL,
        const.StdDev.INTER_EVENT,
        const.StdDev.INTRA_EVENT
    ])
    REQUIRES_SITES_PARAMETERS = set(('vs30', ))
    REQUIRES_RUPTURE_PARAMETERS = set(('mag', ))
    REQUIRES_DISTANCES = set(('rrup', ))

    TABLE_PATHS = os.listdir(NGAEastUSGSGMPE.PATH)
    this_module = os.path.dirname(__file__)
    NGA_BASE_PATH = os.path.join(
        this_module, '..', '..', 'shakemap', 'data', 'nga_east_tables')

    NGA_EAST_USGS_WEIGHT = 0.667
    NGA_EAST_SEED_WEIGHT = 0.333

    NGA_EAST_USGS = pd.read_csv(os.path.join(
        NGA_BASE_PATH, 'nga-east-usgs-weights.dat'))
    NGA_EAST_SEEDS = pd.read_csv(os.path.join(
        NGA_BASE_PATH, 'nga-east-seed-weights.dat'))

    # Sigma models and their weights
    SIGMA_MODS = ["EPRI", "PANEL"]
    SIGMA_WEIGHTS = [0.8, 0.2]

    # -------------------------------------------------------------------------
    # To simplify, use the COLLAPSED branch, but cannot get inter and intra
    # event standard deviations in this case.
    # SIGMA_MODS = ["COLLAPSED"]
    # SIGMA_WEIGHTS = [1.0]

    # Parse the periods of the columns
    per_idx_start = 1
    per_idx_end = -2
    per_list_str = NGA_EAST_USGS.keys().tolist()[per_idx_start:per_idx_end]
    per_array = np.array(
        [float(p.replace('SA', '').replace('P', '.')) for p in per_list_str]
    )

    def __init__(self):
        gmpes = []
        sigma_wts = []
        all_table_paths = []
        for i, sigma_mod in enumerate(self.SIGMA_MODS):
            for table_path in self.TABLE_PATHS:
                gmpe = NGAEastUSGSGMPE(
                    gmpe_table=table_path, sigma_model=sigma_mod)
                gmpes.append(gmpe)
                sigma_wts.append(self.SIGMA_WEIGHTS[i])
                all_table_paths.append(table_path)
        self.gmpes = gmpes
        self.sigma_weights = np.array(sigma_wts)
        self.ALL_TABLE_PATHS = all_table_paths

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        # List of GMPE weights, which is the product of the the branch weights
        # for the seed models vs the NGA East resampled models as well as the
        # weights for the indivudual GMPES as defined by Petersen et al. (2019)
        #
        # Note that the NGA East resampled models are a function of spectral
        # period.
        #
        # NGA East Seeds (1/3)
        #     ├── B_bca10d (0.06633), wts = 0.333 * 0.06633 = 0.02208789
        #     ├── B_ab95 (0.02211), wts = 0.333 * 0.02211 = 0.00736263
        #     ...
        # NGA East Resampled or "USGS" (2/3)
        #     ├── Model 1 (0.1009 for PGA), wts = 0.667 * 0.1009 = 0.0673003
        #     ├── Model 2 (0.1606 for PGA), wts = 0.667 * 0.1606 = 0.1071202
        #     ...
        #
        wts = [0] * len(self.gmpes)

        for i, tp in enumerate(self.ALL_TABLE_PATHS):
            if 'usgs' in tp:
                # Get model number from i-th path using regex
                mod_num = int(re.search(r'\d+', tp).group())
                coefs = np.array(
                    self.NGA_EAST_USGS.iloc[mod_num - 1]
                )
                # Is the IMT PGA, PGA, or SA?
                if imt == IMT.PGA():
                    iweight = coefs[-2]
                elif imt == IMT.PGV():
                    iweight = coefs[-1]
                else:
                    # For SA, need to interpolate; we'll use log-period and
                    # linear-weight interpolation.
                    iweight = np.interp(
                        np.log(imt.period),
                        np.log(self.per_array),
                        coefs[self.per_idx_start:self.per_idx_end]
                    )
                wts[i] = self.NGA_EAST_USGS_WEIGHT * iweight
            else:
                # Strip off the cruft to get the string we need to match
                str_match = tp.replace('nga_east_', '').replace('.hdf5', '')
                matched = self.NGA_EAST_SEEDS[
                    self.NGA_EAST_SEEDS['model'] == str_match]
                if len(matched):
                    iweight = self.NGA_EAST_SEEDS[
                        self.NGA_EAST_SEEDS['model'] == str_match].iloc[0, 1]
                    wts[i] = self.NGA_EAST_SEED_WEIGHT * iweight

        total_gmpe_weights = self.sigma_weights * wts

        if not np.allclose(np.sum(total_gmpe_weights), 1.0):
            raise ValueError('Weights must sum to 1.0.')

        mean = np.full_like(sites.vs30, 0)
        stddevs = []
        for i in range(len(stddev_types)):
            stddevs.append(np.full_like(sites.vs30, 0))
        for i, gm in enumerate(self.gmpes):
            tmean, tstddevs = gm.get_mean_and_stddevs(
                sites, rup, dists, imt, stddev_types)
            mean += tmean * total_gmpe_weights[i]
            for j, sd in enumerate(tstddevs):
                stddevs[j] += sd * total_gmpe_weights[i]

        return mean, stddevs
