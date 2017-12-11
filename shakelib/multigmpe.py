#!/usr/bin/env python

import copy
from importlib import import_module

import numpy as np

from openquake.hazardlib.gsim.base import GMPE
from openquake.hazardlib.imt import PGA, PGV, SA
from openquake.hazardlib import const

from shakelib.conversions.imt.newmark_hall_1982 import NewmarkHall1982
from shakelib.conversions.imc.beyer_bommer_2006 import BeyerBommer2006
from shakelib.sites import Sites


class MultiGMPE(GMPE):
    """
    Implements a GMPE that is the combination of multiple GMPEs.

    To do

        * Allow site to be based on a model that isn't a GMPE (e.g.,
          Borcherdt).

    """

    DEFINED_FOR_TECTONIC_REGION_TYPE = None
    DEFINED_FOR_INTENSITY_MEASURE_TYPES = None
    DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = None
    DEFINED_FOR_STANDARD_DEVIATION_TYPES = None
    REQUIRES_SITES_PARAMETERS = None
    REQUIRES_RUPTURE_PARAMETERS = None
    REQUIRES_DISTANCES = None

    def get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types):
        """
        See superclass `method <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#openquake.hazardlib.gsim.base.GroundShakingIntensityModel.get_mean_and_stddevs>`__.
        """  # noqa

        # Evaluate MultiGMPE:
        lnmu, lnsd = self.__get_mean_and_stddevs(
            sites, rup, dists, imt, stddev_types)

        # Check for large-distance cutoff/weights
        if hasattr(self, 'CUTOFF_DISTANCE'):
            lnmu_large, lnsd_large = self.__get_mean_and_stddevs(
                sites, rup, dists, imt, stddev_types, large_dist=True)
            # Stomp on lnmu and lnsd at large distances
            dist_cutoff = self.CUTOFF_DISTANCE
            lnmu[dists.rjb > dist_cutoff] = lnmu_large[dists.rjb > dist_cutoff]
            for i in range(len(lnsd)):
                lnsd[i][dists.rjb > dist_cutoff] = \
                    lnsd_large[i][dists.rjb > dist_cutoff]

        return lnmu, lnsd

    def __get_mean_and_stddevs(self, sites, rup, dists, imt, stddev_types,
                               large_dist=False):

        # ---------------------------------------------------------------------
        # Sort out which set of weights to use
        # ---------------------------------------------------------------------
        if large_dist is False:
            wts = self.WEIGHTS
        else:
            wts = self.WEIGHTS_LARGE_DISTANCE

        # ---------------------------------------------------------------------
        # Sort out shapes of sites and dists elements
        # ---------------------------------------------------------------------

        shapes = []
        for k, v in sites.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                shapes.append(v.shape)
        for k, v in dists.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                shapes.append(v.shape)

        shapeset = set(shapes)
        if len(shapeset) != 1:
            raise Exception(
                'All sites and dists elements must have same shape.')
        else:
            orig_shape = list(shapeset)[0]

        # Need to turn all 2D arrays into 1D arrays because of
        # inconsistencies in how arrays are handled in OpenQuake.
        for k, v in dists.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                dists.__dict__[k] = np.reshape(dists.__dict__[k], (-1,))
        for k, v in sites.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                sites.__dict__[k] = np.reshape(sites.__dict__[k], (-1,))

        # ---------------------------------------------------------------------
        # These are arrays to hold the weighted combination of the GMPEs
        # ---------------------------------------------------------------------
        lnmu = np.zeros_like(sites.vs30)
        sd_avail = self.DEFINED_FOR_STANDARD_DEVIATION_TYPES
        if not sd_avail.issuperset(set(stddev_types)):
            raise Exception("Requested an unavailable stddev_type.")

        lnsd2 = [np.zeros_like(sites.vs30) for a in stddev_types]

        for i in range(len(self.GMPES)):
            # -----------------------------------------------------------------
            # Loop over GMPE list
            # -----------------------------------------------------------------

            gmpe = self.GMPES[i]

            sites = MultiGMPE.set_sites_depth_parameters(sites, gmpe)

            # -----------------------------------------------------------------
            # Evaluate GMPEs
            # -----------------------------------------------------------------

            gmpe_imts = [imt.__name__ for imt in
                         gmpe.DEFINED_FOR_INTENSITY_MEASURE_TYPES]
            if (isinstance(imt, PGV)) and ("PGV" not in gmpe_imts):
                # -------------------------------------------------------------
                # If IMT is PGV and PGV is not given by the GMPE, then
                # convert from PSA10.
                # -------------------------------------------------------------
                if self.HAS_SITE[i] is True:
                    psa10, psa10sd = gmpe.get_mean_and_stddevs(
                        sites, rup, dists, SA(1.0), stddev_types)
                else:
                    lamps = self.get_site_factors(
                        sites, rup, dists, SA(1.0), default=True)
                    psa10, psa10sd = gmpe.get_mean_and_stddevs(
                        sites, rup, dists, SA(1.0), stddev_types)
                    psa10 = psa10 + lamps

                lmean, lsd = NewmarkHall1982.psa102pgv(psa10, psa10sd[0])
            else:
                if self.HAS_SITE[i] is True:
                    lmean, lsd = gmpe.get_mean_and_stddevs(
                        sites, rup, dists, imt, stddev_types)
                else:
                    lamps = self.get_site_factors(
                        sites, rup, dists, imt, default=True)
                    lmean, lsd = gmpe.get_mean_and_stddevs(
                        sites, rup, dists, imt, stddev_types)
                    lmean = lmean + lamps

            # -----------------------------------------------------------------
            # Convertions due to component definition
            # -----------------------------------------------------------------

            imc_in = gmpe.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT
            imc_out = self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT
            lmean = BeyerBommer2006.ampIMCtoIMC(lmean, imc_in, imc_out, imt)
            for j in range(len(lnsd2)):
                lsd[j] = BeyerBommer2006.sigmaIMCtoIMC(
                    lsd[j], imc_in, imc_out, imt)

            # -----------------------------------------------------------------
            # Compute weighted mean and sd
            # -----------------------------------------------------------------

            lnmu = lnmu + wts[i] * lmean

            # Note: the lnsd2 calculation isn't complete until we drop out of
            # this loop and substract lnmu**2
            for j in range(len(lnsd2)):
                lnsd2[j] = lnsd2[j] + wts[i] * (lmean**2 + lsd[j]**2)

        for j in range(len(lnsd2)):
            lnsd2[j] = lnsd2[j] - lnmu**2

        lnsd = [np.sqrt(a) for a in lnsd2]

        # Undo reshapes of inputs
        for k, v in dists.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                dists.__dict__[k] = np.reshape(dists.__dict__[k], orig_shape)
        for k, v in sites.__dict__.items():
            if (k is not 'lons') and (k is not 'lats'):
                sites.__dict__[k] = np.reshape(sites.__dict__[k], orig_shape)

        # Reshape output
        lnmu = np.reshape(lnmu, orig_shape)
        for i in range(len(lnsd)):
            lnsd[i] = np.reshape(lnsd[i], orig_shape)

        return lnmu, lnsd

    @classmethod
    def from_config(cls, conf, filter_imt=None, verbose=False):
        """
        Construct a MultiGMPE from a config file.

        Args:
            conf (dict): Dictionary of config options.
            filter_imt (IMT): An optional IMT to filter/reweight the GMPE list.
            verbose (bool): Print verbose output for debugging.

        Returns:
            MultiGMPE object.

        """
        IMC = conf['component_modules'][conf['interp']['component']]
        selected_gmpe = conf['modeling']['gmpe']

        if verbose is True:
            print('selected_gmpe: %s' % selected_gmpe)
            print('IMC: %s' % IMC)

        # ---------------------------------------------------------------------
        # Allow for selected_gmpe to be found in either conf['gmpe_sets'] or
        # conf['gmpe_modules'], if it is a GMPE set, then all entries must be
        # either a GMPE or a GMPE set (cannot have a GMPE set that is a mix of
        # GMPEs and GMPE sets).
        # ---------------------------------------------------------------------

        if selected_gmpe in conf['gmpe_sets'].keys():
            selected_gmpe_sets = conf['gmpe_sets'][selected_gmpe]['gmpes']
            gmpe_set_weights = \
                [float(w) for w in conf['gmpe_sets'][selected_gmpe]['weights']]
            if verbose is True:
                print('selected_gmpe_sets: %s' % selected_gmpe_sets)
                print('gmpe_set_weights: %s' % gmpe_set_weights)

            # -----------------------------------------------------------------
            # If it is a GMPE set, does it contain GMPEs or GMPE sets?
            # -----------------------------------------------------------------

            set_of_gmpes = all([s in conf['gmpe_modules'] for s in
                                selected_gmpe_sets])
            set_of_sets = all([s in conf['gmpe_sets'] for s in
                               selected_gmpe_sets])

            if set_of_sets is True:
                mgmpes = []
                for s in selected_gmpe_sets:
                    mgmpes.append(cls.__multigmpe_from_gmpe_set(
                        conf, s, filter_imt=filter_imt, verbose=verbose))
                out = MultiGMPE.from_list(mgmpes, gmpe_set_weights, imc=IMC)
            elif set_of_gmpes is True:
                out = cls.__multigmpe_from_gmpe_set(
                    conf,
                    selected_gmpe,
                    filter_imt=filter_imt,
                    verbose=verbose)
            else:
                raise Exception("%s must consist exclusively of keys in "
                                "conf['gmpe_modules'] or conf['gmpe_sets']"
                                % selected_gmpe)
        elif selected_gmpe in conf['gmpe_modules'].keys():
            modinfo = conf['gmpe_modules'][selected_gmpe]
            mod = import_module(modinfo[1])
            tmpclass = getattr(mod, modinfo[0])
            out = MultiGMPE.from_list([tmpclass()], [1.0], imc=IMC)
        else:
            raise Exception("conf['modeling']['gmpe'] must be a key in "
                            "conf['gmpe_modules'] or conf['gmpe_sets']")

        out.DESCRIPTION = selected_gmpe
        return out

    def __multigmpe_from_gmpe_set(conf, set_name, filter_imt=None,
                                 verbose=False):
        """
        Private method for constructing a MultiGMPE from a set_name.

        Args:
            conf (ConfigObj): A ShakeMap config object.
            filter_imt (IMT): An optional IMT to filter/reweight the GMPE list.
            set_name (str): Set name; must correspond to a key in
                conf['set_name'].

        Returns:
            MultiGMPE.

        """
        IMC = conf['component_modules'][conf['interp']['component']]

        selected_gmpes = conf['gmpe_sets'][set_name]['gmpes']
        selected_gmpe_weights = \
            [float(w) for w in conf['gmpe_sets'][set_name]['weights']]

        # Check for large distance GMPEs
        if 'weights_large_dist' in conf['gmpe_sets'][set_name].keys():
            if not conf['gmpe_sets'][set_name]['weights_large_dist']:
                selected_weights_large_dist = None
            else:
                selected_weights_large_dist = \
                    [float(w) for w in
                     conf['gmpe_sets'][set_name]['weights_large_dist']]
        else:
            selected_weights_large_dist = None

        if 'dist_cutoff' in conf['gmpe_sets'][set_name].keys():
            if np.isnan(conf['gmpe_sets'][set_name]['dist_cutoff']):
                selected_dist_cutoff = None
            else:
                selected_dist_cutoff = \
                    float(conf['gmpe_sets'][set_name]['dist_cutoff'])
        else:
            selected_dist_cutoff = None

        if 'site_gmpes' in conf['gmpe_sets'][set_name].keys():
            if not conf['gmpe_sets'][set_name]['site_gmpes']:
                selected_site_gmpes = None
            else:
                selected_site_gmpes = \
                    conf['gmpe_sets'][set_name]['site_gmpes']
        else:
            selected_site_gmpes = None

        if 'weights_site_gmpes' in conf['gmpe_sets'][set_name].keys():
            if not conf['gmpe_sets'][set_name]['weights_site_gmpes']:
                selected_weights_site_gmpes = None
            else:
                selected_weights_site_gmpes = \
                    conf['gmpe_sets'][set_name]['weights_site_gmpes']
        else:
            selected_weights_site_gmpes = None

        # ---------------------------------------------------------------------
        # Import GMPE modules and initialize classes into list
        # ---------------------------------------------------------------------
        gmpes = []
        for g in selected_gmpes:
            mod = import_module(conf['gmpe_modules'][g][1])
            tmpclass = getattr(mod, conf['gmpe_modules'][g][0])
            gmpes.append(tmpclass())

        # ---------------------------------------------------------------------
        # Filter out GMPEs not applicable to this period
        # ---------------------------------------------------------------------
        if filter_imt is not None:
            filtered_gmpes, filtered_wts = filter_gmpe_list(
                gmpes, selected_gmpe_weights, filter_imt)
        else:
            filtered_gmpes, filtered_wts = gmpes, selected_gmpe_weights

        # ---------------------------------------------------------------------
        # Import site GMPEs
        # ---------------------------------------------------------------------
        if selected_site_gmpes is not None:
            if isinstance(selected_site_gmpes, str):
                selected_site_gmpes = [selected_site_gmpes]
            site_gmpes = []
            for g in selected_site_gmpes:
                mod = import_module(conf['gmpe_modules'][g][1])
                tmpclass = getattr(mod, conf['gmpe_modules'][g][0])
                site_gmpes.append(tmpclass())
        else:
            site_gmpes = None

        # ---------------------------------------------------------------------
        # Filter out site GMPEs not applicable to this period
        # ---------------------------------------------------------------------
        if site_gmpes is not None:
            if filter_imt is not None:
                filtered_site_gmpes, filtered_site_wts = filter_gmpe_list(
                    site_gmpes, selected_weights_site_gmpes, filter_imt)
            else:
                filtered_site_gmpes = copy.copy(site_gmpes)
                filtered_site_wts = copy.copy(selected_weights_site_gmpes)
        else:
            filtered_site_gmpes = None
            filtered_site_wts = None

        # ---------------------------------------------------------------------
        # Construct MultiGMPE
        # ---------------------------------------------------------------------
        if verbose is True:
            print('    filtered_gmpes: %s' % filtered_gmpes)
            print('    filtered_wts: %s' % filtered_wts)

        mgmpe = MultiGMPE.from_list(
            filtered_gmpes, filtered_wts,
            default_gmpes_for_site=filtered_site_gmpes,
            default_gmpes_for_site_weights=filtered_site_wts,
            imc=IMC)

        # ---------------------------------------------------------------------
        # Append large-distance info if specified
        # ---------------------------------------------------------------------
        if selected_dist_cutoff is not None:
            if filter_imt is not None:
                filtered_gmpes_ld, filtered_wts_ld = filter_gmpe_list(
                    gmpes, selected_weights_large_dist, filter_imt)
            else:
                filtered_wts_ld = copy.copy(selected_weights_large_dist)

            mgmpe.CUTOFF_DISTANCE = copy.copy(selected_dist_cutoff)
            mgmpe.WEIGHTS_LARGE_DISTANCE = copy.copy(filtered_wts_ld)

        return mgmpe

    @classmethod
    def from_list(cls, gmpes, weights,
                  imc=const.IMC.GREATER_OF_TWO_HORIZONTAL,
                  default_gmpes_for_site=None,
                  default_gmpes_for_site_weights=None,
                  reference_vs30=760):
        """
        Construct a MultiGMPE instance from lists of GMPEs and weights.

        Args:
            gmpes (list): List of OpenQuake
                `GMPE <http://docs.openquake.org/oq-hazardlib/master/gsim/index.html#built-in-gsims>`__
                instances.

            weights (list): List of weights; must sum to 1.0.

            imc: Requested intensity measure component. Must be one listed
                `here <http://docs.openquake.org/oq-hazardlib/master/const.html?highlight=imc#openquake.hazardlib.const.IMC>`__.
                The amplitudes returned by the GMPEs will be converted to this
                IMT. Default is 'GREATER_OF_TWO_HORIZONTAL', which is used by
                ShakeMap. See discussion in
                `this section <http://usgs.github.io/shakemap/tg_choice_of_parameters.html#use-of-peak-values-rather-than-mean>`__
                of the ShakeMap manual.

            default_gmpes_for_site (list):
                Optional list of OpenQuake GMPE instance to use as a site term
                for any of the GMPEs that do not have a site term.

                Notes:

                    * We do not check for consistency in the reference rock
                      defintion, so the user nees to be aware of this issue and
                      holds responsibiilty for ensuring compatibility.
                    * We check whether or not a GMPE has a site term by c
                      hecking the REQUIRES_SITES_PARAMETERS slot for vs30.

            default_gmpes_for_site_weights: Weights for default_gmpes_for_site.
                Must sum to one and be same length as default_gmpes_for_site.
                If None, then weights are set to be equal.

            reference_vs30:
                Reference rock Vs30 in m/s. We do not check that this matches
                the reference rock in the GMPEs so this is the responsibility
                of the user.

        """  # noqa

        # ---------------------------------------------------------------------
        # Check that GMPE weights sum to 1.0:
        # ---------------------------------------------------------------------

        if np.abs(np.sum(weights) - 1.0) > 1e-7:
            raise Exception('Weights must sum to one.')

        # ---------------------------------------------------------------------
        # Check that length of GMPE weights equals length of gmpe list
        # ---------------------------------------------------------------------

        if len(weights) != len(gmpes):
            raise Exception(
                'Length of weights must match length of GMPE list.')

        # ---------------------------------------------------------------------
        # Check that gmpes is a list of OQ GMPE instances
        # ---------------------------------------------------------------------

        for g in gmpes:
            if not isinstance(g, GMPE):
                raise Exception("\"%s\" is not a GMPE instance." % g)

        self = cls()
        self.GMPES = gmpes
        self.WEIGHTS = weights

        # ---------------------------------------------------------------------
        # Combine the intensity measure types. This is problematic:
        #   - Logically, we should only include the intersection of the sets
        #     of imts for the different GMPEs.
        #   - In practice, this is not feasible because most GMPEs in CEUS and
        #     subduction zones do not have PGV.
        #   - So instead we will use the union of the imts and then convert
        #     to get the missing imts later in get_mean_and_stddevs.
        # ---------------------------------------------------------------------

        imts = [g.DEFINED_FOR_INTENSITY_MEASURE_TYPES for g in gmpes]
        self.DEFINED_FOR_INTENSITY_MEASURE_TYPES = set.union(*imts)

        # ---------------------------------------------------------------------
        # For VirtualIPE class, we also want to know if ALL of the GMPEs are
        # defined for PGV, in which case we will convert from PGV to MI,
        # otherwise use PGA or Sa.
        # ---------------------------------------------------------------------
        haspgv = [PGV in g.DEFINED_FOR_INTENSITY_MEASURE_TYPES for g in gmpes]
        self.ALL_GMPES_HAVE_PGV = all(haspgv)

        # ---------------------------------------------------------------------
        # Store intensity measure types for conversion in get_mean_and_stddevs.
        # ---------------------------------------------------------------------
        self.IMCs = [g.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT for g in gmpes]

        # ---------------------------------------------------------------------
        # Store the component
        # ---------------------------------------------------------------------
        self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT = imc

        # ---------------------------------------------------------------------
        # Intersection of GMPE standard deviation types
        # ---------------------------------------------------------------------
        stdlist = [set(g.DEFINED_FOR_STANDARD_DEVIATION_TYPES) for g in gmpes]
        self.DEFINED_FOR_STANDARD_DEVIATION_TYPES = \
            set.intersection(*stdlist)

        # ---------------------------------------------------------------------
        # Need union of site parameters, but it is complicated by the
        # different depth parameter flavors.
        # ---------------------------------------------------------------------
        sitepars = [g.REQUIRES_SITES_PARAMETERS for g in gmpes]
        self.REQUIRES_SITES_PARAMETERS = set.union(*sitepars)

        # ---------------------------------------------------------------------
        # Construct a list of whether or not each GMPE has a site term
        # ---------------------------------------------------------------------
        self.HAS_SITE = ['vs30' in g.REQUIRES_SITES_PARAMETERS for g in gmpes]

        # ---------------------------------------------------------------------
        # Checks and sort out defaults
        # ---------------------------------------------------------------------

        # things to check if default_gmpes_for_site is provided
        if default_gmpes_for_site is not None:
            # check that default_gmpe_for_site are OQ GMPEs or None
            for g in default_gmpes_for_site:
                if not isinstance(g, GMPE):
                    raise Exception("\"%s\" is not a GMPE instance." % g)

            # apply default weights if necessary
            if default_gmpes_for_site_weights is None:
                n = len(default_gmpes_for_site)
                default_gmpes_for_site_weights = [1 / n] * n

        # Things to check if one or more GMPE does not have a site term
        if not all(self.HAS_SITE):
            # Raise an exception if no default site is provided
            if default_gmpes_for_site is None:
                raise Exception('Must provide default_gmpes_for_site if one or'
                                ' more GMPE does not have site term.')

            # If weights are unspecified, use equal weight
            if default_gmpes_for_site_weights is None:
                default_gmpes_for_site_weights = \
                    [1 / len(default_gmpes_for_site)] * \
                    len(default_gmpes_for_site)

            # check that length of default_gmpe_for_site matches length of
            # default_gmpe_for_site_weights
            if len(default_gmpes_for_site_weights) != \
               len(default_gmpes_for_site):
                raise Exception('Length of default_gmpes_for_site_weights '
                                'must match length of default_gmpes_for_site '
                                'list.')

            # check weights sum to one if needed
            if not all(self.HAS_SITE):
                if np.sum(default_gmpes_for_site_weights) != 1.0:
                    raise Exception('default_gmpes_for_site_weights must sum'
                                    ' to one.')

        # Note: if ALL of the GMPEs do not have a site term (requiring Vs30),
        #       then REQUIRES_SITES_PARAMETERS for the MultiGMPE will not
        #       include Vs30 even though it will be needed to compute the
        #       default site term. So if the site checks have passed to this
        #       point, we should add Vs30 to the set of required site pars:
        self.REQUIRES_SITES_PARAMETERS = set.union(
            self.REQUIRES_SITES_PARAMETERS, set(['vs30']))

        self.DEFAULT_GMPES_FOR_SITE = default_gmpes_for_site
        self.DEFAULT_GMPES_FOR_SITE_WEIGHTS = default_gmpes_for_site_weights
        self.REFERENCE_VS30 = reference_vs30

        # ---------------------------------------------------------------------
        # Union of rupture parameters
        # ---------------------------------------------------------------------
        ruppars = [g.REQUIRES_RUPTURE_PARAMETERS for g in gmpes]
        self.REQUIRES_RUPTURE_PARAMETERS = set.union(*ruppars)

        # ---------------------------------------------------------------------
        # Union of distance parameters
        # ---------------------------------------------------------------------
        distpars = [g.REQUIRES_DISTANCES for g in gmpes]
        self.REQUIRES_DISTANCES = set.union(*distpars)

        return self

    def get_site_factors(self, sites, rup, dists, imt, default=False):
        """
        Method for computing site amplification factors from the defalut GMPE
        to be applied to GMPEs which do not have a site term.

        **NOTE** Amps are calculated in natural log units and so the ln(amp)
        is returned.

        Args:
            sites (SitesContext): Instance of SitesContext.
            rup (RuptureContext): Instance of RuptureContext.
            dists (DistancesContext): Instance of DistancesContext.
            imt: An instance openquake.hazardlib.imt.
            default (bool): Boolean of whether or not to return the
                amplificaiton factors for the gmpes or default_gmpes_for_site.
                This argument is primarily only intended to be used internally
                for when we just need to access the default amplifications to
                apply to those GMPEs that do not have site terms.

        Returns:
            Site amplifications in natural log units.
        """

        # ---------------------------------------------------------------------
        # Make reference sites context
        # ---------------------------------------------------------------------

        ref_sites = copy.deepcopy(sites)
        ref_sites.vs30 = np.ones_like(sites.vs30) * self.REFERENCE_VS30

        # ---------------------------------------------------------------------
        # If default True, construct new MultiGMPE with default GMPE/weights
        # ---------------------------------------------------------------------
        if default is True:
            tmp = MultiGMPE.from_list(
                self.DEFAULT_GMPES_FOR_SITE,
                self.DEFAULT_GMPES_FOR_SITE_WEIGHTS,
                self.DEFINED_FOR_INTENSITY_MEASURE_COMPONENT)

        # ---------------------------------------------------------------------
        # If default False, just use self
        # ---------------------------------------------------------------------
        else:
            tmp = self

        lmean, lsd = tmp.get_mean_and_stddevs(
                sites, rup, dists, imt,
                list(tmp.DEFINED_FOR_STANDARD_DEVIATION_TYPES))
        lmean_ref, lsd = tmp.get_mean_and_stddevs(
                ref_sites, rup, dists, imt,
                list(tmp.DEFINED_FOR_STANDARD_DEVIATION_TYPES))

        lamps = lmean - lmean_ref

        return lamps

    @staticmethod
    def set_sites_depth_parameters(sites, gmpe):
        """
        Need to select the appropriate z1pt0 value for different GMPEs.
        Note that these are required site parameters, so even though
        OQ has these equations built into the class in most cases.
        I have submitted an issue to OQ requesting subclasses of these
        methods that do not require the depth parameters in the
        SitesContext to make this easier.

        Args:
            sites:1 An OQ sites context.
            gmpe: An OQ GMPE instance.

        Returns:
            An OQ sites context with the depth parameters set for the
            requested GMPE.
        """

        sites = Sites._addDepthParameters(sites)

        if gmpe == 'AbrahamsonEtAl2014()':
            sites.z1pt0 = sites.z1pt0_ask14_cal
        if gmpe == 'ChiouYoungs2014()':
            # Also BooreEtAl2014() if using subclass with depth parameter
            sites.z1pt0 = sites.z1pt0_cy14_cal
        if gmpe == 'CampbellBozorgnia2014()':
            sites.z2pt5 = sites.z2pt5_cb14_cal
        if gmpe == 'ChiouYoungs2008()':
            sites.z1pt0 = sites.z1pt0_cy08
        if gmpe == 'CampbellBozorgnia2008()':
            sites.z2pt5 = sites.z2pt5_cb07

        return sites


def filter_gmpe_list(gmpes, wts, imt):
    """
    Method to remove GMPEs from the GMPE list that are not applicable
    to a specific IMT. Rescales the weights to sum to one.

    Args:
        gmpes (list): List of GMPE instances.
        wts (list): List of floats indicating the weight of the GMPEs.
        imt (IMT): OQ IMT to filter GMPE list for.

    Returns:
        tuple: List of GMPE instances and list of weights.

    """
    if wts is None:
        n = len(gmpes)
        wts = [1 / n] * n

    per_max = [np.max(get_gmpe_sa_periods(g)) for g in gmpes]
    per_min = [np.min(get_gmpe_sa_periods(g)) for g in gmpes]
    if imt == PGA():
        sgmpe = [g for g in gmpes if imt in
                 get_gmpe_coef_table(g).non_sa_coeffs]
        swts = [w for g, w in zip(gmpes, wts) if imt in
                get_gmpe_coef_table(g).non_sa_coeffs]
    elif imt == PGV():
        sgmpe = []
        swts = []
        for i in range(len(gmpes)):
            if (imt in get_gmpe_coef_table(gmpes[i]).non_sa_coeffs) or\
               (per_max[i] >= 1.0 and per_min[i] <= 1.0):
                sgmpe.append(gmpes[i])
                swts.append(wts[i])
    else:
        per = imt.period
        sgmpe = []
        swts = []
        for i in range(len(gmpes)):
            if (per_max[i] >= per and per_min[i] <= per):
                sgmpe.append(gmpes[i])
                swts.append(wts[i])

    if len(sgmpe) == 0:
        raise Exception('No applicable GMPEs from GMPE list for %s' % imt)

    # Scale weights to sum to one
    swts = np.array(swts)
    swts = swts / np.sum(swts)

    return sgmpe, swts


def get_gmpe_sa_periods(gmpe):
    """
    Method to extract the SA periods defined by a GMPE.

    Args:
        gmpe (GMPE): A GMPE instance.

    Retunrs:
        list: List of periods.

    """
    ctab = get_gmpe_coef_table(gmpe).sa_coeffs
    ilist = list(ctab.keys())
    per = [i.period for i in ilist]
    return per


def get_gmpe_coef_table(gmpe):
    """
    Method for finding the (or "a") GMPE table.

    Notes:

      *  The reason for the complexity here is that there can be multiple
         coefficient tables, and some of them may not have the sa_coeffs
         attribute, which is the main reason for getting the table.
      *  We are also assuming that if there are more than one  coefficient
         table, the range of periods will be the same across all of the
         tables.

    Args:
        gmpe (GMPE): An OQ GMPE instance.

    Returns:
        The associated coefficient table.

    """
    stuff = gmpe.__dir__()
    coef_list = [s for s in stuff if 'COEFFS' in s]
    for coef_sel in coef_list:
        cobj = getattr(gmpe, coef_sel)
        if "sa_coeffs" in cobj.__dir__():
            return cobj
    raise Exception("GMPE %s does not contain sa_coeffs attribute." % gmpe)
