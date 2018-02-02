"""
Module implements BooreKishida2017 class to convert between various
horizontal intensity measure components.
"""
import logging
import pkg_resources
import glob
import os.path

import numpy as np
import pandas as pd
from openquake.hazardlib import const
from openquake.hazardlib.const import IMC
from openquake.hazardlib.imt import PGA, PGV, SA


class BooreKishida2017(object):
    """
    This class implements the Boore and Kishida (2017) conversions for
    horizontal intensity measure components. 

    This class explicitly supports the following subset of the conversions
    provided by B&K (specified as OpenQuake IMCs):

        - RotD50 <=> GMRotI50
        - RotD50 <=> AVERAGE_HORIZONTAL (i.e., Geometric Mean "as recorded")
        - RotD100 <=> RotD50
        - RotD50 <=> GREATER_OF_TWO_HORIZONTAL
        - RotD100 <=> GREATER_OF_TWO_HORIZONTAL
        - GMRotI50 <=> GREATER_OF_TWO_HORIZONTAL
        - AVERAGE_HORIZONTAL <=> GREATER_OF_TWO_HORIZONTAL

    Not supported are IMCs for which there is no OpenQuake equivalent; 
    also not supported are chained conversions (that is, we do not 
    support conversions from A to C when conversions from A to B and 
    B to C are available; these must be done in two steps. For IMCs 
    not explicitly
    supported by B&K, we assume the IMC is equivalent to the geometric
    mean (which B&K call GM_AR), based on our reading of
    Beyer & Bommer (2006).

    References

        Boore, D.M. and T. Kishida (2017). Relations between some
        horizontal-component ground-motion intensity measures used
        in practice. Bulletin of the Seismological Society of 
        America, 107(1), 334-343, doi: 10.1785/0120160250.

        Beyer, K., & Bommer, J. J. (2006). Relationships between median values
        and between aleatory variabilities for different definitions of the
        horizontal component of motion. Bulletin of the Seismological Society
        of America, 96(4A), 1512-1522.
    """

    def __init__(self, imc_in, imc_out):
        """
        Args:
            imc_in (OpenQuake IMC): The IMC of the ground motions to be
                converted.
            imc_out (OpenQuake IMC): The IMC to which the ground motions
                are to be converted.
        """

        filename, forward = self._imcPairToFile(imc_in, imc_out)
        if filename is None:
            raise ValueError("Can't find a conversion file for %s and %s" %
                             (imc_in, imc_out))
        self.forward = forward
        if filename == 'Null':
            # Null conversion -- imc_in and imc_out are either identical
            # or at least functionally equivalent
            self.pars = None
        else:
            self.pars = pd.read_csv(filename)

    def convertAmps(self, imt, amps, rrups, mag):
        """
        Return an array of amps converted from one IMC to another.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            amps (array): A numpy array of the (logged) ground motions
                to be converted. 
            rrups (array): A numpy array of the same shape as amps, 
                containing the rupture distances of the ground motions.
            mag (float): The earthquake magnitude.

        Returns:
            array: A numpy array of converted ground motions (logged).

        Raises:
            ValueError: If the IMT is not an allowed type.
        """

        if self.pars is None:
            # imc_in and imc_out are the same
            return amps.copy()

        (sigma, c0, r1, m1, m2) = self._getParamsFromIMT(imt)

        rrups_clipped = np.clip(rrups, 1e-2, 400)
        if mag < 2:
            mag = 2.0
        elif mag > 9:
            mag = 9.0
        ln_ratio = c0 + r1 * np.log(rrups_clipped / 50) + \
            m1 * (mag - 5.5) + m2 * (mag - 5.5)**2
        #
        # The B&K file naming convention has things like D100D50, which
        # means the parameters give the (log) ratio of RotD100/RotD50,
        # but we use the convention that RotD100 would be the input IMC
        # and RotD50 would be the output IMC, so we reverse the sense
        # of the conversion here.
        #
        if self.forward:
            amps = amps - ln_ratio
        else:
            amps = amps + ln_ratio

        return amps

    def convertStddevs(self, imt, stddevs, rrups, mag):
        """
        Return an array of standard deviations converted from one IMC 
        to another.

        Note that the action of this method is to always increase the
        input standard deviations. Thus, while converting from one IMC
        to another and then back again will yield the original ground
        motions via convertAmps(), the standard deviations will be 
        inflated by both conversions via this method.

        Args:
            imt (OpenQuake IMT): The intensity measure type of the input
                ground motions. Valid IMTs are PGA, PGV, and SA.
            stddevs (array): A numpy array of the standard deviations of
                the logged ground motions.
            rrups (array): Ignored by this method.
            mag (float): Ignored by this method.

        Returns:
            array: A numpy array of converted standard deviations.

        Raises:
            ValueError: If the IMT is not an allowed type.
        """

        if self.pars is None:
            # imc_in and imc_out are the same
            return stddevs.copy()

        (sigma, c0, r1, m1, m2) = self._getParamsFromIMT(imt)

        stddevs = np.sqrt(stddevs**2 + sigma**2)

        return stddevs

    def _getParamsFromIMT(self, imt):
        """
        Helper function to return (possibly interpolated) conversion
        parameters for a given IMT.
        """

        if imt == PGA():
            sigma = self.pars['sigma'][0]
            c0 = self.pars['c0smooth'][0]
            r1 = self.pars['r1smooth'][0]
            m1 = self.pars['m1smooth'][0]
            m2 = self.pars['m2smooth'][0]
        elif imt == PGV():
            sigma = self.pars['sigma'][1]
            c0 = self.pars['c0smooth'][1]
            r1 = self.pars['r1smooth'][1]
            m1 = self.pars['m1smooth'][1]
            m2 = self.pars['m2smooth'][1]
        elif 'SA' in imt:
            imt_per = imt.period
            pa = self.pars['per'][2:]
            sigma = np.interp(imt_per, pa, self.pars['sigma'][2:])
            c0 = np.interp(imt_per, pa, self.pars['c0smooth'][2:])
            r1 = np.interp(imt_per, pa, self.pars['r1smooth'][2:])
            m1 = np.interp(imt_per, pa, self.pars['m1smooth'][2:])
            m2 = np.interp(imt_per, pa, self.pars['m2smooth'][2:])
        else:
            raise ValueError("Unknown IMT: %s" % str(imt))

        return (sigma, c0, r1, m1, m2)

    @staticmethod
    def _imcPairToFile(imc_in, imc_out):
        """
        Helper function to find the name of the file representing
        the conversion.

        Returns:
            (str, bool): The filename and a boolean 'forward' indicating
            whether the conversion should be done in the forward (True)
            or inverse (False) direction. If filename is None, then no
            appropriate conversion file coule be found; if it is the
            string 'Null', then imc_in and imc_out evaluated to be the
            same.
        """

        datadir = pkg_resources.resource_filename('shakelib.conversions.imc',
                                                  'data')
        conv_files = glob.glob(os.path.join(datadir, '*.csv'))
        stub1 = BooreKishida2017._imcToFilestr(imc_in)
        stub2 = BooreKishida2017._imcToFilestr(imc_out)
        if stub1 == stub2:
            # No conversion necessary
            return ('Null', True)
        #
        # Look for the conversion from imc_in -> imc_out
        #
        stub = stub1 + stub2
        filelist = [name for name in conv_files if stub in name]
        if len(filelist) == 1:
            return (filelist[0], True)
        #
        # Now try the conversion from imc_out -> imc_in
        #
        stub = stub2 + stub1
        filelist = [name for name in conv_files if stub in name]
        if len(filelist) == 1:
            return (filelist[0], False)
        #
        # Can't find anything
        #
        return (None, None)

    @staticmethod
    def _imcToFilestr(oq_imc):
        """
        Helper function to convert an OpenQuake IMC into part of the
        Boore & Kishida file name.
        """

        if oq_imc == IMC.RotD50:
            return 'D50'
        elif oq_imc == IMC.RotD100:
            return 'D100'
        elif oq_imc == IMC.GMRotI50:
            return 'GM50'
        elif oq_imc == IMC.AVERAGE_HORIZONTAL or \
                oq_imc == IMC.HORIZONTAL or \
                oq_imc == IMC.RANDOM_HORIZONTAL or \
                oq_imc == IMC.MEDIAN_HORIZONTAL:
            return 'GMAR'
        elif oq_imc == IMC.GREATER_OF_TWO_HORIZONTAL:
            return 'Larger'
        else:
            #
            # For less common IMCs, Beyer & Bommer (2006) found most
            # of them to be more or less equivalent to geometric mean
            #
            logging.warn("Can't handle IMC %s, using GMAR" % oq_imc)
            return 'GMAR'
