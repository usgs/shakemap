"""
Module implements BeyerBommer2006 class to convert between various
horizontal intensity measure components.
"""
# Third party imports
from openquake.hazardlib import const
import numpy as np

# Local imports
from shakelib.conversions.convert_imc import ComponentConverter


class BeyerBommer2006(ComponentConverter):
    """
    Implements conversion for various "Intensity Measure
    Components" (IMCs) per Beyer and Bommer (2006).

    IMC equivalencies:

    +---------------------------+------------------------+
    | OpenQuake                 | Beyer & Bommer         |
    +===========================+========================+
    | AVERAGE_HORIZONTAL        | GMxy (Geometric mean)  |
    +---------------------------+------------------------+
    | HORIZONTAL                | Random?                |
    +---------------------------+------------------------+
    | MEDIAN_HORIZONTAL         | AMxy (Arithmetic mean) |
    +---------------------------+------------------------+
    | GMRotI50                  | GMRotI50               |
    +---------------------------+------------------------+
    | RotD50                    | GMRotD50               |
    +---------------------------+------------------------+
    | RANDOM_HORIZONTAL         | Random                 |
    +---------------------------+------------------------+
    | GREATER_OF_TWO_HORIZONTAL | Env_xy                 |
    +---------------------------+------------------------+
    | VERTICAL                  | ---                    |
    +---------------------------+------------------------+

    Notes

        - AVERAGE_HORIZONTAL is the "reference" type.
        - The OQ IMC "HORIZONAL" indicates that the horizontal IMC category
          may be ambiguous. In these cases, we are assuming that it is a
          random horizontal component as a default.
        - Assumes ALL unknown IMC types are AVERAGE_HORIZONTAL

    References

        Beyer, K., & Bommer, J. J. (2006). Relationships between median values
        and between aleatory variabilities for different definitions of the
        horizontal component of motion. Bulletin of the Seismological Society
        of America, 96(4A), 1512-1522.
        `[link] <http://www.bssaonline.org/content/96/4A/1512.short>`__

    """
    def __init__(self, imc_in, imc_out):
        super().__init__()
        self.imc_in = imc_in
        self.imc_out = imc_out
        # c12 = median ratio
        # c34 = std. of log ratio
        # R = ratio of sigma_pga values
        self.__pga_pgv_col_names = ['c12', 'c34', 'R']
        self.__sa_col_names = ['c1', 'c2', 'c3', 'c4', 'R']
        # Possible conversions
        self.conversion_graph = {'Greater of two horizontal': set([
                'Median horizontal',
                'Average Horizontal (GMRotI50)',
                'Average Horizontal (RotD50)',
                'Random horizontal',
                'Horizontal',
                'Average horizontal']),
         'Median horizontal': set([
                'Greater of two horizontal',
                'Average Horizontal (GMRotI50)',
                'Average Horizontal (RotD50)',
                'Random horizontal',
                'Horizontal',
                'Average horizontal']),
         'Average Horizontal (GMRotI50)': set([
                'Greater of two horizontal',
                'Median horizontal',
                'Average Horizontal (RotD50)',
                'Random horizontal',
                'Horizontal',
                'Average horizontal']),
         'Average Horizontal (RotD50)': set([
                'Greater of two horizontal',
                'Median horizontal',
                'Average Horizontal (GMRotI50)',
                'Random horizontal',
                'Horizontal',
                'Average horizontal']),
         'Random horizontal': set([
                'Greater of two horizontal',
                'Median horizontal',
                'Average Horizontal (GMRotI50)',
                'Average Horizontal (RotD50)',
                'Horizontal',
                'Average horizontal']),
        'Horizontal': set([
               'Greater of two horizontal',
               'Median horizontal',
               'Average Horizontal (GMRotI50)',
               'Average Horizontal (RotD50)',
               'Random horizontal',
               'Average horizontal']),
        'Average horizontal': set([
               'Greater of two horizontal',
               'Median horizontal',
               'Average Horizontal (GMRotI50)',
               'Average Horizontal (RotD50)',
               'Random horizontal',
               'Horizontal']),}
        # Check if any imc values are unknown. If they are, convert
        # to AVERAGE_HORIZONTAL
        self.checkUnknown()
        # Get shortest conversion "path" between imc_in and imc_out
        self.path = self.getShortestPath(self.conversion_graph,
                self.imc_in, self.imc_out)
        # Set dictionaries of constants
        self.__pga_dict = {
            const.IMC.MEDIAN_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
            const.IMC.GMRotI50:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
            const.IMC.RotD50:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
            const.IMC.RANDOM_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
            const.IMC.HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
            const.IMC.GREATER_OF_TWO_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
        self.__pgv_dict = {
            const.IMC.MEDIAN_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.01, 1.00]))),
            const.IMC.GMRotI50:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
            const.IMC.RotD50:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.02, 1.00]))),
            const.IMC.RANDOM_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
            const.IMC.HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.0, 0.07, 1.03]))),
            const.IMC.GREATER_OF_TWO_HORIZONTAL:
                dict(list(zip(self.__pga_pgv_col_names, [1.1, 0.05, 1.02])))}
        self.__sa_dict = {
            const.IMC.MEDIAN_HORIZONTAL:
                dict(list(zip(self.__sa_col_names, [1.0, 1.0, 0.01, 0.02, 1.00]))),
            const.IMC.GMRotI50:
                dict(list(zip(self.__sa_col_names, [1.0, 1.0, 0.03, 0.04, 1.00]))),
            const.IMC.RotD50:
                dict(list(zip(self.__sa_col_names, [1.0, 1.0, 0.02, 0.03, 1.00]))),
            const.IMC.RANDOM_HORIZONTAL:
                dict(list(zip(self.__sa_col_names, [1.0, 1.0, 0.07, 0.11, 1.05]))),
            const.IMC.HORIZONTAL:
                dict(list(zip(self.__sa_col_names, [1.0, 1.0, 0.07, 0.11, 1.05]))),
            const.IMC.GREATER_OF_TWO_HORIZONTAL:
                dict(list(zip(self.__sa_col_names, [1.1, 1.2, 0.04, 0.07, 1.02])))}

    def convertAmpsOnce(self, imt, amps, rrups=None, mag=None):
        """
        Returns amps converted from one IMC to another.

        **Important**:

            - Assumes the input amps are in natural log (not linear) space
            - IMC type 'VERTICAL' is not supported

        Args:
            imt (IMT): OpenQuake IMT of the input amps (must be one of PGA,
                PGV, or SA).
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`

            amps (array): Numpy array of ground motion amplitudes.
            rrups (array): A numpy array of the same shape as amps,
                containing the rupture distances of the ground motions.
                Ignored by this method.
            mag (float): The earthquake magnitude. Default is None.
                 Ignored by this method.

        Returns:
            array: Numpy array of amps converted from imc_in to imc_out.
        """

        denom = self.__GM2other(imt, self.imc_in)
        numer = self.__GM2other(imt, self.imc_out)

        return amps + np.log(numer / denom)

    def convertSigmasOnce(self, imt, sigmas):
        """
        Returns standard deviations converted from one IMC to another.

        **Important**:

            - Assumes the input sigmas are in natural log space
            - IMC types 'VERTICAL' and 'HORIZONTAL' are not supported

        Args:
            imt (IMT): OpenQuake IMT of the input sigmas (must be one of PGA,
                 PGV, or SA) `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__
            sigmas (array): Numpy array of standard deviations.

        Returns:
            array: Numpy array of standard deviations converted from imc_in to
                imc_out.

        """
        # ---------------------------------------------------------------------
        # Take input sigma to geometric mean sigma.
        # Solve eqn 8 for sigma_GM2, which is sigma_logSa_GM
        # This is the sigma converted to geometric mean
        # (i.e., reference component).
        # ---------------------------------------------------------------------
        R, sig_log_ratio = self.__GM2otherSigma(imt, self.imc_in)
        sigma_GM2 = (sigmas**2 - sig_log_ratio**2) / R**2

        # ---------------------------------------------------------------------
        # Evaluate equation 8 to go from GM to requested IMC
        # ---------------------------------------------------------------------
        R, sig_log_ratio = self.__GM2otherSigma(imt, self.imc_out)
        sigma_out = np.sqrt(sigma_GM2 * R**2 + sig_log_ratio**2)

        return sigma_out

    def __GM2other(self, imt, imc):
        """
        Helper function to extract coefficients from the parameters for
        converting the median ground motions.

        Args:
            imt (IMT): OQ intensity measure type.
            imc (IMC): OQ Intensity measure component.

        Returns:
            float: Median ratios. This is directly from Table 2 for PGA and
                PGV, and computed from coefficients in Table 3 along with
                eqn 10 for Sa.

        Raises:
            ValueError if unknown IMT specified.
        """
        if imc == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            return 1.

        self._verifyConversion(imc)

        if 'PGA' in imt:
            return self.__pga_dict[imc]['c12']
        elif 'PGV' in imt:
            return self.__pgv_dict[imc]['c12']
        elif 'SA' in imt:
            pp = imt.period
            if pp <= 0.15:
                return self.__sa_dict[imc]['c1']
            elif pp < 0.8:
                c1 = self.__sa_dict[imc]['c1']
                c2 = self.__sa_dict[imc]['c2']
                return c1 + (c2 - c1) * np.log(pp / 0.15) / np.log(0.8 / 0.15)
            elif pp <= 5.0:
                return self.__sa_dict[imc]['c2']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c2
                return self.__sa_dict[imc]['c2']
        else:
            raise ValueError('unknown IMT %r' % imt)

    def __GM2otherSigma(self, imt, imc):
        """
        Helper function to extract coefficients from the parameters for
        converting standard deviations.

        Args:
            imt (IMT): OQ intensity measure type.
            imc (IMC): OQ intensity measure component.

        Returns:
            tuple: Coefficients: R (ratio of sigma values), standard deviation
            of sigma ratios.

        Raises:
            ValueError if unknown IMT specified.
        """

        if imc == const.IMC.AVERAGE_HORIZONTAL:
            # The amps are already in the B&B "reference" type ("GM", i.e.,
            # geometric mean)
            return 1., 0.

        self._verifyConversion(imc)

        if 'PGA' in imt:
            return self.__pga_dict[imc]['R'],\
                self.__pga_dict[imc]['c34']
        elif 'PGV' in imt:
            return self.__pgv_dict[imc]['R'],\
                self.__pgv_dict[imc]['c34']
        elif 'SA' in imt:
            R = self.__sa_dict[imc]['R']
            pp = imt.period
            if pp <= 0.15:
                return R, self.__sa_dict[imc]['c3']
            elif pp < 0.8:
                c3 = self.__sa_dict[imc]['c3']
                c4 = self.__sa_dict[imc]['c4']
                return R, c3 + (c4 - c3) * np.log(pp / 0.15) /\
                    np.log(0.8 / 0.15)
            elif pp <= 5.0:
                return R, self.__sa_dict[imc]['c4']
            else:
                # Not sure what's right here; should probably raise an error
                # but for now let's just use c4
                return R, self.__sa_dict[imc]['c4']
        else:
            raise ValueError('unknown IMT %r' % imt)

    def _verifyConversion(self, imc_in, imc_out=None):
        """
        Helper method to ensure that the conversion is possible.

        Args:
            imc_in (IMC): OpenQuake IMC type of the input amp array.
            imc_out (IMC): Desired OpenQuake IMC type of the output amps.
                Default is None. Ignored by this method.

        Raises:
            ValueError if imc_in is not valid.
        """
        if imc_in != const.IMC.GREATER_OF_TWO_HORIZONTAL and \
           imc_in != const.IMC.MEDIAN_HORIZONTAL and \
           imc_in != const.IMC.GMRotI50 and \
           imc_in != const.IMC.RotD50 and \
           imc_in != const.IMC.RANDOM_HORIZONTAL and \
           imc_in != const.IMC.HORIZONTAL:
            raise ValueError('unknown IMC %r' % imc_in)
