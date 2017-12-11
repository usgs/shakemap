import numpy as np


class Boore2003(object):
    """
    Implements the correlation model proposed in the appendix of Boore et al.
    (2003).

    To do
        - Inherit from SpatialCorrelation class.

    References:
        Boore, D. M., Gibbs, J. F., Joyner, W. B., Tinsley, J. C., & Ponti, D.
        J. (2003). Estimated ground motion from the 1994 Northridge,
        California, earthquake at the site of the Interstate 10 and La Cienega
        Boulevard bridge collapse, West Los Angeles, California. Bulletin of
        the Seismological Society of America, 93(6), 2737-2751.
        `[link] <http://www.bssaonline.org/content/93/6/2737.short>`__
    """
    @staticmethod
    def getSpatialCorrelation(dists, imt='PGA'):
        """
        Method for evalulating spatial correlation model.

        Args:
            dists (ndarray): Numpy array of distances (km).
            imt (string): Openquake intensity measure type instance string.
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__
                This model was developed specifically for PGA, and so this is
                the default value.

        Returns:
            Numpy array of correlation values.
        """  # noqa
        if imt != 'PGA':
            raise Exception('PGA is the only supported IMT.')
        return np.exp(-1.0 * np.sqrt(0.6 * dists))
