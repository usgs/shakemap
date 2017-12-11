import numpy as np
import numexpr as ne


class GodaAtkinson2010(object):
    """
    Imlements the Goda and Atkinson (2010) spatial correlation model for IMTs.

    To do
        - Inherit from SpatialCorrelation class.

    References:
        Goda, K., & Atkinson, G. M. (2010). Intraevent spatial correlation of
        ground-motion parameters using SK-net data. Bulletin of the
        Seismological Society of America, 100(6), 3055-3067.
        `[link] <http://www.bssaonline.org/content/100/6/3055.short>`__
    """
    @staticmethod
    def getSpatialCorrelation(dists, imt):
        """
        Method for evalulating spatial correlation model.

        Args:
            dists (ndarray): Numpy array of distances (km).
            imt (IMT): Openquake intensity measure type instance.
                `[link] <http://docs.openquake.org/oq-hazardlib/master/imt.html>`__

        Returns:
            ndarray: Numpy array of correlation values.
        """  # noqa
        if 'PGA' in imt:
            alpha = 0.060
            beta = 0.283
            gamma = 5.0
        elif 'PGV' in imt:
            # Here we use the average values because there is no PGV in G&A
            alpha = 0.054
            beta = 0.319
            gamma = 5.0
        elif 'SA' in imt:
            pp = imt.period
            if pp == 0.1:
                alpha = 0.062
                beta = 0.276
                gamma = 5.0
            elif pp == 0.2:
                alpha = 0.073
                beta = 0.248
                gamma = 5.0
            elif pp == 0.3:
                alpha = 0.086
                beta = 0.219
                gamma = 5.0
            elif pp == 0.5:
                alpha = 0.073
                beta = 0.248
                gamma = 5.0
            elif pp == 1.0:
                alpha = 0.051
                beta = 0.329
                gamma = 5.0
            elif pp == 2.0:
                alpha = 0.061
                beta = 0.421
                gamma = 3.035
            elif pp == 3.0:
                alpha = 0.092
                beta = 0.671
                gamma = 1.189
            elif pp == 5.0:
                alpha = 0.071
                beta = 0.741
                gamma = 1.201
            else:
                # Here we use the average values because we don't have terms
                # for this period
                alpha = 0.054
                beta = 0.319
                gamma = 5.0
        else:
            # Again we use the average values because we don't know the IMT
            alpha = 0.054
            # Variable used within ne.evaluate, not seen by linter
            beta = 0.319  # noqa
            gamma = 5.0
            pass
        nal = -1.0 * alpha  # noqa
        gm1 = gamma - 1.0  # noqa
#        cor = 1.0 - np.sqrt(1.0 - np.maximum(
#                gamma * np.exp(nal * np.power(dists, beta)) - gm1,
#                0))
#        tmp = ne.evaluate("gamma * exp(nal * dists**beta) - gm1")
#        tmp[tmp < 0] = 0
#        cor = ne.evaluate("1.0 - sqrt(1.0 - tmp)")
        cor = ne.evaluate(
            "1.0 - sqrt(1.0 - (gamma * exp(nal * dists**beta) - gm1))")
#        cor[cor < 0] = 0
        return np.clip(cor, 0, 1)
