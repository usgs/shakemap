
class BommerAlarcon2006(object):

    """
    Class for conversion between PGV (units of cm/s) and PSA05 (units of g)
    by Bommer and Alarcon (2006).

    - PSA05 stands for spectral acceleration with oscillator period of 0.5 sec
    - PGV is peak ground velocity. 

    To do
        - Inherit from ConvertIMT class. 

    References:
        Bommer, J. J., & Alarcon, J. E. (2006). The prediction and use of peak 
        ground velocity. Journal of Earthquake Engineering, 10(01), 1-31.
        `[link] <http://www.worldscientific.com/doi/abs/10.1142/S1363246906002463>`__
    """

    __vfact = 1.0 / (20.0) * 100.0 * 9.81

    @staticmethod
    def pgv2psa05(pgv):
        """
        Convert PGV in cm/s to PSA05 in g.
        **Important:** PGV must be linear units.

        :param pgv:
            Numpy array or float of PGV values; linear units.
        :returns:
            Numpy array or float of PSA05 (spectral acceleration with oscillator
            period of 0.5 sec) converted from PGV.
        """
        return pgv / BommerAlarcon2006.__vfact

    @staticmethod
    def psa052pgv(psa05):
        """
        Convert PSA05 (spectral acceleration with oscillator period of 0.5 sec)
        in g to PGV cm/s.
        **Important:** PSA10 must be linear units.

        :param psa05:
            Numpy array or float of PSA05 values; linear units.
        :returns:
            Numpy array or float of PGV converted from psa05.
        """
        return psa05 * BommerAlarcon2006.__vfact

    @staticmethod
    def getVfact():
        """
        :returns: 
            The Bommer and Alarcon (2006) conversion factor. 
        """
        return BommerAlarcon2006.__vfact
