
class NewmarkHall1982(object):

    """
    Class for conversion between PGA and PSA10 by Newmark and Hall (1982). 

    - PSA10 stands for spectral acceleration with oscillator period of 1.0 sec
    - PGV is peak ground velocity. 

    The conversion factor is 37.27*2.54. 
    Note that 2.54 is the conversion factor to convert from cm/s to in/s and 
    37.27 is derived from SA(f) = 1.65*(2*pi*V*f)/386.09, where:

    - SA(f) is spectral acceleration at frequency f (in g)
    - f is the frequency of interest
    - V is the velocity
    - 1.65 is the N&H amplification factor for velocity at 5% damping
    - 386.09 is the acceleration of gravity in inches per second per g

    Todo: 
        Inherit from ConvertIMT class. 
    
    References: 
        Newmark, N. M., & Hall, W. J. (1982). Earthquake spectra and design. 
        Earthquake Engineering Research Institute, El Cerrito, California. 
    """

    __vfact = 37.27 * 2.54


    @staticmethod
    def pgv2psa10(pgv):
        """
        Convert PGV in cm/s to PSA10 (spectral acceleration with oscillator
        period of 1.0 sec) in g.

        **Important:** PGV must be linear units. 

        :param pgv:
            Numpy array or float of PGV values; linear units.
        :returns:
            Values converted to PSA10.
        """
        return pgv / NewmarkHall1982.__vfact

    @staticmethod
    def psa102pgv(psa10):
        """
        Convert PSA10 (spectral acceleration with oscillator period of 1.0 sec)
        in g to PGV cm/s.

        **Important:** PSA10 must be linear units. 

        :param psa10:
            Numpy array or float of PSA10 values; linear units.
        :returns:
            Values converted to PGV.
        """
        return psa10 * NewmarkHall1982.__vfact

    @staticmethod
    def getVfact():
        """
        :returns: 
            The Newmark and Hall (1982) conversion factor. 
        """
        return NewmarkHall1982.__vfact

