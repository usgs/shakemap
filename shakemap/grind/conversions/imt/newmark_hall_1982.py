import numpy as np

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

    The sigma factor was computed from an average sigma value determined
    by plotting PGV/PSA10 versus Distance for earthquakes with magnitudes
    greater than or equal to 5.0.
    
    To do
        - Inherit from ConvertIMT class. 

    References: 
        Newmark, N. M., & Hall, W. J. (1982). Earthquake spectra and design. 
        Earthquake Engineering Research Institute, El Cerrito, California. 
    """

    __vfact = 37.27 * 2.54
    __vsigma = 1.6315154

    @staticmethod
    def pgv2psa10(pgv, sigmaGlin):
        """
        Convert PGV in cm/s to PSA10 (spectral acceleration with oscillator
        period of 1.0 sec) in g.

        **Important:** PGV and sigma must be linear units. 

        :param pgv:
            Numpy array or float of PGV values; linear units (cm/s).
        :param sigmaGlin:
            Numpy array or float of standard deviation of PGV from a GMPE;
            linear units.
        :returns:
            Values converted to PSA10 and total standard deviation.
        """
        pgv = NewmarkHall1982.__vfact * NewmarkHall1982.__vfact
        sigmaG = np.log(sigmaGlin)
        sigmaTot = np.sqrt(((sigmaG) ** 2) + ((NewmarkHall1982.__vsigma) ** 2))
        return pgv, exp(sigmaTot)
    
    
    @staticmethod
    def psa102pgv(psa10, sigmaGlin):
        """
        Convert PSA10 (spectral acceleration with oscillator period of 1.0 sec)
        in g to PGV cm/s.

        **Important:** PSA10 and sigma must be linear units. 

        :param psa10:
            Numpy array or float of PSA10 values; linear units (%g).
        :param sigmaGlin:
            Numpy array or float of standard deviation of PGV from a GMPE;
            linear units.
        :returns:
            Values converted to PGV and total standard deviation.
        """
        pgv = psa10 * NewmarkHall1982.__vfact
        sigmaG = np.log(sigmaGlin)
        sigmaTot = np.sqrt(((sigmaG) ** 2) + ((NewmarkHall1982.__vsigma) ** 2))
        return pgv, exp(sigmaTot)

    
    @staticmethod
    def getVfact():
        """
        :returns: 
            The Newmark and Hall (1982) conversion factor. 
        """
        return NewmarkHall1982.__vfact
