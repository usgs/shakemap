"""
Implements the GMICE of Worden et al., 2013, BSSA, 102(1), pp. 204-221
"""

import numpy as np

class WGRW12(object):
      #--------------------------------------------------------------------------
      #
      # MMI = c2->C1 + c2->C2 * log(Y)  for log(Y) <= c2->T1
      # MMI = C1 + C2 * log(Y)          for c2->T1 < log(Y) <= T1
      # MMI = C3 + C4 * log(Y)          for log(Y) > T1
      #
      # or
      #
      # MMI = c2->C1 + c2->C2 * log(Y) + C5 + C6 * log(D) + C7 * M   
      #                            for log(Y) <= c2->T1
      # MMI = C1 + C2 * log(Y) + C5 + C6 * log(D) + C7 * M for c2->T1 < log(Y) <= T1
      # MMI = C3 + C4 * log(Y) + C5 + C6 * log(D) + C7 * M for log(Y) > T1
      #
      # Limit the distance residuals to between 10 and 300 km.
      # Limit the magnitude residuals to between M3.0 and M7.3.
      #
      #--------------------------------------------------------------------------

    __constants = { \
        'pga'   : { 'C1' :  1.78, 'C2' :  1.55, 'C3' : -1.60, 'C4' : 3.70, 'C5'   : -0.91, \
                    'C6' :  1.02, 'C7' : -0.17, 'T1' :  1.57, 'T2' : 4.22, 'SMMI' : 0.66,  \
                    'SPGM' : 0.35 }, \
        'pgv'   : { 'C1' :  3.78, 'C2' :  1.47, 'C3' :  2.89, 'C4' : 3.16, 'C5'   :  0.90, \
                    'C6' :  0.00, 'C7' : -0.18, 'T1' :  0.53, 'T2' : 4.56, 'SMMI' : 0.63,  \
                    'SPGM' : 0.38 }, \
        'psa03' : { 'C1' :  1.26, 'C2' :  1.69, 'C3' : -4.15, 'C4' : 4.14, 'C5'   : -1.05, \
                    'C6' :  0.60, 'C7' :  0.00, 'T1' :  2.21, 'T2' : 4.99, 'SMMI' : 0.82,  \
                    'SPGM' : 0.44 }, \
        'psa10' : { 'C1' :  2.50, 'C2' :  1.51, 'C3' :  0.20, 'C4' : 2.90, 'C5'   :  2.27, \
                    'C6' : -0.49, 'C7' : -0.29, 'T1' :  1.65, 'T2' : 4.98, 'SMMI' : 0.75,  \
                    'SPGM' : 0.47 }, \
        'psa30' : { 'C1' :  3.81, 'C2' :  1.17, 'C3' :  1.99, 'C4' : 3.01, 'C5'   :  1.91, \
                    'C6' : -0.57, 'C7' : -0.21, 'T1' :  0.99, 'T2' : 4.96, 'SMMI' : 0.89,  \
                    'SPGM' : 0.64 } \
    }

    __constants2 = { \
        'pga'   : { 'C1' :  1.71, 'C2' :  2.08, 'T1' :  0.14, 'T2' : 2.0 }, \
        'pgv'   : { 'C1' :  4.62, 'C2' :  2.17, 'T1' : -1.21, 'T2' : 2.0 }, \
        'psa03' : { 'C1' :  1.15, 'C2' :  1.92, 'T1' :  0.44, 'T2' : 2.0 }, \
        'psa10' : { 'C1' :  2.71, 'C2' :  2.17, 'T1' : -0.33, 'T2' : 2.0 }, \
        'psa30' : { 'C1' :  7.35, 'C2' :  3.45, 'T1' : -1.55, 'T2' : 2.0 }  \
    }

    def getMIfromGM(self, amps, imt, dists=None, mag=None):
        """ 
        Function getMIfromGM
        Put some documentation here.
        """
        c, c2 = self.__getConsts(imt)

        if dists is not None and mag is not None:
            doresid = True
            ldd = np.log10(np.clip(dists, 10, 300))
            if mag < 3.0:
                mag = 3.0
            elif mag > 7.3:
                mag = 7.3
        else:
            doresid = False

        #
        # Convert (for accelerations) from %g to cm/s^2
        # then take the log10
        #
        if 'PGV' not in imt:
            units = 9.81
        else:
            units = 1.0
        lamps = np.log10(amps * units)
        mmi = np.zeros_like(amps)
        #
        # This is the MMI 1 to 2 range that is discussed in the paper but not 
        # specifically implemented
        #
        idx = lamps < c2['T1']
        mmi[idx] = c2['C1'] + c2['C2'] * lamps[idx]
        #
        # This is the lower segment of the bi-linear fit
        #
        idx = (lamps >= c2['T1']) & (lamps < c['T1'])
        mmi[idx] = c['C1'] + c['C2'] * lamps[idx]
        #
        # This is the upper segment of the bi-linear fit
        #
        idx = lamps >= c['T1']
        mmi[idx] = c['C3'] + c['C4'] * lamps[idx]

        if doresid:
            mmi += c['C5'] + c['C6'] * ldd[idx] + c['C7'] * mag

        mmi = np.clip(mmi, 1.0, 10.0)
        return mmi

    def getGMfromMI(self, mmi, imt, dists=None, mag=None):
        """ 
        Function getGMfromMI
        Put some documentation here.
        """
        c, c2 = self.__getConsts(imt)

        if dists is not None and mag is not None:
            doresid = True
            ldd = np.log10(np.clip(dists, 10, 300))
            if mag < 3.0:
                mag = 3.0
            elif mag > 7.3:
                mag = 7.3
        else:
            doresid = False

        if doresid:
            mmi -= c['C5'] + c['C6'] * ldd[idx] + c['C7'] * mag

        pgm = np.zeros_like(mmi)

        #
        # MMI 1 to 2
        #
        idx = mmi < 2.0
        pgm[idx] = np.power(10, (mmi[idx] - c2['C1']) / c2['C2'])
        #
        # Lower segment of bi-linear relationship
        #
        idx = (mmi >= 2.0) & (mmi < c['T2'])
        pgm[idx] = np.power(10, (mmi[idx] - c['C1']) / c['C2'])
        #
        # Upper segment of bi-linear relationship
        #
        idx = mmi >= c['T2']
        pgm[idx] = np.power(10, (mmi[idx] - c['C3']) / c['C4'])
        if 'PGV' not in imt:
            units = 9.81
        else:
            units = 1.0
        pgm /= units

        return pgm
        
    def getGM2MIsd(self):
        """ Return dictionary of GM to MI sigmas (in MMI units)"""
        return { 'pga'   : self.__constants['pga']['SMMI'],
                 'pgv'   : self.__constants['pgv']['SMMI'],
                 'psa03' : self.__constants['psa03']['SMMI'],
                 'psa10' : self.__constants['psa10']['SMMI'],
                 'psa30' : self.__constants['psa30']['SMMI'] }
        
    def getMI2GMsd(self):
        """ Return dictionary of MI to GM sigmas (in linear PGM units)"""
        return { 'pga'   : 10**self.__constants['pga']['SPGM'],
                 'pgv'   : 10**self.__constants['pgv']['SPGM'],
                 'psa03' : 10**self.__constants['psa03']['SPGM'],
                 'psa10' : 10**self.__constants['psa10']['SPGM'],
                 'psa30' : 10**self.__constants['psa30']['SPGM'] }

    def getName(self):
        return 'Worden et al. (2012)'

    def getScale(self):
        return 'scale_wgrw12.ps'

    def getMinMax(self):
        return (1.0, 10.0)

    def getDistanceType(self):
        return 'rrup'
          
    def __getConsts(self, imt):
        """ Helper function to get the constants """

        if 'PGA' in imt:
            c = self.__constants['pga']
            c2 = self.__constants2['pga']
        elif 'PGV' in imt:
            c = self.__constants['pgv']
            c2 = self.__constants2['pgv']
        elif 'SA' in imt:
            pp = imt.period
            if pp == 0.3:
                c = self.__constants['psa03']
                c2 = self.__constants2['psa03']
            elif pp == 1.0:
                c = self.__constants['psa10']
                c2 = self.__constants2['psa10']
            elif pp == 3.0:
                c = self.__constants['psa30']
                c2 = self.__constants2['psa30']
            else:
                raise ValueError("Unknown SA period: %f" % pp)
        else:
            raise ValueError("Unknown IMT %r" % imt)
        return (c, c2)
