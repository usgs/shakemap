#!/usr/bin/env python

#stdlib imports

#third party imports
import numpy as np

#local imports
from gmpe import GMPE
from distance import getCartesianDistance
import BB


DEFAULT_RREF = 1
DEFAULT_MREF = 4.5

def repmat(a, m, n):
    """
    Mimic the behavior of Matlab's repmat() function.
    @param a: 2D numpy array.
    @param m: Desired number of output rows.
    @param n: Desired number of output columns.
    @return: 2D numpy array consisting of input array replicated M rows by N columns.
    """
    if a.ndim == 1:
        a = np.array([a])
    (origrows, origcols) = a.shape
    rows = origrows * m
    cols = origcols * n
    b = a.reshape(1,a.size).repeat(m, 0).reshape(rows, origcols).repeat(n, 0)
    return b.reshape(rows, cols)

class BA08(GMPE):
    def __init__(self,sourceparams,fault,sitecorrect,mechanism=None,doPSA=True,verbose=True):
        #We can use the superclass constructor here - this is the weird syntax to do that, apparently
        #super(BA08,self,sourceparams,fault,sitecorrect,doPSA).__init__()
        GMPE.__init__(self,sourceparams,fault,sitecorrect,doPSA)
        self.__setup() #establish all of our constant values
        if mechanism is None:
            self.mechanism = 'ALL'
            print 'BA08: source mechanism "%s".' % self.mechanism
        if mechanism not in self.MECHANISM.keys():
            if verbose:
                print 'BA08: Unknown source mechanism "%s". Using default.' % mechanism
            self.mechanism = 'ALL'
        self.fault = fault
        self.__calculateDistances()
        self.biasApplied = False #we optionally apply the bias in the super class after initialization
        
    def __calculateDistances(self):
        if self.fault is None:
            #get the latitude and longitude of every cell in the grid
            ny = self.nrows
            nx = self.ncols
            dx = self.xres
            dy = self.yres
            ulx = self.xmin
            uly = self.ymax
            rowlat = uly - np.arange(0,ny).reshape(-1,1) * dy
            collon = ulx + np.arange(0,nx) * dx
            lat = repmat(rowlat,1,nx)
            lon = repmat(collon,ny,1)
            self.distance = getCartesianDistance(lat,lon,0,self.lat,self.lon,0)/1000.0 #cartesian distance, km
        else:
            raise Exception,"Distance to fault not yet supported!"

    def __setup(self):
        """
        Set up all of our tables of constants for use later in calculations.
        """
        self.bias = {'pgv':1.0,
                     'pga':1.0,
                     'psa03':1.0,
                     'psa10':1.0,
                     'psa30':1.0}
        
        dsckeys = ['c1','c2','c3','h']
        dscvalues = [[-0.87370, 0.10060, -0.00334, 2.54],#pgv
                     [-0.66050, 0.11970, -0.01151, 1.35],#pga
                     [-0.55430, 0.01955, -0.00750, 2.14],#psa03
                     [-0.81830, 0.10270, -0.00334, 2.54],#psa10
                     [-0.78440, 0.07282, -0.00191, 2.83]]#psa30
        self.DSC = {'pgv':dict(zip(dsckeys,dscvalues[0][:])),
                    'pga':dict(zip(dsckeys,dscvalues[1][:])),
                    'psa03':dict(zip(dsckeys,dscvalues[2][:])),
                    'psa10':dict(zip(dsckeys,dscvalues[3][:])),
                    'psa30':dict(zip(dsckeys,dscvalues[4][:]))}
                     
        msckeys = ['e1', 'e2', 'e3', 'e4', 'e5', 'e6', 'e7', 'Mh']
        mscvalues = [[ 5.00121,  5.04727,  4.63188,  5.08210, 0.18322, -0.12736, 0.00000, 8.50], #pgv
                     [-0.53804, -0.50350, -0.75472, -0.50970, 0.28805, -0.10164, 0.00000, 6.75], #pga
                     [ 0.43825,  0.44516,  0.25356,  0.51990, 0.64472, -0.15694, 0.10601, 6.75], #psa03
                     [-0.46896, -0.43443, -0.78465, -0.39330, 0.67880, -0.18257, 0.05393, 6.75], #psa10
                     [-1.82979, -1.74690, -2.22584, -1.91814, 0.77966, -0.45384, 0.67466, 6.75]] #psa30
        self.MSC = {'pgv':dict(zip(msckeys,mscvalues[0][:])),
                    'pga':dict(zip(msckeys,mscvalues[1][:])),
                    'psa03':dict(zip(msckeys,mscvalues[2][:])),
                    'psa10':dict(zip(msckeys,mscvalues[3][:])),
                    'psa30':dict(zip(msckeys,mscvalues[4][:]))}

        aukeys = ['sigma','tauU','sigmaTU','tauM','sigmaTM']
        auvalues = [[0.500, 0.286, 0.576, 0.256, 0.560], #pgv
                    [0.502, 0.265, 0.566, 0.260, 0.564], #pga
                    [0.546, 0.272, 0.608, 0.269, 0.608], #psa03
                    [0.573, 0.318, 0.654, 0.302, 0.647], #psa10
                    [0.566, 0.410, 0.700, 0.401, 0.695]] #psa30
        self.AU = {'pgv':dict(zip(aukeys,auvalues[0][:])),
                   'pga':dict(zip(aukeys,auvalues[1][:])),
                   'psa03':dict(zip(aukeys,auvalues[2][:])),
                   'psa10':dict(zip(aukeys,auvalues[3][:])),
                   'psa30':dict(zip(aukeys,auvalues[4][:]))}

        #Period-dependent site amplification coefficients
        pdskeys = ['blin','b1','b2']
        pdsvalues = [[-0.600, -0.500, -0.06],
                     [-0.360, -0.640, -0.14],
                     [-0.440, -0.520, -0.14],
                     [-0.700, -0.440,  0.00],
                     [-0.740, -0.340,  0.00]]
        self.PDSAC = {'pgv':dict(zip(pdskeys,pdsvalues[0][:])),
                      'pga':dict(zip(pdskeys,pdsvalues[1][:])),
                      'psa03':dict(zip(pdskeys,pdsvalues[2][:])),
                      'psa10':dict(zip(pdskeys,pdsvalues[3][:])),
                      'psa30':dict(zip(pdskeys,pdsvalues[4][:]))}

        # Period-independent site amplification coefficients
        self.PISAC = {'a1':0.03,
                      'pga_low':0.06,
                      'a2':0.09,
                      'V1':180.0,
                      'V2':300.0,
                      'Vref':760.0}
        mekeys = ['U','S','N','R']
        mevalues = [[ 1,0,0,0 ], #ALL - unspecified
                    [ 0,1,0,0 ], #SS - strike-slip
                    [ 0,0,1,0 ], #NM - normal
                    [ 0,0,0,1 ]] #RS - Reverse (thrust fault)
        self.MECHANISM = {'ALL':dict(zip(mekeys,mevalues[0][:])),
                          'SS':dict(zip(mekeys,mevalues[1][:])),
                          'NM':dict(zip(mekeys,mevalues[2][:])),
                          'RS':dict(zip(mekeys,mevalues[3][:]))}
        
    def calculateAmplitudes(self,Vs30,eps=0.0):
        """
        Implement the ground motion prediction equation as described in Boore and Atkinson 2008.
        @param Vs30: 2D numpy array containing Vs30 values (must be same shape as internal data grids)
        @keyword eps: Epsilon value to use if you need to produce NON-median ground motions.
        """
        mech = self.MECHANISM[self.mechanism]
        M = self.mag
        for gm in ['pgv','pga','psa03','psa10','psa30']:
            bias = self.bias[gm]
            dsc = self.DSC[gm]
            msc = self.MSC[gm]
            au = self.AU[gm]
            pds = self.PDSAC[gm]
            
            #Get the distance scaling coefficients
            c1,c2,c3,h = dsc['c1'],dsc['c2'],dsc['c3'],dsc['h']
            
            Rref = DEFAULT_RREF
            Mref = DEFAULT_MREF
            R = np.sqrt(np.power(self.distance,2) + np.power(h,2)) #Eq 4
            Fd = (c1 + c2*(self.mag-Mref))*np.log(R/Rref) + c3*(R-Rref) #Eq 3

            #Get the magnitude scaling coefficients
            e1,e2,e3,e4,e5,e6,e7,Mh = msc['e1'],msc['e2'],msc['e3'],msc['e4'],msc['e5'],msc['e6'],msc['e7'],msc['Mh']

            #Get the mechanism "dummy" variables
            U,SS,NS,RS = mech['U'],mech['S'],mech['N'],mech['R']

            if M <= Mh:
                Fm = e1*U + e2*SS + e3*NS + e4*RS + e5*(M-Mh) + np.power(e6*(M-Mh),2) #Eq 5a
            else:
                Fm = e1*U + e2*SS + e3*NS + e4*RS + e7*(M-Mh) #Eq 5b

            #site amplification
            #Get the period dependent site amplification coefficients
            blin,b1,b2 = pds['blin'],pds['b1'],pds['b2']

            #Get the period independent site amplification coefficients
            a1,a2,pga_low,V1,V2,Vref = self.PISAC['a1'],self.PISAC['a2'],self.PISAC['pga_low'],self.PISAC['V1'],self.PISAC['V2'],self.PISAC['Vref']
            
            Flin = blin * np.log(Vs30/Vref) #Eq 7
            pga4nl = np.exp(Fm + Fd) #Eq 1, where Fs = 0 and eps = 0
            dx = np.log(a2/a1) #Eq 11
            bnl = np.zeros(self.pga.shape)
            bnl[Vs30 > V1] = b1 #Eq 13a
            vidx1 = ((Vs30 > V1) & (Vs30 <= V2))
            vidx2 = ((Vs30 > V2) & (Vs30 < Vref))
            bnl[vidx1] = (b1-b2)*np.log(Vs30[vidx1]/V2)/(np.log(V1/V2) + b2) #Eq 13b
            bnl[vidx2] = b2*np.log(Vs30[vidx2]/Vref)/np.log(V2/Vref) #Eq 13c
            bnl[Vref <= Vs30] = 0.0 #Eq 13d
            dy = bnl*np.log(a2/pga_low) #Eq 12
            c = (3*dx - bnl*dx)/np.power(dx,2) #Eq 9
            d = (-1*(2*dy - bnl*dx))/np.power(dx,3) #Eq 10
            

            #vectorized assignment here
            Fnl = np.zeros(self.pga.shape)
            pidx1 = (pga4nl <= a1) #all the indices in pga4nl <= a1
            pidx2 = ((pga4nl > a1) & (pga4nl <= a2))
            pidx3 = (pga4nl > a2)
            Fnl[pidx1] = bnl[pidx1]*np.log(pga_low/0.1) #Eq 8a
            #Eq 8b
            Fnl[pidx2] = bnl[pidx2]*np.log(pga_low/0.1) + c[pidx2]*np.power(np.log(pga4nl[pidx2]/a1),2) + d[pidx2]*np.power(np.log(pga4nl[pidx2]/a1),3)
            Fnl[pidx3] = bnl[pidx3] * np.log(pga4nl[pidx3]/0.1) #Eq 8c

            Fs = Flin + Fnl
            
            sigmaT = self.getUncertaintyTerm(au,gm)
            if np.isnan(sigmaT):
                sigmaT = 0.0

            #All ground motions are calculated the same way, so we can dynamically generate the line of code
            #and use eval().  Hopefully not too confusing.
            gm_mult = BB.getRand2Max_pgm(gm) #get the multiplier to convert random component GMs to peak GMs.
            if gm == 'pga':
                self.pga = np.exp(Fm + Fd + Fs + eps*sigmaT)*gm_mult
            elif gm == 'pgv':
                self.pgv = np.exp(Fm + Fd + Fs + eps*sigmaT)*gm_mult
            elif gm == 'psa03':
                self.psa03 = np.exp(Fm + Fd + Fs + eps*sigmaT)*gm_mult
            elif gm == 'psa10':
                self.psa10 = np.exp(Fm + Fd + Fs + eps*sigmaT)*gm_mult
            else: #psa30
                self.psa30 = np.exp(Fm + Fd + Fs + eps*sigmaT)*gm_mult

            
    def getUncertaintyTerm(self,au,gm):
        """
        Determine the uncertainty terms in the GMPE equation (sigma and tau)
        @param au: Dictionary containing keys: ['sigma','tauU','sigmaTU','tauM','sigmaTM']
        @return: SigmaT as defined in equation 1
        """
        if self.fault is not None:
            if self.biasApplied:
                sigma = au['sigma']
            else:
                sigma = au['sigmaTM']
            tau = au['tauM']
        else:
            if self.biasApplied:
                sigma = au['sigma']
            else:
                sigma = au['sigmaTU']
            tau = au['tauU']
        sigmaT = np.sqrt(np.power(sigma,2) + np.power(tau,2))
        #modify as per BB06...
        sigmaT = BB.getRand2Max_sigma(gm,sigmaT)
        return sigmaT


if __name__ == '__main__':
    sparams = {'lat':32.1234,'lon':-119.1234,'depth':24,'mag':7.5}
    fault = None
    scorrect = None
    myba = BA08(sparams,fault,scorrect)
    vs30 = np.ones(myba.pga.shape)*625.0
    myba.calculateAmplitudes(vs30)
    print 'PGA Range:'
    print myba.pga.min(),myba.pga.max()
    print 'PGV Range:'
    print myba.pgv.min(),myba.pgv.max()
    print 'PSA03 Range:'
    print myba.psa03.min(),myba.psa03.max()
    print 'PSA10 Range:'
    print myba.psa10.min(),myba.psa10.max()
    print 'PSA30 Range:'
    print myba.psa30.min(),myba.psa30.max()
        
            
                                                           
            
