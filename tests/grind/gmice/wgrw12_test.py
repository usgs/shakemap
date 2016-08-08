
import numpy as np

import shakemap.grind.gmice.wgrw12 as w12

def test_wgrw12():

    w12().getMIfromGM(amps, imt, dists=None, mag=None)
    
    w12().getGMfromMI(mmi, imt, dists=None, mag=None)
    
    w12().getGM2MIsd()
    
    w12().getMI2GMsd()
    
    w12().getName()
    
    w12().getScale()
    
    w12().getMinMax()
    
    w12().getDistanceType()
    
    



