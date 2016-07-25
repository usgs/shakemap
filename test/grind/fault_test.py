#!/usr/bin/env python

#stdlib imports
import os.path
import sys
import io

#third party
import numpy as np

import pytest

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.fault import Fault
from shakemap.utils.exception import ShakeMapException

def test_northridge():
    #this should fail!
    fault_text = """
    # Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """
    cbuf = io.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
    quad = fault.getQuadrilaterals()[0]
    topdist = quad[0].distance(quad[1])
    fault.multiplyFaultLength(2.0)    
    quad = fault.getQuadrilaterals()[0]
    topdist2 = quad[0].distance(quad[1])
    x = 1
    
def test_correct():
    #this fault should parse correctly
    fault_text = """#SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, and W. Page (2002). The Surface Rupture and Slip Distribution of the 17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. Seism. Soc. Am. 92, 43-60.
    40.70985 29.33760 0
    40.72733 29.51528 0
    40.72933 29.51528 20
    40.71185 29.33760 20
    40.70985 29.33760 0
    >
    40.70513 29.61152 0
    40.74903 29.87519 0
    40.75103 29.87519 20
    40.70713 29.61152 20
    40.70513 29.61152 0
    >
    40.72582 29.88662 0
    40.72336 30.11126 0
    40.73432 30.19265 0
    40.73632 30.19265 20
    40.72536 30.11126 20
    40.72782 29.88662 20
    40.72582 29.88662 0
    >
    40.71210 30.30494 0
    40.71081 30.46540 0
    40.70739 30.56511 0
    40.70939 30.56511 20
    40.71281 30.46540 20
    40.71410 30.30494 20
    40.71210 30.30494 0
    >
    40.71621 30.57658 0
    40.70068 30.63731 0
    40.70268 30.63731 20
    40.71821 30.57658 20
    40.71621 30.57658 0
    >
    40.69947 30.72900 0
    40.79654 30.93655 0
    40.79854 30.93655 20
    40.70147 30.72900 20
    40.69947 30.72900 0
    >
    40.80199 30.94688 0
    40.84501 31.01799 0
    40.84701 31.01799 20
    40.80399 30.94688 20
    40.80199 30.94688 0"""

    cbuf = io.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)

def test_incorrect():
    fault_text = """# Source: Ji, C., D. V. Helmberger, D. J. Wald, and K.-F. Ma (2003). Slip history and dynamic implications of the 1999 Chi-Chi, Taiwan, earthquake, J. Geophys. Res. 108, 2412, doi:10.1029/2002JB001764.
    24.27980 120.72300	0 
    24.05000 121.00000	17
    24.07190 121.09300	17
    24.33120 121.04300	17
    24.27980 120.72300	0 
    >   
    24.27980 120.72300	0
    23.70000 120.68000	0
    23.60400 120.97200	17
    24.05000 121.00000	17
    24.27980 120.72300	0
    >
    23.60400 120.97200	17 
    23.70000 120.68000	0 
    23.58850 120.58600	0
    23.40240 120.78900	17
    23.60400 120.97200	17"""

    cbuf = io.StringIO(fault_text)
    with pytest.raises(ShakeMapException):
        fault = Fault.readFaultFile(cbuf)

def test_trace():
    xp0 = [0.0]
    xp1 = [0.0]
    yp0 = [0.0]
    yp1 = [0.05]
    zp = [0.0]
    widths = [10.0]
    dips = [45.0]

    fault = Fault.fromTrace(xp0,yp0,xp1,yp1,zp,widths,dips,reference='From J Smith, (personal communication)')
    fstr = io.StringIO()
    fault.writeFaultFile(fstr)
    print(fstr.getvalue())
    
    xp0 = [-121.81529,-121.82298]
    xp1 = [-121.82298,-121.83068]
    yp0 = [37.73707,37.74233]
    yp1 = [37.74233,37.74758]
    zp = [10,15]
    widths = [15.0,20.0]
    dips = [30.0,45.0]
    fault = Fault.fromTrace(xp0,yp0,xp1,yp1,zp,widths,dips,reference='From J Smith, (personal communication)')
