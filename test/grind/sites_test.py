#!/usr/bin/env python

#stdlib imports
import sys
import os.path

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

#local imports
from shakemap.grind.sites import Sites

#local imports
from shakemap.utils.exception import ShakeMapException
        
def test_sites():
    vs30file = None
    print('Testing creation of Sites object with Vs30 file of %s...' % vs30file)
    cx = -118.2
    cy = 34.1
    dx = 0.0083
    dy = 0.0083
    xspan = 3.0
    yspan = 3.0
    mysite = Sites.createFromCenter(cx,cy,xspan,yspan,dx,dy,vs30File=vs30file,
                                    padding=True,resample=False)
    sc = mysite.getSitesContext()
    
    cx = -118.2
    cy = 83
    dx = 0.0083
    dy = 0.0083
    xspan = 3.0
    yspan = 3.0
    mysite = Sites.createFromCenter(cx,cy,xspan,yspan,dx,dy,vs30File=vs30file,padding=True,resample=False)

    xmin = 116.234
    xmax = 120.876
    ymin = 20.12345
    ymax = 24.75435
    dx = 0.0083
    dy = 0.0083
    mysite = Sites.createFromBounds(xmin,xmax,ymin,ymax,dx,dy,vs30File=vs30file,padding=False,resample=False)
    print('Passed creation of Sites object with Vs30 file of %s' % vs30file)

if __name__ == '__main__':
    test_sites()

