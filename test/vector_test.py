#!/usr/bin/env python

#stdlib imports
import os.path
import sys

#third party
import numpy as np

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.vector import Vector

def test():
    print('Testing Vector class...')
    a = Vector(1,1,1)
    b = Vector(2,2,2)
    c = Vector(1,1,1)
    np.testing.assert_almost_equal(a.getArray(),np.array([1,1,1]))
    assert a == c
    alen = a.mag()
    np.testing.assert_almost_equal(alen,1.73205,decimal=5)
    anorm = a.norm()
    bnorm = b.norm()
    assert anorm == bnorm
    acrossb = a.cross(b)
    assert acrossb == Vector(0,0,0)
    adotb = a.dot(b)
    assert adotb == 6
    aplusb = a + b
    print('Passed Vector class tests.')
    
if __name__ == '__main__':
    test()
    
