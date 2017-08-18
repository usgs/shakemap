#!/usr/bin/env python
import os.path

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

from shakemap.mapping.sm_colordict import SM_colordict

def test_colordict():

    smcd = SM_colordict()

    cmap = smcd.getSMcolordict()

    assert cmap['blue'][5] == [0.5, 0.5764705882352941, 0.5764705882352941]
    assert cmap['green'][5] == [0.5, 1.0, 1.0] 
    assert cmap['red'][5] == [0.5, 0.47843137254901963, 0.47843137254901963]

if __name__ == '__main__':
    test_colordict()
