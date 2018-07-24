#!/usr/bin/env python

import os.path

from shakelib.utils.containers import ShakeMapOutputContainer
from shakemap.utils.macros import get_macros

def test_macros():
    homedir = os.path.dirname(os.path.abspath(__file__))  # where is this script?
    shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
    out_file = os.path.join(shakedir, 'tests', 'data',
                            'containers', 'northridge',
                            'shake_result.hdf')
    container = ShakeMapOutputContainer.load(out_file)
    info = container.getMetadata()
    macros = get_macros(info)
    test_dict = {'LAT': '34.213',
                 'MAG': '6.7',
                 'DATETIME': '1994-01-17T12:30:55Z',
                 'DEP': 'ci3144585',
                 'LOC': '1km NNW of Reseda, CA',
                 'NETID': '',
                 'DATE': 'Jan 01, 1994',
                 'LON': '-118.537',
                 'TIME': '12:30:55',
                 'EVENTID': 'ci3144585'}
    assert macros == test_dict

if __name__ == '__main__':
    test_macros()
