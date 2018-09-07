#!/usr/bin/env python

import os.path

from impactutils.io.smcontainers import ShakeMapOutputContainer
from shakemap.utils.macros import get_macros


def test_macros():
    homedir = os.path.dirname(os.path.abspath(__file__))
    shakedir = os.path.abspath(os.path.join(homedir, '..', '..', '..'))
    out_file = os.path.join(shakedir, 'tests', 'data',
                            'containers', 'northridge',
                            'shake_result.hdf')
    container = ShakeMapOutputContainer.load(out_file)
    info = container.getMetadata()
    macros = get_macros(info)
    test_dict = {'LON': '-118.5357',
                 'LAT': '34.213',
                 'DATE': 'Jan 01, 1994',
                 'DEP': 'northridge',
                 'MAG': '6.7',
                 'LOC': 'Northridge',
                 'NETID': 'ci',
                 'DATETIME': '1994-01-17T12:30:55.000000Z',
                 'EVENTID': 'northridge',
                 'TIME': '12:30:55'}
    assert macros == test_dict


if __name__ == '__main__':
    test_macros()
