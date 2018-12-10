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
                 'DATE': 'Jan 17, 1994',
                 'DEP': '18.0',
                 'MAG': '6.7',
                 'LOC': 'Northridge',
                 'NETID': 'ci',
                 'DATETIME': '1994-01-17T12:30:55.000000Z',
                 'EVENTID': 'northridge',
                 'VERSION': '1',
                 'PRODUCT_CODE': 'northridge',
                 'TIME': '12:30:55'}
    assert list(sorted(macros.keys())) == list(sorted(test_dict.keys()))
    for key, mvalue in macros.items():
        tvalue = test_dict[key]
        print('Testing key %s: %s vs %s.' % (key, mvalue, tvalue))
        assert mvalue == tvalue


if __name__ == '__main__':
    test_macros()
