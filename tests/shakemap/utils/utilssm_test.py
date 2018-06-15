#!/usr/bin/env python

import os.path
import pytest
from _pytest.monkeypatch import MonkeyPatch

from configobj import ConfigObj

from shakemap.utils.utils import (get_network_name,
                                  migrate_gmpe,
                                  set_gmpe,
                                  get_object_from_config,
                                  query_yes_no)
from shakemap.utils.config import get_data_path

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))


def test_network():
    assert get_network_name('foo') == 'unknown'
    assert get_network_name(
        'us') == 'USGS National Earthquake Information Center, PDE'


def test_migrate():
    with pytest.raises(OSError):
        new_gmpe, ref = migrate_gmpe('CY08')

    data_path = get_data_path()
    cfile = os.path.join(data_path, 'migrate.conf')
    config = ConfigObj(cfile)

    with pytest.raises(KeyError):
        new_gmpe, ref = migrate_gmpe('NotAGMPE', config)

    new_gmpe, ref = migrate_gmpe('CY08', config)
    assert new_gmpe == 'CY14'

    new_gmpe, ref = migrate_gmpe('Kanno2006', config)
    assert new_gmpe == 'Kea06s'

    new_gmpe, ref = migrate_gmpe('MA2005', config)
    assert new_gmpe == 'ASK14'


def test_set_gmpe():
    config = {'modeling': {},
              'gmpe_sets': {}}
    cfg = set_gmpe('ASK14', config, 'TestEvent')
    assert cfg['modeling']['gmpe'] == 'gmpe_TestEvent_custom'
    assert cfg['gmpe_sets']['gmpe_TestEvent_custom']['gmpes'] == ['ASK14']
    assert cfg['gmpe_sets']['gmpe_TestEvent_custom']['weights'] == [1.0]


def test_get_object_from_config():
    config = {'modeling': {'gmice': 'WGRW12'},
              'gmice_modules': {'WGRW12': ['WGRW12', 'shakelib.gmice.wgrw12']}}
    gmice = get_object_from_config('gmice', 'modeling', config)
    assert gmice.getName() == 'Worden et al. (2012)'


def test_query_yes_no(monkeypatch):

    monkeypatch.setattr('builtins.input', lambda: 'y')
    torf = query_yes_no('Question?', default='yes')
    assert torf is True

    monkeypatch.setattr('builtins.input', lambda: 'n')
    torf = query_yes_no('Question?', default='yes')
    assert torf is False

    monkeypatch.setattr('builtins.input', lambda: '')
    torf = query_yes_no('Question?', default='yes')
    assert torf is True

    monkeypatch.setattr('builtins.input', lambda: 'ye')
    torf = query_yes_no('Question?', default=None)
    assert torf is True

    monkeypatch.setattr('builtins.input', lambda: 'no')
    torf = query_yes_no('Question?', default='no')
    assert torf is False

    with pytest.raises(ValueError):
        torf = query_yes_no('Question?', default='Hello')


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_network()
    test_migrate()
    test_set_gmpe()
    test_get_object_from_config()
    test_query_yes_no(MonkeyPatch())
