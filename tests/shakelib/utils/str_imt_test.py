#!/usr/bin/env python

from shakelib.utils.imt_string import oq_to_file, file_to_oq


def test_imt():
    """Test the imt string functions.

    """
    assert oq_to_file('SA(1.0)') == 'psa1p0'
    assert oq_to_file('SA(0.3)') == 'psa0p3'
    assert oq_to_file('SA(15.0)') == 'psa15p0'
    assert oq_to_file('SA(3)') == 'psa3p0'
    assert oq_to_file('SA(.5)') == 'psa0p5'

    try:
        _ = oq_to_file('SA()')
    except ValueError as ve:
        assert 1 == 1

    assert file_to_oq('psa1p0') == 'SA(1.0)'
    assert file_to_oq('psa0p3') == 'SA(0.3)'
    assert file_to_oq('psa15p0') == 'SA(15.0)'

    try:
        _ = file_to_oq('psa2')
    except ValueError as ve:
        assert 1 == 1

    try:
        _ = file_to_oq('psa2p')
    except ValueError as ve:
        assert 1 == 1

    # Test that a fileimt that is the same as ['PGA', 'PGV', 'MMI']
    # is simply returned
    assert file_to_oq('pga') == 'PGA'
    assert file_to_oq('pgv') == 'PGV'
    assert file_to_oq('mmi') == 'MMI'


if __name__ == '__main__':
    test_imt()
