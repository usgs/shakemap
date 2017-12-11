#!/usr/bin/env python

from shakelib.utils.imt_string import oq_to_file, file_to_oq

def test_imt():
    """Test the imt string functions.

    """
    assert oq_to_file('SA(1.0)') == 'PSA1p0'
    assert oq_to_file('SA(0.3)') == 'PSA0p3'
    assert oq_to_file('SA(15.0)') == 'PSA15p0'
    assert oq_to_file('SA(3)') == 'PSA3p0'
    assert oq_to_file('SA(.5)') == 'PSA0p5'

    try:
        _ = oq_to_file('SA()')
    except ValueError as ve:
        assert 1==1
    
    assert file_to_oq('PSA1p0') == 'SA(1.0)'
    assert file_to_oq('PSA0p3') == 'SA(0.3)'
    assert file_to_oq('PSA15p0') == 'SA(15.0)'

    try:
        _ = file_to_oq('PSA2')
    except ValueError as ve:
        assert 1==1

    try:
        _ = file_to_oq('PSA2p')
    except ValueError as ve:
        assert 1==1
    

if __name__ == '__main__':
    test_imt()
