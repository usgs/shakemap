#!/usr/bin/env python

import os.path
import tempfile
import textwrap
import sys

from configobj import ConfigObj

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))

from validate import ValidateError

import shakelib.utils.config as config

def test_config():

    #
    # Should test get_data_path() and get_config_path() here
    # but I don't know how that stuff will work on Travis,
    # and we haven't run shakeprofile on Travis... Probably
    # ought to look into how to do that, and maybe add arguments
    # to the functions to use a specific path when '~' isn't 
    # meaningful.
    #

    #
    # Some stuff we just call and see if it bombs out
    #
    junk = config.get_custom_validator()

    #
    # annotatedfloat_type()
    #
    try:
        res = config.annotatedfloat_type('4.0')
        assert isinstance(res, float)
        assert res == 4.0
    except ValidateError:
        assert False
    try:
        res = config.annotatedfloat_type('4.0d')
        assert isinstance(res, float)
        assert res == 4.0
    except ValidateError:
        assert False
    try:
        res = config.annotatedfloat_type('4.0m')
        assert isinstance(res, float)
        assert res == 4.0 / 60.0
    except ValidateError:
        assert False
    try:
        res = config.annotatedfloat_type('4.0c')
        assert isinstance(res, float)
        assert res == 4.0 / 3600.0
    except ValidateError:
        assert False
    try:
        res = config.annotatedfloat_type('4.0caweoifaw')
        assert False
    except ValidateError:
        assert True
    try:
        res = config.annotatedfloat_type('')
        assert False
    except ValidateError:
        assert True
    #
    # weight_list()
    #
    try:
        res = config.weight_list(['0.2', '0.3', '0.5'], min=0)
        assert isinstance(res, list)
        assert res == [0.2, 0.3, 0.5]
    except ValidateError:
        assert False
    try:
        res = config.weight_list('None', min=0)
        assert isinstance(res, list)
        assert res == []
    except ValidateError:
        assert False
    try:
        res = config.weight_list('[]', min=0)
        assert isinstance(res, list)
        assert res == []
    except ValidateError:
        assert False
    try:
        res = config.weight_list('[]', min=1)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.weight_list('[None]', min=1)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.weight_list(['0.2', '0.3', '0.5'], min=4)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.weight_list(['0.2', '0.3', '0.5'], min=3)
        assert isinstance(res, list)
        assert res == [0.2, 0.3, 0.5]
    except ValidateError:
        assert False
    try:
        res = config.weight_list(['-0.2', '0.3', '0.5'], min=3)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.weight_list(['0.1', '0.3', '0.5'], min=3)
        assert False
    except ValidateError:
        assert True
    #
    # gmpe_list()
    #
    try:
        res = config.gmpe_list('[]', min=0)
        assert isinstance(res, list)
        assert res == []
    except ValidateError:
        assert False
    try:
        res = config.gmpe_list('[]', min=1)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.gmpe_list('thing1', min=0)
        assert isinstance(res, list)
        assert res == ['thing1']
    except ValidateError:
        assert False
    try:
        res = config.gmpe_list(['thing1'], min=0)
        assert isinstance(res, list)
        assert res == ['thing1']
    except ValidateError:
        assert False
    try:
        res = config.gmpe_list(['thing1', 'thing2'], min=0)
        assert isinstance(res, list)
        assert res == ['thing1', 'thing2']
    except ValidateError:
        assert False
    try:
        res = config.gmpe_list(['thing1', 'thing2'], min=3)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.gmpe_list(7, min=0)
        assert False
    except ValidateError:
        assert True
    #
    # extent_list()
    #
    try:
        res = config.extent_list('[]')
        assert isinstance(res, list)
        assert res == []
    except ValidateError:
        assert False
    try:
        res = config.extent_list(7)
        assert False
    except ValidateError:
        assert True
    try:
        res = config.extent_list(['-20.0', '-10.0', '20.0'])
        assert False
    except ValidateError:
        assert True
    try:
        res = config.extent_list(['-20.0', '-10.0', '20.0', 'thing'])
        assert False
    except ValidateError:
        assert True
    try:
        res = config.extent_list(['-20.0', '-10.0', '20.0', '1000.0'])
        assert False
    except ValidateError:
        assert True
    try:
        res = config.extent_list(['-20.0', '-10.0', '20.0', '10.0'])
        assert isinstance(res, list)
        assert res == [-20.0, -10.0, 20.0, 10.0]
    except ValidateError:
        assert False
    #
    # file_type()
    #
    try:
        res = config.file_type('None')
        assert isinstance(res, str)
        assert not res
    except ValidateError:
        assert False
    try:
        res = config.file_type('/home/xxxyyyzzz/awefawe')
        assert False
    except ValidateError:
        assert True
    try:
        res = config.file_type(os.path.abspath(__file__))
        assert isinstance(res, str)
        assert res == os.path.abspath(__file__)
    except ValidateError:
        assert False
    #
    # directory_type()
    #
    try:
        res = config.directory_type('None')
        assert isinstance(res, str)
        assert not res
    except ValidateError:
        assert False
    try:
        res = config.directory_type('/home/xxxyyyzzz/awefawe')
        assert False
    except ValidateError:
        assert True
    try:
        res = config.directory_type(os.path.dirname(os.path.abspath(__file__)))
        assert isinstance(res, str)
        assert res == os.path.dirname(os.path.abspath(__file__))
    except ValidateError:
        assert False

if __name__ == '__main__':
    test_config()
