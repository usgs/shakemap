#!/usr/bin/env python

import os
import os.path
import pytest
import copy

from configobj import ConfigObj
from validate import ValidateError

import shakemap.utils.config as config

homedir = os.path.dirname(os.path.abspath(__file__))
shakedir = os.path.abspath(os.path.join(homedir, '..', '..'))


def test_config():

    install_dir, data_dir = config.get_config_paths()
    #
    # get_logger()
    #
    log_file = os.path.join(data_dir, 'nc72282711', 'shake.log')
    if os.path.isfile(log_file):
        os.remove(log_file)
    logger = config.get_logger('nc72282711', log_file=True,
                               log_option='debug')
    logger.debug('xyxyxyzz')
    with open(log_file, 'r') as log_fd:
        line = log_fd.readline()
        assert 'xyxyxyzz' in line
    os.remove(log_file)

    logger = config.get_logger('nc72282711', log_option='quiet')
    logger = config.get_logger('nc72282711')
    logger = config.get_logger('nc72282711', log_option='debug')

    #
    # Some stuff we just call and see if it bombs out
    #
    mydatapath = config.get_data_path()
    myinstall, mydata = config.get_config_paths()

    myspec = config.get_configspec()
    myvalid = config.get_custom_validator()

    c1 = ConfigObj(os.path.join(mydatapath, "model.conf"),
                   configspec=myspec)
    c2 = ConfigObj(os.path.join(mydatapath, "modules.conf"),
                   configspec=myspec)
    c3 = ConfigObj(os.path.join(mydatapath, "gmpe_sets.conf"),
                   configspec=myspec)
    c4 = ConfigObj(os.path.join(mydatapath, "northridge_model.conf"),
                   configspec=myspec)
    c5 = ConfigObj(os.path.join(mydatapath, "products.conf"),
                   configspec=myspec)
    c1.merge(c2)
    c1.merge(c3)
    c1.merge(c4)
    c1.merge(c5)

    results = c1.validate(myvalid, preserve_errors=True)

    assert isinstance(results, bool) and results

    config.check_config(c1, logger)
    #
    # Break the config
    #
    ctest = copy.deepcopy(c1)
    ctest['modeling']['ccf'] = 'NotACCF'
    with pytest.raises(ValidateError):
        config.check_config(ctest, logger)
    ctest = copy.deepcopy(c1)
    ctest['modeling']['ipe'] = 'NotAnIPE'
    with pytest.raises(ValidateError):
        config.check_config(ctest, logger)
    ctest = copy.deepcopy(c1)
    ctest['modeling']['gmice'] = 'NotAGMICE'
    with pytest.raises(ValidateError):
        config.check_config(ctest, logger)
    ctest = copy.deepcopy(c1)
    ctest['modeling']['gmpe'] = 'NotAGMPE'
    with pytest.raises(ValidateError):
        config.check_config(ctest, logger)

    ctest = copy.deepcopy(c1)
    ctest['modeling']['gmpe'] = 47
    results = ctest.validate(myvalid, preserve_errors=True)
    assert isinstance(results, dict)
    with pytest.raises(RuntimeError):
        config.config_error(ctest, results)

    ctest = copy.deepcopy(c1)
    del ctest['interp']
    results = ctest.validate(myvalid, preserve_errors=True)
    assert isinstance(results, dict)
    with pytest.raises(RuntimeError):
        config.config_error(ctest, results)

    #
    # Test the profile checker
    #
    ctest = ConfigObj()
    ctest['profiles'] = {'prof1': {'data_path': '/xyz/zzsx/zz',
                                   'install_path': '/xyz/zzsx/zz'},
                         'prof2': {'data_path': data_dir,
                                   'install_path': install_dir}}
    ct1 = config.check_profile_config(ctest)
    assert 'prof1' not in list(ct1['profiles'].keys())
    # os.remove(config_file)
    #
    # annotatedfloat_type()
    #
    res = config.annotatedfloat_type('4.0')
    assert isinstance(res, float)
    assert res == 4.0
    res = config.annotatedfloat_type('4.0d')
    assert isinstance(res, float)
    assert res == 4.0
    res = config.annotatedfloat_type('4.0m')
    assert isinstance(res, float)
    assert res == 4.0 / 60.0
    res = config.annotatedfloat_type('4.0c')
    assert isinstance(res, float)
    assert res == 4.0 / 3600.0
    with pytest.raises(ValidateError):
        res = config.annotatedfloat_type('4.0caweoifaw')
    with pytest.raises(ValidateError):
        res = config.annotatedfloat_type('')
    #
    # weight_list()
    #
    res = config.weight_list(['0.2', '0.3', '0.5'], min=0)
    assert isinstance(res, list)
    assert res == [0.2, 0.3, 0.5]
    res = config.weight_list('None', min=0)
    assert isinstance(res, list)
    assert res == []
    res = config.weight_list('[]', min=0)
    assert isinstance(res, list)
    assert res == []
    res = config.weight_list(['0.2', '0.3', '0.5'], min=3)
    assert isinstance(res, list)
    assert res == [0.2, 0.3, 0.5]
    with pytest.raises(ValidateError):
        res = config.weight_list([], min=1)
    with pytest.raises(ValidateError):
        res = config.weight_list('[]', min=1)
    with pytest.raises(ValidateError):
        res = config.weight_list('[None]', min=1)
    with pytest.raises(ValidateError):
        res = config.weight_list('None', min=1)
    with pytest.raises(ValidateError):
        res = config.weight_list(['0.2', '0.3', '0.5'], min=4)
    with pytest.raises(ValidateError):
        res = config.weight_list(['-0.2', '0.3', '0.5'], min=3)
    with pytest.raises(ValidateError):
        res = config.weight_list(['0.1', '0.3', '0.5'], min=3)
    #
    # gmpe_list()
    #
    res = config.gmpe_list('[]', min=0)
    assert isinstance(res, list)
    assert res == []
    with pytest.raises(ValidateError):
        res = config.gmpe_list('[]', min=1)
    res = config.gmpe_list('thing1', min=0)
    assert isinstance(res, list)
    assert res == ['thing1']
    res = config.gmpe_list(['thing1'], min=0)
    assert isinstance(res, list)
    assert res == ['thing1']
    res = config.gmpe_list(['thing1', 'thing2'], min=0)
    assert isinstance(res, list)
    assert res == ['thing1', 'thing2']
    with pytest.raises(ValidateError):
        res = config.gmpe_list(['thing1', 'thing2'], min=3)
    with pytest.raises(ValidateError):
        res = config.gmpe_list(7, min=0)
    with pytest.raises(ValidateError):
        res = config.gmpe_list([7], min=0)
    #
    # extent_list()
    #
    res = config.extent_list('[]')
    assert isinstance(res, list)
    assert res == []
    res = config.extent_list([])
    assert isinstance(res, list)
    assert res == []
    with pytest.raises(ValidateError):
        res = config.extent_list(7)
    with pytest.raises(ValidateError):
        res = config.extent_list(['-20.0', '-10.0', '20.0'])
    with pytest.raises(ValidateError):
        res = config.extent_list(['-20.0', '-10.0', '20.0', 'thing'])
    with pytest.raises(ValidateError):
        res = config.extent_list(['-20.0', '-10.0', '20.0', '1000.0'])
    res = config.extent_list(['-20.0', '-10.0', '20.0', '10.0'])
    assert isinstance(res, list)
    assert res == [-20.0, -10.0, 20.0, 10.0]
    #
    # file_type()
    #
    res = config.file_type('None')
    assert isinstance(res, str)
    assert not res
    with pytest.raises(ValidateError):
        res = config.file_type('/home/xxxyyyzzz/awefawe')
    res = config.file_type(os.path.abspath(__file__))
    assert isinstance(res, str)
    assert res == os.path.abspath(__file__)
    #
    # directory_type()
    #
    res = config.directory_type('None')
    assert isinstance(res, str)
    assert not res
    with pytest.raises(ValidateError):
        res = config.directory_type('/home/xxxyyyzzz/awefawe')
    res = config.directory_type(os.path.dirname(os.path.abspath(__file__)))
    assert isinstance(res, str)
    assert res == os.path.dirname(os.path.abspath(__file__))
    #
    # status_string()
    #
    res = config.status_string('', min=1)
    assert res == 'automatic'
    res = config.status_string('automatic', min=1)
    assert res == 'automatic'
    with pytest.raises(ValidateError):
        res = config.status_string('thing', min=1)
    #
    # cfg_float_list()
    #
    res = config.cfg_float_list(['2.0', '3.0', '4.0'])
    assert res == [2.0, 3.0, 4.0]
    res = config.cfg_float_list('2.0')
    assert res == [2.0]
    with pytest.raises(ValidateError):
        res = config.cfg_float_list('')
    with pytest.raises(ValidateError):
        res = config.cfg_float_list({})
    with pytest.raises(ValidateError):
        res = config.cfg_float_list({'a': 'b'})
    with pytest.raises(ValidateError):
        res = config.cfg_float_list([])
    with pytest.raises(ValidateError):
        res = config.cfg_float_list('thing')
    #
    # cfg_float()
    #
    res = config.cfg_float('2.0')
    assert res == 2.0
    with pytest.raises(ValidateError):
        res = config.cfg_float(['2.0'])
    with pytest.raises(ValidateError):
        res = config.cfg_float('')
    with pytest.raises(ValidateError):
        res = config.cfg_float('None')
    with pytest.raises(ValidateError):
        res = config.cfg_float('thing')


if __name__ == '__main__':
    os.environ['CALLED_FROM_PYTEST'] = 'True'
    test_config()
