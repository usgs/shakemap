#!/usr/bin/env python

# stdlib imports
import os.path
import io

# third party libraries
from configobj import ConfigObj
from validate import Validator, VdtTypeError, VdtParamError


def getCustomValidator():
    fdict = {
        'file_type': file_type,
        'directory_type': directory_type,
        'annotatedfloat_type': annotatedfloat_type,
        'gmpe_type': gmpe_type, }

    validator = Validator(fdict)
    return validator


def createDefaultConfig(configspec, includeComments=False):
    configfolder = os.path.join(os.path.expanduser('~'), '.shakemap')
    if not os.path.isdir(configfolder):
        os.mkdir(configfolder)
    outfile = os.path.join(configfolder, 'config.ini')
    config = ConfigObj(configspec=configspec, stringify=False)
    validator = getCustomValidator()
    config.validate(validator, copy=True)
    lines = config.write()
    f = open(outfile, 'wt')
    # TODO - should we write out the docs for the parameters that have
    # default=None?
    for line in lines:
        if line.strip().startswith('#') and not includeComments:
            continue
        parts = line.split('=')
        # this is a hack because I can't figure out what to do with floats with default value of None.
        # tried stringify=False/True, still get errors when I try to validate
        if len(parts) > 1 and parts[1].strip() == '""':
            continue
        if not len(line.strip()):
            continue
        f.write(line + '\n')
    f.close()
    return outfile


def whatIs(configspec, param):
    config = ConfigObj(configspec=configspec)
    validator = getCustomValidator()
    config.validate(validator, copy=True)
    comment = '%s is not a value found in %s.' % (param, configspec)
    for sectionkey in config.sections:
        section = config[sectionkey]
        for key in list(section.keys()):
            if param.lower() == key.lower():
                comment = 'Section [%s] has that option:' % (
                    sectionkey) + '\n'.join(section.comments[key])
    return comment


def annotatedfloat_type(value):
    try:
        out = float(value)
    except:
        if value.find('c') < 0 and value.find('m') < 0:
            raise VdtTypeError(value)
    out = []
    if value.find('c') > 0:
        try:
            out = float(value.replace('c', '')) / 3600.0
        except:
            raise VdtTypeError(value)
    if value.find('m') > 0:
        try:
            out = float(value.replace('m', '')) / 60.0
        except:
            raise VdtTypeError(value)
    return out


def gmpe_type(value, *args):
    if len(args) != len(value):
        raise VdtParamError('gmpe', value)
    out = []
    try:
        minmag = float(value[0])
        maxmag = float(value[1])
        mindep = float(value[2])
        maxdep = float(value[3])
        if not (minmag >= 0.0 and minmag <= 9.9):
            raise VdtTypeError(value)
        if not (maxmag >= 0.0 and maxmag <= 9.9 and maxmag > minmag):
            raise VdtTypeError(value)
        if not (mindep >= 0.0 and mindep <= 10000.0):
            raise VdtTypeError(value)
        if not (maxdep >= 0.0 and maxdep <= 10000.0 and maxdep > mindep):
            raise VdtTypeError(value)
        out += [minmag, maxmag, mindep, maxdep]
    except:
        raise VdtParamError('gmpe', value)
    return out


def file_type(value):
    if not os.path.isfile(value):
        raise VdtTypeError(value)
    return value


def directory_type(value):
    if not os.path.isdir(value):
        raise VdtTypeError(value)
    return value


def validate(configspec, configfile, macros=None):
    '''return a validated config object.
    '''
    # first, replace all the macros if we have the values
    if macros is not None:
        configstr = open(configfile, 'rt').read()
        for key, value in macros.items():
            macro = '<' + key.upper() + '>'
            configstr = configstr.replace(macro, value)
        configfile = io.StringIO(configstr)
    config = ConfigObj(configfile, configspec=configspec)
    validator = getCustomValidator()
    result = config.validate(validator)
    # for rkey,rvalue in result.iteritems():
    #     if not rvalue and rkey.find('transfer') > -1:
    #         if config[rkey]
    if result == True:
        return config
    else:
        return False
