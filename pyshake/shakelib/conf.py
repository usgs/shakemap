#!/usr/bin/env python

#stdlib imports
import os.path

#third party libraries
from configobj import ConfigObj
from validate import Validator,VdtTypeError

def getCustomValidator():
    fdict = {
        'file_type': file_type,
        'annotatedfloat_type':annotatedfloat_type,}
    validator = Validator(fdict)
    return validator

def createDefaultConfig(configspec,includeComments=False):
    configfolder = os.path.join(os.path.expanduser('~'),'.shakemap')
    if not os.path.isdir(configfolder):
        os.mkdir(configfolder)
    outfile = os.path.join(configfolder,'config.ini')
    config = ConfigObj(configspec=configspec,stringify=False)
    validator = getCustomValidator()
    config.validate(validator,copy=True)
    lines = config.write()
    f = open(outfile,'wt')
    #TODO - should we write out the docs for the parameters that have default=None?
    for line in lines:
        if line.strip().startswith('#') and not includeComments:
            continue
        parts = line.split('=')
        #this is a hack because I can't figure out what to do with floats with default value of None.
        #tried stringify=False/True, still get errors when I try to validate
        if len(parts) > 1 and parts[1].strip() == '""': 
            continue
        if not len(line.strip()):
            continue
        f.write(line+'\n')
    f.close()
    return outfile

def whatIs(configspec,param):
    config = ConfigObj(configspec=configspec)
    validator = getCustomValidator()
    config.validate(validator,copy=True)
    comment = '%s is not a value found in %s.' % (param,configspec)
    for sectionkey in config.sections:
        section = config[sectionkey]
        for key in section.keys():
            if param.lower() == key.lower():
                comment = 'Section [%s] has that option:' % (sectionkey) +'\n'.join(section.comments[key])
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
            out = float(value.replace('c',''))/3600.0
        except:
            raise VdtTypeError(value)
    if value.find('m') > 0:
        try:
            out = float(value.replace('m',''))/60.0
        except:
            raise VdtTypeError(value)
    return out

def file_type(value):
    #value="file_type('/Users/mhearne/src/python/testme.py')"
    sidx = value.find('(')+1
    eidx = value.find(')')
    pat = '''['"](.*?)['"]'''
    fname = re.findall(pat,value[sidx:eidx])[0]
    if not os.path.isfile(fname):
        raise VdtTypeError(fname)
    return fname

def validate(configspec,configfile):
    config = ConfigObj(configfile,configspec=configspec)
    validator = getCustomValidator()
    result = config.validate(validator)
    if result == True:
        return True
    else:
        return False

if __name__ == '__main__':
    homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    configspec = os.path.join(homedir,'configspec.ini')
    configfile = createDefaultConfig(configspec)

    print whatIs(configspec,'pgm2mi')

    print validate(configspec,configfile)
