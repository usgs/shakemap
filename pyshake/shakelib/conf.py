#!/usr/bin/env python

#stdlib imports
import os.path

#third party libraries
from configobj import ConfigObj
from validate import Validator

def createDefaultConfig(configspec,includeComments=False):
    configfolder = os.path.join(os.path.expanduser('~'),'.shakemap')
    if not os.path.isdir(configfolder):
        os.mkdir(configfolder)
    outfile = os.path.join(configfolder,'config.ini')
    config = ConfigObj(configspec=configspec,stringify=False)
    validator = Validator()
    config.validate(validator,copy=True)
    lines = config.write()
    f = open(outfile,'wt')
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
    validator = Validator()
    config.validate(validator,copy=True)
    comment = '%s is not a value found in %s.' % (param,configspec)
    for sectionkey in config.sections:
        section = config[sectionkey]
        for key in section.keys():
            if param.lower() == key.lower():
                comment = 'Section [%s] has that option:' % (sectionkey) +'\n'.join(section.comments[key])
    return comment

def validate(configspec,configfile):
    config = ConfigObj(configfile,configspec=configspec)
    validator = Validator()
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
