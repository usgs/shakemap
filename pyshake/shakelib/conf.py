#!/usr/bin/env python

#stdlib imports
import os.path
import textwrap
import tempfile

#third party
import numpy as np

CLASSES = ['GrindConfig','MappingConfig']
TYPES = {'int':int,'float':float,'string':(str,unicode),'file':(str,unicode),
         'list':list,'dict':dict,'bool':bool}

class AnnotatedValue(object):
    def __init__(self,valuestr):
        #split the valuestr into a number piece and a string piece
        #string pieces can be 'c', 'm', 'k' (second, minute, kilometer)
        if valuestr.find('c') > -1:
            self.value = float(valuestr.replace('c','')) * 1/3600.0
        elif valuestr.find('m') > -1:
            self.value = float(valuestr.replace('m','')) * 1/60.0
        else:
            self.value = float(valuestr) #value is in decimal degrees
            
    def getValue(self):
        return self.value

class ConfigException(Exception):
    """
    Class to represent errors in the Fault class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def getClass(classname):
    myclass = None
    if classname == 'GrindConfig':
        myclass = GrindConfig()
    elif classname == 'MappingConfig':
        myclass = MappingConfig()
    return myclass
    
def getConfig(configfile,section=None):
    validate(configfile)
    variables = {}
    myglobals = {}
    execfile(configfile,myglobals,variables)
    if section is None:
        return variables
    if section not in variables.keys():
        raise LookupError('Section "%s" not found in config file.' % section)
    return variables[section]
    
def validate(configfile):
    missingMessage = 'Missing variables:'
    wrongTypeMessage = 'Incorrect types:'
    outsideMessage = 'Out of range:'

    #repeat for all subclasses of Config (possible to loop over these?)
    hasError = False
    for classname in CLASSES:
        config = getClass(classname)
        tmissing,twrongtypes,toutsiderange = config.validate(configfile)
        if len(tmissing) or len(twrongtypes) or len(toutsiderange):
            hasError = True
        for missing in tmissing:
            missingMessage += '\t%s : %s\n' % (classname,missing)
        for wrongtype in twrongtypes:
            wrongTypeMessage += '\t%s : %s\n' % (classname,wrongtype)
        for outside in toutsiderange:
            outsideMessage += '\t%s : %s\n' % (classname,outside)

    if hasError:
        raise ConfigException(missingMessage+wrongTypeMessage+outsideMessage)
    
class Config(object):
    def validate(self,configfile):
        tkeys = dir(self)
        keys = [] #the list of attribute dictionaries we explicitly created (not builtins)
        for key in tkeys:
            if (key.startswith('__') and key.endswith('__')) or not key.endswith('_var'):
                continue
            #the dictionaries have _var on the end of them, we're looking for bare names
            keys.append(key.replace('_var','')) 

        variables = {} #this will have a dictionary of classes in it
        myglobals = {}
        execfile(configfile,myglobals,variables)
        classname = type(self).__name__
        classobject = variables[classname]
        extras = []
        missing = []
        missingdoc = []
        wrongtypes = []
        outsiderange = []
        notsupported = []
        for variable in keys: #loop over the things we know are supposed to be in the config file
            varinfo = getattr(self,variable+'_var')
            isRequired = varinfo.has_key('required') and varinfo['required']

            if variable not in classobject.__dict__.keys() and isRequired:
                missing.append(variable)
                continue

            #is the type of the variable in the config correct?
            vtype = varinfo['type']
            varvalue = classobject.__dict__[variable]
            if vtype in ['int','string','list','dict']:
                if not isinstance(varvalue,TYPES[vtype]):
                    wrongtypes.append(variable)
                    continue
            elif vtype == 'float':
                try:
                    float(varvalue)
                except ValueError:
                    wrongtypes.append(variable)
                    continue
            elif vtype == 'file':
                if not os.path.isfile(varvalue):
                    wrongtypes.append(variable)
                    continue
            elif vtype == 'enum':
                if varvalue not in varinfo['supported']:
                    wrongtypes.append(variable)
                    continue
            elif vtype == 'annotatedval':
                if not isinstance(varvalue,TYPES['string']):
                    wrongtypes.append(variable)
                    continue
                try:
                    aval = AnnotatedValue(varvalue)
                except:
                    wrongtypes.append(variable)
                    continue

            #check range of value
            if varinfo.has_key('range'):
                if len(varinfo['range']) == 1:
                    vmin = varinfo['range'][0]
                    vmax = np.inf
                else:
                    vmin = varinfo['range'][0]
                    vmax = varinfo['range'][1]
                if vtype == 'annotatedval':
                    varvalue = AnnotatedValue(varvalue).getValue()
                if varvalue < vmin or varvalue > vmax:
                    outsiderange.append(variable)

        return (missing,wrongtypes,outsiderange)

class GrindConfig(Config):
    ampfactor_file_var = {'value':None,
                          'type':'file',
                          'required':False,
                          'doc':'''use Borcherdt-style site correction tables
    File containing the Ma and Mv factors for short- and mid-period amps
    as a function of input acceleration in g. the path is relative to 
    $shake_home; lib/sitecorr/Borcherdt94.dat is the default. Old-style
    Borcherdt tables will not work with the new site amplification functions.
    Do not specify an old-style table below unless you run grind with -oldsc.
    
    If the user calls grind with -oldsc, the old-style Borcherdt table MUST
    be specified below, e.g.:
    
    ampfactor_file : lib/sitecorr/site_corr_cdmg.dat

    The structure of the old-style table is very different from the new
    table. See site_corr_cdmg.dat as an example, and src/lib/SiteCorrGrd.pm
    for a more detailed explanation of the structure.'''}

    smVs30default_var = {'value':None,
                         'type':'float',
                         'required':True,
                         'default': 686.0,
                         'range':[0.0,1500.0],
                         'fmt':'%.1f',
                         'doc':'''This parameter sets the base site velocity (Vs30) for which 
    the GMPEs will attempt to produce amps. If you are using Borcherdt-style 
    amplifications, it should be set to the velocity of the "rock" site
    class (i.e., the one that generates amplifications of unity (1.0)).
    If you use the GMPEs' native site corrections (i.e., you run grind with
    the -gmpesc flag), this value doesn't really matter much -- it should 
    just be set to something sane, or left alone. The default is 686 (m/s).'''}

    use_gmpe_sc_var = {'value':None,
                       'type':'bool',
                       'required':True,
                       'default':False,
                       'doc':'''Has the same effect as calling grind with -nativesc or 
    -gmpesc from the command line. Which is to say that grind will apply
    the site amplification factors as defined by the GMPE (see 'gmpe', below)
    rather than the Borcherdt-style corrections. Calling grind with -gmpesc
    or -nativesc will force the use of GMPE-native site corrections regardless
    of the value set here. There are two acceptable arguments: 'true' and 
    'false' ('false' is the default).
    
    Example:
    
    use_gmpe_sc : true'''}

    bad_station_var = {'value':None,
                       'type':'dict',
                       'required':False,
                       'doc':'''specifies stations to flag as bad under certain circumstances
    the format of the statement is:
    
    bad_station : code mag start_date-[end_date [mag start_date-[...]]]

    Where 'code' is the station code, 'mag' is the event magnitude cutoff 
    below which the station is considered bad, 'start_date' is the event 
    date to begin applying the cutoff, and 'end_date' is the event date 
    at which the cutoff no longer applies; dates are given in the yyyymmdd
    format, and are UTC (i.e. GMT) dates; a missing end date implys dates 
    inclusive of the current date; multiple 'mag start_date-end_date' groups 
    are allowed, e.g.:
    
    bad_station : BC3 3.8 19990101-19990407 2.2 19990407-
    
    In the above example, the station 'BC3' will be flagged as bad for
    events smaller than 3.8 from January 1, 1999 to April 7, 1999 and
    for events smaller than 2.2 from April 7, 1999 to the present.  The
    station will not be flagged for events before January 1, 1999'''}

    bias_norm_var = {'value': None,
                     'type':'enum',
                     'required':True,
                     'default':'l1',
                     'supported':['l1','l2'],
                     'doc':'''acceptable values are 'l2' (for least squares) or 'l1' (for
    absolute deviation); the default is 'l1'.'''}

    basin_module_var = {'value':None,
                        'type':'string',
                        'required':False,
                        'doc':'''specifies the module to use for performing basin depth
    corrections when grind is called with -basement.  This module should
    reside in the Basin subdirectory of the library modules.  The default
    is Field2000.
    Example:
           basin_module : Field2000'''}

    x_grid_interval_var = {'value':None,
                           'type':'annotatedval',
                           'required':True,
                           'range':[0.0,10.0],
                           'doc':'''<floating point value>[<units>]
    Where <units> is one of:
    ' ' => decimal degrees (no unit given)
    m  => arc minutes
    c  => arc seconds
    
    e.g.:
    x_grid_interval : 30c
    y_grid_interval : .5m

    sets the x and y grid sizes to 30 arc seconds.

    x_grid_interval and y_grid_interval specify the output grid spacing, and
    input grids are resampled as needed.'''}

class MappingConfig(Config):
    organization_var = {'value':None,
                    'required':True,
                    'type':'string',
                    'doc':'''Specify the name of the map-producing organization for inclusion in
    map titles; don't get carried away here, there isn't much room; statement
    format is:
    
    organization : Name of Organization

    Example:

    organization : USGS'''}

    mi_xres = {'value':None,
               'required':True,
               'default':'30c',
               'type':'annotatedval',
               'doc':'''<floating point value><units>

    Where <units> is one of:
    ' ' => degrees (i.e. no units specified)
    'm' => minutes
    'c' => seconds
    
    e.g.:
    mi_xres   : 0.5m
    mi_yres   : 0.5m
    mi_xhires : 0.05m
    mi_yhires : 0.05m
    
    The default value is '30c' for low-resolution and 3c for high-resolution.'''}

def _test_validate():
    config = '''
    from pyshake.shakelib.conf import GrindConfig,MappingConfig
    GrindConfig.ampfactor_file = '/Users/mhearne/src/python/testme.py'
    GrindConfig.ampfactor_file_doc = 'This is a docstring'

    GrindConfig.smVs30default = 686
    GrindConfig.smVs30default_doc = 'This is a docstring'

    GrindConfig.use_gmpe_sc = True
    GrindConfig.use_gmpe_sc_doc = 'This is a docstring'

    GrindConfig.bad_station = {'8016':(9.9,'19990101-')}
    GrindConfig.bad_station_doc = 'This is a docstring'

    GrindConfig.bias_norm = 'l2'
    GrindConfig.bias_norm_doc = 'This is a docstring'

    GrindConfig.basin_module = 'Field2000'
    GrindConfig.basin_module_doc = 'This is a docstring'

    GrindConfig.x_grid_interval = '30c'
    GrindConfig.basin_module_doc = 'This is a docstring'

    MappingConfig.organization = 'USGS'
    MappingConfig.mi_xres = '10m'
    '''
    configfile = None
    try:
        th,configfile = tempfile.mkstemp()
        os.close(th)
        f = open(configfile,'wt')
        f.write(textwrap.dedent(config))
        f.close()
        #gc = GrindConfig()
        #missing,wrongtypes,outsiderange = gc.validate(configfile)
        validate(configfile)
        grindconfig = getConfig(configfile,'GrindConfig')
        mapconfig = getConfig(configfile,'MappingConfig')
        print grindconfig
        print mapconfig
    except Exception,eobj:
        print 'Error! - %s' % str(eobj)
        if configfile is not None and os.path.isfile(configfile):
            os.remove(configfile)

if __name__ == '__main__':
    _test_validate()

    
