#!/usr/bin/env python

import os.path
import json
import pprint
import sys
import argparse
import cmd
from collections import OrderedDict

TYPES = ['bool','int','float','string','list','dict','file']
NUMTYPES = ['int','float']
MAXITER = 3
DEBUG = True

# TODO:
#  - Figure out what to do about list and dict values
#    - Don't allow them
#    - allow user to edit them as strings, save them using eval()
#    - pass to friend functions from other modules (grind, mapping, etc.) => PREFERRED
#    - allow configuration of simple dicts (one deep) and lists
#  - Decide what to do about booleans (default range from False to True, make int?)
#  - What to do about deprecated settings?
#    - Add deprecated flag in translator dictionaries?
#    - Just leave them out of program dictionary => PREFERRED
#  - What to do about NEW fields that aren't in the old file
#    - pass them to configure() methods, get them added into new data file => PREFERRED
#  - Is the caption enough information to get users to understand what the option is for?
#    - Bruce thinks not, have a separate doc string with examples and the like
#  - What to do about files that have macros in them like <SHAKE_HOME>? 
#    - Treat them as strings?
#  - List and dict values translated by separate modules could be passed all the config lines
#    associated with the values, then can handle anything that requires maintaining state.
#  - Separate modules would all have methods defined translate() and configure().
#  - Should those modules contain the data structures defined below in this code, or should
#    those data structures be saved in data files somewhere on the system?
#    Argument for putting them in code: Code is in source control, config files should not 
#    be, or at least not in the same source control system as the code.
#  - Should we track changes to config files?  How?
#  - Should have an option to generate config file from scratch using just defaults - argument 
#    for every parameter having a default value (or prompting user for those that don't).
#  - For some of these config options, the order in which they appear matters (grind->gmpe is an example).
#    - Use OrderDicts to solve this problem

# MAPPING = {'organization':{'ctype':'string','default':None,'caption':'Name of organization'},
#            'measure_unit':{'ctype':'string','default':'inch','caption':'Postscript units'},
#            'mi_xres':{'ctype':'string','default':'30c','caption':'Mapping X resolution'},
#            'mi_yres':{'ctype':'string','default':'30c','caption':'Mapping Y resolution'},
#            'network_color':{'ctype':'dict','default':None,'caption':'Network plotting color and width'},
#            'station_size':{'ctype':'float','default':0.12,'caption':'Station plotting size'},
#            'contour_color':{'ctype':'string','default':'white','caption':'Contour color'},
#            'ncontours':{'ctype':'int','default':8,'caption':'Number of contours to render'},
#            'contour_width':{'ctype':'int','default':5,'caption':'Contour width (points)'},
#            'highway_color':{'ctype':'string','default':'grey','caption':'Highway color'},
#            'highway_width':{'ctype':'int','default':6,'caption':'Highway width (points)'},
#            'fault_color':{'ctype':'string','default':'darkred','caption':'Fault color'},
#            'border_width':{'ctype':'int','default':5,'caption':'Border width (points)'},
#            'water_color':{'ctype':'string','default':'lightblue','caption':'Water color'},
#            'epicenter_size':{'ctype':'float','default':0.30,'caption':'Size of epicenter (points)'},
#            'main_title_size':{'ctype':'int','default':15,'caption':'Title size (points)'},
#            'sub_title_size':{'ctype':'int','default':12,'caption':'Subtitle size (points)'},
#            'city_dot_size':{'ctype':'float','default':0.08,'caption':'City dot size (points)'},
#            'city_name_size':{'ctype':'int','default':12,'caption':'City name size (points)'},
#            'tv_highway_width':{'ctype':'int','default':8,'caption':'Highway width for TV map (points)'},
#            'tv_highway_color':{'ctype':'string','default':'lightgrey','caption':'Highway color for TV map (points)'},
#            'tv_border_width':{'ctype':'int','default':8,'caption':'Border width for TV map (points)'},
#            'ff_width':{'ctype':'int','default':12,'caption':'Finite fault width (points)'},
#            'ff_color':{'ctype':'string','default':'darkgrey','caption':'Finite fault color'},
#            'xorig': {'ctype':'float','default':1.0,'caption':'Position of the lower left corner of the maps'},
#            'yorig': {'ctype':'float','default':1.0,'caption':'Position of the lower left corner of the maps'},
#            'width': {'ctype':'float','default':6.5,'caption':'Width of the maps'},
#            'tvwidth': {'ctype':'float','default':10.0,'caption':'Width of the maps'},
#            'map_axes':{'ctype':'string','default':'a60mf30m/a30mf30mWSen','caption':'GMT map axes definitions'},
#            'map_data_dir':{'ctype':'string','default': '<SHAKE_HOME>/lib/mapping','caption':'Directory for output map data'},
#            'topo_cmap': {'ctype':'string','default':'<MAP_DIR>/tan.cpt','caption':'Colormap file for topo data'},
#            'ii_cmap': {'ctype':'string','default':'<MAP_DIR>/tan.cpt','caption':'Colormap file for intensity data'},
#            'ii_tvmap_cmap': {'ctype':'string','default':'<MAP_DIR>/tan.cpt','caption':'Colormap file for intensity data (TV)'},
#            'sd_cmap': {'ctype':'string','default':'<MAP_DIR>/tan.cpt','caption':'Colormap file for intensity data (TV)'},
#            sd_cmap         : <MAP_DIR>/sd.cpt
# map_roads	    : <MAP_DIR>/ca_roads.xy
# map_faults          : <MAP_DIR>/global_faults.txt
# map_cities_label    : <MAP_DIR>/us_cities.txt
# map_cities_on_pgm   : false
# big_cities		: dummy
# big_cities_label	: <MAP_DIR>/ca_bigcities_label.txt
# very_big_cities		: dummy
# very_big_cities_label	: <MAP_DIR>/ca_verybigcities_label.txt
# use_gmt_coast : true
# faults_on_mi		: true
# stations_on_mi		: true
# filled_stations_on_mi 	: false
# noscenariosplash: false

           

GRIND = {'smVs30default':{'ctype':'float','default':686,'crange':[0,1500],'caption':'Vs30 default value'},
         'bad_station':{'ctype':'dict','caption':'Stations known to have bad data'},
         'gmpe':{'ctype':'dict','caption':'Ground Motion Prediction Equations (magmin,magmax,depthmin,depthmax)'},
         'outlier_deviation_level' : {'ctype':'int','default':3,'crange':[1,5],'caption':'Outlier deviation level'},
         'outlier_max_mag': {'ctype':'float','default':7.0,'crange':[0,9.9],'caption':'Maximum magnitude value for outliers'},
         'bias_norm': {'ctype':'string','default':'l1','caption':'Bias norm'},
         'bias_max_range': {'ctype':'int','default':120,'crange':[0,1000],'caption':'Maximum bias range'},
         'bias_min_stations':{'ctype':'int','default':6,'crange':[0,100],'caption':'Minimum number of stations required to compute bias'},
         'bias_max_mag': {'ctype':'float','default':7.0,'crange':[0,9.9],'caption':'Maximum magnitude value for bias'},
         'bias_max_bias': {'ctype':'float','default':2.0,'crange':[-10.0,10.0],'caption':'Maximum allowed bias'},
         'bias_min_bias': {'ctype':'float','default':-2.0,'crange':[-10.0,10.0],'caption':'Minimum allowed bias'},
         'bias_log_amp': {'ctype':'bool','default':'false','caption':'?????'},
         'direct_patch_size': {'ctype':'int','default':1000,'crange':[0,10000],'caption':'Patch size for ??'},
         'fwdata_file': {'ctype':'string','caption':'Fault width data file'},
         'topobin' : {'ctype':'string','caption':'Program used to calculate topography'}}

def getValue(ctype,tvalue):
    value = None
    if ctype == 'string':
        value = tvalue
    elif ctype == 'int':
        try:
            value = int(tvalue)
        except:
            pass
    elif ctype == 'float':
        try:
            value = float(tvalue)
        except:
            pass
    elif ctype == 'list':
        value = tvalue
    elif ctype == 'dict':
        parts = tvalue.split()
        key = parts[0]
        dlist = parts[1:]
        value = (key,dlist)
    else:
        value = tvalue
    return value

def translate(configfile,jsonfile):
    programs = {'grind':GRIND}
    program = None
    pkey = None
    for key in programs.keys():
        if configfile.find(key) > -1:
            program = programs[key]
            pkey = key
            break
    if program is None:
        raise Exception("Program not found in list")
    f = open(configfile,'rt')
    #dicts is a dictionary of dict types
    dicts = {}
    #lists is a dictionary of list types
    lists = {}
    #configs is a dictionary of simple types (strings, ints, floats, files, bools)
    configs = {}
    #loop over lines.  left side of : is the key, right side is the value
    #if key is known to be a dictionary, then assume the right side consists of a dictionary entry with key 
    #of the first part of the list, and the value is the remainder of the list.  Save the dictionary entry in 
    #global dictionary of dictionaries, then go through those at the end.
    #if key is known to be a list, then the right side should be appended to the list for that key
    for line in f.readlines():
        line = line.strip()
        if line.startswith('#'):
            continue
        if line == '':
            continue
        configkey,tvalue = line.split(':')
        configkey = configkey.strip()
        tvalue = tvalue.strip()
        if configkey not in program.keys():
            print 'Key %s not found in program %s' % (configkey,pkey)
            continue
        
        ctype = program[configkey]['ctype']
        if ctype == 'dict':
            value = getValue(ctype,tvalue)
            if value is None:
                print 'Could not parse value %s from line %s in %s' % (tvalue,line,configfile)
                continue
            if dicts.has_key(configkey):
                dicts[configkey][value[0]] = value[1]
            else:
                dicts[configkey] = OrderedDict([(value[0],value[1])])
        elif ctype == 'list':
            if lists.has_key(configkey):
                lists[configkey].append(value)
            else:
                lists[configkey] = [value]
        else:
            value = getValue(ctype,tvalue)
            if value is None:
                print 'Could not parse value %s from line %s in %s' % (tvalue,line,configfile)
                continue
            cdict = program[configkey].copy()
            cdict['value'] = value
            configs[configkey] = cdict
        
            
    f.close()
    for configkey,value in dicts.iteritems():
        cdict = program[configkey]
        cdict['value'] = value
        configs[configkey] = cdict
    for configkey,value in lists.iteritems():
        cdict = program[configkey]
        cdict['value'] = value
        configs[configkey] = cdict

    configtop = {pkey:configs}
    config = Config(configtop)
    config.save(jsonfile)
        
                           

def validateType(ctype,value):
    if ctype == 'int' and not isinstance(value,int):
        raise TypeError('Value of "%s" must be an integer' % (str(value)))
    if ctype == 'bool' and not isinstance(value,str) and value.lower() not in ['true','false']:
        raise TypeError('Value of "%s" must be an integer' % (str(value)))
    if ctype == 'float' and (not isinstance(value,(int,float))):
        raise TypeError('Value of "%s" must be a number' % (str(value)))
    if ctype == 'string' and (not isinstance(value,(str,unicode))):
        raise TypeError('Value of "%s" must be unicode or ASCII string.' % (str(value)))
    if ctype == 'list' and not isinstance(value,list):
        raise TypeError('Value of "%s" must be a list.' % (str(value)))
    if ctype == 'dict' and not isinstance(value,dict):
        raise TypeError('Value of "%s" must be a dictionary.' % (str(value)))
    if ctype == 'file':
        if not isinstance(value,str) and not isinstance(value,unicode):
            raise TypeError('Value of "%s" is not a file.' % (str(value)))
        isfile = os.path.isfile(value)
        if not isfile:
            raise TypeError('Value of "%s" is not a file.' % (str(value)))

def validateRange(crange,value):
    if value < crange[0] or value > crange[1]:
        raise ValueError('Input value %f is outside accepted range of %f to %f' % (value,crange[0],crange[1]))
                         
class ConfigItem(object):
    def __init__(self,name,caption,value,ctype='string',crange=None,default=None):
        self.name = name
        self.caption = caption
        if ctype not in TYPES:
            raise LookupError('Input type must be one of: %s' % str(TYPES))
        self.ctype = ctype
        if crange is not None and ctype not in NUMTYPES:
            raise TypeError('Input range only applies to config items of types: %s' % str(NUMTYPES))
        if crange is not None and (not isinstance(crange,(list,tuple)) or len(crange) != 2):
            raise TypeError('Range must be a list or tuple with 2 elements')
        self.crange = crange
        try:
            validateType(ctype,value)
        except TypeError,exc:
            raise exc
        if crange is not None:
            try:
                validateRange(self.crange,value)
            except ValueError,exc:
                raise exc
        self.value = value
        if default is not None:
            try:
                validateType(ctype,default)
            except TypeError,exc:
                raise exc
            if crange is not None:
                try:
                    validateRange(self.crange,default)
                except ValueError,exc:
                    raise exc
            self.default = default
        else:
            self.default = None
        return

    def getValue(self):
        return self.value
    
    def getName(self):
        return self.name
    
    def getCaption(self):
        return self.caption
    
    def getDefault(self):
        return self.default

    def getRange(self):
        return self.crange
    
    def getType(self):
        return self.ctype

    def printItem(self,ntab=0):
        print '  '*ntab + 'Config Option: %s' % self.name
        print '  '*(ntab+1) + 'Description: %s' % self.caption
        print '  '*(ntab+1) + 'Value: %s' % str(self.value)
        print '  '*(ntab+1) + 'Type: %s' % self.ctype
        print '  '*(ntab+1) + 'Default: %s' % self.default
        print '  '*(ntab+1) + 'Range: %s' % self.crange

    def configure(self):
        if self.ctype == 'string':
            fmt = 'Enter a new value for %s.  Current/Default: %s/[%s]: '
            prompt = fmt % (self.name,self.value,self.default)
            resfmt = 'Setting value of %s to %s'
        elif self.ctype == 'file':
            fmt = 'Enter a new value for %s.  Current/Default: %s/[%s]: '
            prompt = fmt % (self.name,self.value,self.default)
            resfmt = 'Setting value of %s to %s'
        elif self.ctype == 'bool':
            fmt = 'Enter a new value for %s.  Current/Default: %s/[%s]: '
            prompt = fmt % (self.name,self.value,self.default)
            resfmt = 'Setting value of %s to %s'
        elif self.ctype == 'int':
            if self.crange is not None:
                fmt = 'Enter a new value for %s between %i and %i.  Current/Default: %i/[%i]: '
                prompt = fmt % (self.name,self.crange[0],self.crange[1],self.value,self.default)
            else:
                fmt = 'Enter a new value for %s.  Current/Default: %i/[%i]: '
                prompt = fmt % (self.name,self.value,self.default)
            resfmt = 'Setting value of %s to %i'
        elif self.ctype == 'float':
            if self.crange is not None:
                fmt = 'Enter a new value for %s between %f and %f.  Current/Default: %f/[%f]: '
                prompt = fmt % (self.name,self.crange[0],self.crange[1],self.value,self.default)
            else:
                fmt = 'Enter a new value for %s.  Current/Default: %f/[%f]: '
                prompt = fmt % (self.name,self.value,self.default)
            resfmt = 'Setting value of %s to %f'
        else:
            raise TypeError('Editing of dicts and lists not currently supported')
        niter = 0
        while niter < MAXITER:
            niter += 1
            answer = raw_input(prompt)
            if answer.strip() == '':
                answer = self.default
            if self.ctype == 'int':
                try:
                    answer = int(answer)
                except:
                    print 'Expecting an integer value: Try again.'
                    continue
            if self.ctype == 'float':
                try:
                    answer = float(answer)
                except:
                    print 'Expecting a floating point value: Try again.'
                    continue
            try:
                validateType(self.ctype,answer)
            except TypeError,obj:
                print 'Expecting a value of type %s.  Try again.' % self.ctype
            print resfmt % (self.name,answer)
            self.value = answer
            break
        return

class ConfigSection(object):
    def __init__(self,configdict):
        if not isinstance(configdict,dict):
            raise TypeError('Input to contructor is not a dictionary object')
        self.configdict = {}
        for name,item in configdict.iteritems():
            if not isinstance(item,dict):
                raise TypeError('Item %s in input list is not a dictionary' % name)
            if not item.has_key('caption') or not item.has_key('value'):
                raise TypeError('Item %s in input list is missing a required field ("caption" or "value")' % name)
            caption = item['caption']
            value = item['value']
            ctype = 'string'
            if item.has_key('ctype'):
                ctype = item['ctype']
            crange = None
            if item.has_key('crange'):
                crange = item['crange']
            default = None
            if item.has_key('default'):
                default = item['default']
            configitem = ConfigItem(name,caption,value,ctype=ctype,crange=crange,default=default)
            self.configdict[name] = configitem


    def printSection(self,parameter=None,ntab=0):
        if parameter is not None:
            if parameter not in self.configdict.keys():
                raise KeyError('Parameter "%s" not found in config' % parameter)
            self.configdict[parameter].printItem(ntab=ntab+1)
            print
            return

        for name in self.configdict.keys():
            self.configdict[name].printItem(ntab=ntab+1)
            print

    def getParamNames(self):
        return self.configdict.keys()

    def getConfigItem(self,param):
        if not self.configdict.has_key(param):
            raise KeyError('ConfigSection does not have a parameter called %s' % param)
        return self.configdict[param]
            
    def configure(self,parameter=None):
        if parameter is not None:
            if parameter not in self.configdict.keys():
                raise KeyError('Parameter "%s" not found in config' % parameter)
            self.configdict[parameter].configure()
            

class Config(object):
    def __init__(self,configdict=None):
        if configdict is not None:
            self.configdict = {}
            for section,configsection in configdict.iteritems():
                self.configdict[section] = ConfigSection(configsection)

    def getSectionNames(self):
        return self.configdict.keys()

    def getSection(self,section):
        if not self.configdict.has_key(section):
            raise KeyError('Config does not have a section called %s' % param)
        return self.configdict[section]
                
    def printConfig(self,section=None,parameter=None):
        if section is not None:
            if section not in self.configdict.keys():
                raise KeyError('Section "%s" not found in config file.' % section)
            configsection = self.configdict[section]
            configsection.printSection(parameter=parameter,ntab=0)
        else:
            for section,configsection in self.configdict.iteritems():
                print 'Section %s:' % section
                configsection.printSection(ntab=1)

    def configure(self,section=None,parameter=None,configfile=None):
        if section is not None:
            if section not in self.configdict.keys():
                raise KeyError('Section "%s" not found in config file.' % section)
            configsection = self.configdict[section]
            configsection.configure(parameter=parameter)
        else:
            for section,configsection in self.configdict.iteritems():
                print 'Configuring %s:' % section
                configsection.configure(parameter=parameter)
        if configfile is not None:
            self.save(configfile)
            
    def save(self,configfile):
        #reconstruct the dictionary
        cdict = {}
        for section in self.configdict.keys():
            configsection = self.configdict[section]
            sectiondict = {}
            for name,configitem in configsection.configdict.iteritems():
                value = configitem.getValue()
                name = configitem.getName()
                caption = configitem.getCaption()
                default = configitem.getDefault()
                crange = configitem.getRange()
                ctype = configitem.getType()
                itemdict = {'name':name,'caption':caption,'value':value,'ctype':ctype}
                if default is not None:
                    itemdict['default'] = default
                if crange is not None:
                    itemdict['crange'] = crange
                sectiondict[name] = itemdict.copy()
            cdict[section] = sectiondict
        f = open(configfile,'wt')
        json.dump(cdict,f)
        f.close()
        print 'Saved configuration to %s' % configfile

class ConfigCommand(cmd.Cmd):
    
    def setConfig(self,configsection):
        self.section = configsection

    def do_quit(self,line):
        "Exit the program"
        return True

    def do_EOF(self,line):
        "Exit the program"
        return True
    
    def do_list(self,foo):
        """
        List all of the parameters for the given program.
        """
        ncols = 4
        paramlist = self.section.getParamNames()
        maxlen = max([len(p) for p in paramlist]) + 2
        nrows = len(paramlist)/ncols
        if len(paramlist) % ncols:
            nrows += 1
        print
        for i in range(0,nrows):
            row = []
            if len(paramlist) >= ncols:
                tcols = ncols
                for j in range(0,ncols):
                    row.append(paramlist.pop(0))
            else:
                tcols = len(paramlist)
                row = paramlist[:]
            fstr = '%-' + '%is' % maxlen
            fmt = fstr*tcols
            rowstr = fmt % tuple(row)
            print rowstr
        print

    def do_info(self,param):
        """List parameter information"""
        print
        self.section.printSection(parameter=param)

    def do_config(self,param):
        """Configure parameter information"""
        ci = self.section.getConfigItem(param) 
        if ci.getType() not in ['dict','list']:
            ci.configure()
        else:
            print 'Configuration of dictionaries and lists not yet supported.'

    def complete_config(self,text,line,begidx,endidx):
        params = self.section.getParamNames()
        if not text:
            completions = params[:]
        else:
            completions = [f for f in params if f.startswith(text)]
        return completions
            
    def complete_info(self,text,line,begidx,endidx):
        params = self.section.getParamNames()
        if not text:
            completions = params[:]
        else:
            completions = [f for f in params if f.startswith(text)]
        return completions

    
def main(args):
    if args.translate:
        translate(args.file,args.translate)
        sys.exit(0)
    if args.validate:
        try:
            configdict = json.load(open(args.file,'rt'))
        except Exception,obj:
            print 'Failed to load JSON format config file into data structure: "%s"' % obj.message
            sys.exit(1)
        try:
            config = Config(configdict)
        except Exception,obj:
            print 'Error validating configuration data: "%s"' % obj.message
            sys.exit(1)
        print '%s is a valid ShakeMap config file.' % args.file
        sys.exit(0)

    
    configdict = json.load(open(args.file,'rt'))
    config = Config(configdict)

    mycmd = ConfigCommand()
    mycmd.setConfig(config.getSection('grind'))
    mycmd.cmdloop('? to get help, Ctrl-D to quit')
    print
    sys.exit(0)
    
    section = None
    param = None
    if args.section:
        section = args.section
    if args.param:
        param = args.param
    if args.configure:
        if args.debug:
            config.configure(section=section,parameter=param)
        else:
            if DEBUG:
                cpath,cfile = os.path.split(args.file)
                cbase,cext = os.path.splitext(cfile)
                newfile = os.path.join(cpath,cbase+'_mod'+cext) 
            else:
                newfile = args.file
            config.configure(section=section,parameter=param,configfile=newfile)
        sys.exit(0)
    config.printConfig(section=section,parameter=param)

        
if __name__ == '__main__':
    desc = """Inspect or modify ShakeMap configuration file.
    """
    parser = argparse.ArgumentParser(description=desc)
    #positional arguments
    parser.add_argument('file',help='Input config JSON file')
    #optional arguments
    parser.add_argument('-s','--section', dest='section', 
                        help='List or configure section SECTION')
    # parser.add_argument('-p','--param', dest='param', 
    #                     help='List or configure parameter PARAMETER')
    # parser.add_argument('-c','--configure',  action='store_true',default=False,
    #                     help='configure sections and parameters')
    # parser.add_argument('-d','--debug',  action='store_true',default=False,
    #                     help='Debug ')
    parser.add_argument('-t','--translate', dest='translate', metavar='JSONFILE',
                        help='Translate ShakeMap 3.x config files into json format JSONFILE')
    parser.add_argument('-v','--validate',  action='store_true',default=False,
                        help='Validate input JSON format config file')
    pargs = parser.parse_args()
    main(pargs)

