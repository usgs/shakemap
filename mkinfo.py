#!/usr/bin/env python

from xml.dom import minidom
import sys
import os.path
import glob
import re
import datetime

#Mon Jun 4 05:35:03 MDT 2012
TIMEFMT = '%a %b %Y %H:%M:%S'

def getEventTime(event):
    year = int(event.getAttribute('year'))
    month = int(event.getAttribute('month'))
    day = int(event.getAttribute('day'))
    hour = int(event.getAttribute('hour'))
    minute = int(event.getAttribute('minute'))
    second = int(event.getAttribute('second'))
    etime = datetime.datetime(year,month,day,hour,minute,second).strftime(TIMEFMT)
    return etime + str(event.getAttribute('timezone'))

def parseStations(stationdom):
    nsensors = 0
    ndyfi = 0
    nmacro = 0
    sources = []
    netids = []
    for station in stationdom.getElementsByTagName('station'):
        source = station.getAttribute('source').lower()
        netid = station.getAttribute('netid').lower()
        if source not in sources:
            sources.append(source)
        if source not in netids:
            netid.append(netid)
        if netid == 'ciim' or source == 'dyfi':
            ndyfi += 1
            continue
        if netid == 'intensity':
            nmacro += 1
            continue
        nsensors += 1
    return nsensors,ndyfi,nmacro,sources,netids
        

def parseInfo(infodom):
    infodict = {}
    descdict = {}
    for tag in infodom.getElementsByTagName('tag'):
        key = tag.getAttribute('name')
        value = tag.getAttribute('value')
        desc = tag.getAttribute('desc')
        infodict[key] = value
        descdict[key] = desc
    return (infodict,descdict)

if __name__ == '__main__':

    helptxt = """
    ################################################################################
    # Program     : mkinfo
    # Description : Program to create summary info.html file, pulling information from various 
    #               ShakeMap input/output files.
    # Options     :
    #     Flag       Arg                            Description
    #     ---------- ---------                      -----------------------------------------------------
    #    -event      event_id                       Specifies the id of the event to process
    #     ---------- ---------                      -----------------------------------------------------
    #    -help                                      Print program documentation and quit.
    ################################################################################
    """
    
    allflags = ['-event','-help']
    #Make sure the INI file exists first
    homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    args = sys.argv[1:]

    if '-help' in args or len(args) == 0:
        print helptxt
        sys.exit(0)

    if '-event' in args:
        idx = args.index('-event')
        if idx < len(args)-1 and not args[idx+1].startswith('-'):
            eventID = args[idx+1]
        else:
            print 'No option argument found for "-event". Exiting.\n\n'
            print helptxt
            sys.exit(0)
    
        
    
    eventfolder = os.path.join(homedir,'../data/',eventID)
    #eventfolder = os.path.join('/opt/local/ShakeMap/data/',eventID)

    if not os.path.isdir(eventfolder):
        print 'Could not locate event directory %s. Exiting.' % eventfolder
        print helptxt
        sys.exit(0)
    
    #output file
    infohtml = os.path.join(eventfolder,'output','info.html')

    #input files
    infoxml = os.path.join(eventfolder,'output','info.xml')
    eventxml = os.path.join(eventfolder,'input','event.xml')
    zonexml = os.path.join(eventfolder,'zoneconfig','zoneinfo.xml')
    gridxml = os.path.join(eventfolder,'output','grid.xml')
    faultfile = glob.glob(os.path.join(eventfolder,'input','*fault*.txt'))
    if len(faultfile):
        faultfile = faultfile[0]
    else:
        faultfile = 'None'
    stationxml = os.path.join(eventfolder,'output','stationlist.xml')

    #info, event, grid files must exist
    reqfiles = [infoxml,eventxml,gridxml]
    optfiles = [zonexml,stationxml]
    missing = []
    for afile in reqfiles:
        if not os.path.isfile(afile):
            missing.append(afile)
    if len(missing):
        print 'Missing one of the following required input files!'
        for m in missing:
            print '\t' + m
        sys.exit(0)

    #check to see if optional files exist
    optfiles = [zonexml,stationxml]

    
    infodom = minidom.parse(infoxml)
    eventdom = minidom.parse(eventxml)
    zonedom = minidom.parse(zonexml)
    stationdom = minidom.parse(stationxml)
    griddom = minidom.parse(gridxml)

    #get stuff from the event.xml DOM
    event = eventdom.getElementsByTagName('earthquake')[0]
    eventid = event.getAttribute('network')+event.getAttribute('id')
    etime = getEventTime(event)
    
    eqinfo = zonedom.getElementsByTagName('eqinfo')[0]
    location = event.getAttribute('locstr')

    nsensors,ndyfi,nmacro,sources,netids = parseStations(stationdom)

    infodict,descdict = parseInfo(infodom)
    processtime = infodict['grind_time'].strip()
    grid = griddom.getElementsByTagName('shakemap_grid')[0]
    shakeversion = grid.getAttribute('shakemap_version')
    codeversion = grid.getAttribute('code_version')
    eventel = griddom.getElementsByTagName('grid_specification')[0]
    latspan = float(eventel.getAttribute('lat_max')) - float(eventel.getAttribute('lat_min'))
    lonspan = float(eventel.getAttribute('lon_max')) - float(eventel.getAttribute('lon_min'))
    latspace = float(eventel.getAttribute('nominal_lat_spacing'))
    lonspace = float(eventel.getAttribute('nominal_lon_spacing'))
    nlat = float(eventel.getAttribute('nlat'))
    nlon = float(eventel.getAttribute('nlon'))
    lat = float(eqinfo.getAttribute('lat'))
    lon = float(eqinfo.getAttribute('lon'))
    depth = float(eqinfo.getAttribute('depth'))
    mag = float(eqinfo.getAttribute('magnitude'))
    focal = zonedom.getElementsByTagName('focalMechanism')[0].firstChild.data
    mechsource = zonedom.getElementsByTagName('momentTensorSource')[0].firstChild.data
    feregion = int(zonedom.getElementsByTagName('feregion')[0].firstChild.data)
    fename = zonedom.getElementsByTagName('fename')[0].firstChild.data
    regime = zonedom.getElementsByTagName('earthquakeType')[0].firstChild.data
    eq2 = zonedom.getElementsByTagName('equations')[0].getAttribute('eq2')
    eq3a = zonedom.getElementsByTagName('equations')[0].getAttribute('eq3a')
    eq3b = zonedom.getElementsByTagName('equations')[0].getAttribute('eq3b')
    gmpe = infodict['GMPE']
    if infodict.has_key('IPE'):
        ipe = infodict['IPE']
    else:
        ipe = 'None'
    pgm2mi = infodict['pgm2mi']
    mi2pgm = infodict['mi2pgm']
    iscale = infodict['miscale']
    pga_max = float(infodict['pga_max'])
    pgv_max = float(infodict['pgv_max'])
    mmi_max = float(infodict['mi_max'])
    psa03_max = float(infodict['psa03_max'])
    psa30_max = float(infodict['psa30_max'])
    psa10_max = float(infodict['psa10_max'])
    map_bound = infodict['map_bound']
    try:
        meanerror = float(infodict['mean_uncertainty'])
    except ValueError:
        meanerror = float('nan')
    grade = infodict['grade']
    distance_used = infodict['median_dist']
    biaslist = [float(b) for b in infodict['bias'].split()]
    biaskeystr = descdict['bias']
    pattern = "\\(.+?\\)"
    m = re.search(pattern,biaskeystr)
    biaskeys = biaskeystr[m.start()+1:m.end()-1].strip().split()
    biasdict = dict(zip(biaskeys,biaslist))
    biasnuggets = []
    for key,value in biasdict.iteritems():
        biasnuggets.append('%s = %.2f' % (key,value))
    mmibias = float(infodict['mi_bias'])
    
    
    #open the output HTML file
    f = open(infohtml,'wt')
    f.write('<html>\n')
    f.write('<header><title>ShakeMap Information for Event %s</title></head>\n' % eventid)

    #event/source info table
    f.write('<table border="1" width="800">\n')
    f.write('<tr><td>\n')
    f.write('\t<table border="1" width="400">\n')
    f.write('\t<tr><th bgcolor="#8db3e2">Event/Source Info (Event parameters and ShakeMap run information)</th></tr>\n')
    f.write('\t<tr><td>ID: %s</td></tr>\n' % eventid)
    f.write('\t<tr><td>Date/Time: %s</td></tr>\n' % etime)
    f.write('\t<tr><td>Location: %s</td></tr>\n' % location)
    f.write('\t<tr><td>Process Time: %s</td></tr>\n' % processtime)
    f.write('\t<tr><td>ShakeMap Version: %s</td></tr>\n' % shakeversion)
    f.write('\t<tr><td>Code Version: %s</td></tr>\n' % codeversion)
    f.write('\t<tr><td>Hypocenter and Magnitude: lat=%.4f lon=%.4f depth=%.1f magnitude=%.1f</td></tr>\n' % (lat,lon,depth,mag))
    f.write('\t<tr><td>Focal Mechanism: %s</td></tr>\n' % focal)
    f.write('\t<tr><td>Mechanism Source: %s</td></tr>\n' % mechsource)
    f.write('\t<tr><td>Flinn Engdahl Region: %i - "%s"</td></tr>\n' % (feregion,fename))
    f.write('\t<tr><td>Tectonic Regime: %s</td></tr>\n' % regime)
    f.write('\t<tr><td>Finite Fault Model: %s</td></tr>\n' % faultfile)
    f.write('\t<tr><td>Finite Fault Reference: %s</td></tr>\n' % faultfile)
    f.write('\t</table>\n')
    f.write('</td>\n')

    #Mapping info table
    f.write('<td>\n')
    f.write('\t<table border="1" width="400">\n')
    f.write('\t<tr><th bgcolor="#8db3e2">Mapping Info (Mapping and Grid Parameters)</th></tr>\n')
    f.write('\t<tr><td>Latitude Span: %.4f Longitude Span: %.4f</td></tr>\n' % (latspan,lonspan))
    f.write('\t<tr><td>Map Bound: %s</td></tr>\n' % map_bound)
    f.write('\t<tr><td>Latitude Grid Interval: %.4f</td></tr>\n' % latspace)
    f.write('\t<tr><td>Longitude Grid Interval: %.4f</td></tr>\n' % lonspace)
    f.write('\t<tr><td>Grid nx: %i</td></tr>\n' % nlon)
    f.write('\t<tr><td>Grid ny: %i</td></tr>\n' % nlat)
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t</table>\n')
    f.write('</td></tr>\n')
    

    #gmpe info table
    f.write('<tr><td>\n')
    f.write('\t<table border="1" width="400">\n')
    f.write('\t<tr><th bgcolor="#8db3e2">GMPE Info (Ground Motion Estimation Information)</th></tr>\n')
    f.write('\t<tr><td>Ground Motion Prediction Equation: %s</td></tr>\n' % gmpe)
    f.write('\t<tr><td>Intensity Prediction Equation: %s</td></tr>\n' % ipe)
    f.write('\t<tr><td>Ground Motion to Intensity: %s</td></tr>\n' % pgm2mi)
    f.write('\t<tr><td>Inverse Intensity Function: %s</td></tr>\n' % mi2pgm)
    f.write('\t<tr><td>Intensity Scale: %s</td></tr>\n' % iscale)
    f.write('\t<tr><td>Site Correction Regime: %s</td></tr>\n' % '')
    f.write('\t<tr><td>Site Correction File: %s</td></tr>\n' % '')
    f.write('\t<tr><td>Median Distance Used?: %s</td></tr>\n' % distance_used)
    f.write('\t<tr><td>Magnitude Bias: %s </td></tr>\n' % ','.join(biasnuggets))
    f.write('\t<tr><td>Intensity Magnitude Bias: %.2f </td></tr>\n' % mmibias)
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t<tr><td>&nbsp;</td></tr>\n')
    f.write('\t</table>\n')
    f.write('</td>\n')    

    #Station info table
    f.write('<td>\n')
    f.write('\t<table border="1" width="400">\n')
    f.write('\t<tr><th bgcolor="#8db3e2">Mapping Info (Mapping and Grid Parameters)</th></tr>\n')
    f.write('\t<tr><td>Seismic Data Source(s): %s</td></tr>\n' % ','.join(sources))
    f.write('\t<tr><td>Seismic Network ID(s): %s</td></tr>\n' % ','.join(netids))
    f.write('\t<tr><td>Number of Seismic Stations: %i</td></tr>\n' % nsensors)
    f.write('\t<tr><td>Number of DYFI Stations: %i</td></tr>\n' % ndyfi)
    f.write('\t<tr><td>Number of Macroseismic Stations: %i</td></tr>\n' % nmacro)
    f.write('\t<tr><td>PGA Max: %.2f</td></tr>\n' % pga_max)
    f.write('\t<tr><td>PGV Max: %.2f</td></tr>\n' % pgv_max)
    f.write('\t<tr><td>MMI Max: %.2f</td></tr>\n' % mmi_max)
    f.write('\t<tr><td>PSA 0.3 sec Max: %.2f</td></tr>\n' % psa03_max)
    f.write('\t<tr><td>PSA 10 sec Max: %.2f</td></tr>\n' % psa10_max)
    f.write('\t<tr><td>PSA 3 sec Max: %.2f</td></tr>\n' % psa30_max)
    f.write('\t<tr><td>Mean Uncertainty: %.2f</td></tr>\n' % meanerror)
    f.write('\t<tr><td>Uncertainty Grade: %s</td></tr>\n' % grade)
    f.write('\t<tr><td>Does Focal Mechanism Satisfy Interface Conditions?: %s</td></tr>\n' % eq2)
    f.write('\t<tr><td>Is Depth Within Interface Depth Interval?: %s</td></tr>\n' % eq3a)
    f.write('\t<tr><td>Is Depth Within Intraslab Depth Interval?: %s</td></tr>\n' % eq3b)
    f.write('\t</table>\n')
    f.write('</td></tr>\n')
    
    f.write('</html>\n')
    f.close()

    infodom.unlink()
    eventdom.unlink()
    zonedom.unlink()
    stationdom.unlink()
