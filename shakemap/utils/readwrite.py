#!/usr/bin/env python

import numpy as np
from lxml import etree


def read_nshmp_fault_xml(file):
    """
    Method for reading the XML format used by the NSHMP for faults. 
    `[example] <https://github.com/usgs/nshmp-model-cous-2014/blob/master/Western%20US/Fault/Geologic%20Model%20Full%20Rupture.xml>`__

    :param file:
        An XML file in the format used for the 2014 USGS NSHMP. 

    :returns:
        A list with length equal to the number of faults; 
        each element contains the a dictionary of fault information, which is 
        stored in a dictionary with two elements: "Geometry" and 
        "IncrementalMfd". Geometry is another dictionary, which contains trace
        information and other paramters such as rake. IncrementalMfd is a list
        of dictionaries. 

    """

    tree = etree.parse(file)
    root = tree.getroot()
    sources = root.findall('Source')
    srclist = []

    for s in sources:
        ch = s.getchildren()

        #-----------------------------------------------
        # Child dictionary, includes Geometry dictionary
        # and IncrementalMfd list of dictionaries
        #-----------------------------------------------
        chd = {'Geometry':{}, 
               'IncrementalMfd':[]}

        for a in ch:
            if a.tag == 'Geometry':
                # Attributes
                for k,v in a.attrib.iteritems():
                    chd['Geometry'].update({k: v})

                # Children (Trace, LowerTrace)
                t = a.getchildren()

                for b in t:
                    # Turn text into arrays
                    tracedata = b.text
                    tlines = tracedata.split('\n')
                    lons = np.array([])
                    lats = np.array([])
                    deps = np.array([])

                    for line in tlines: 
                        if line is not '':
                            lon, lat, dep = line.split(',')
                            lons = np.append(lons, float(lon))
                            lats = np.append(lats, float(lat))
                            deps = np.append(deps, float(dep))
                    chd['Geometry'].update({b.tag: {'lon':lons,
                                                    'lat':lats,
                                                    'dep':deps}})

            elif a.tag == 'IncrementalMfd':
                # These only contain attributes
                imfd = {}

                for k,v in a.attrib.iteritems():
                    imfd.update({k: v})
                chd['IncrementalMfd'].append(imfd)

        srclist.append(chd)

    return srclist

def read_nshmp_grid_xml(file):
    """
    Method for reading the XML format used by the NSHMP for gridded seismicity.  
    `[example] <https://github.com/usgs/nshmp-model-cous-2014/blob/master/Central%20%26%20Eastern%20US/Grid/rlme/Charlevoix%20Seismic%20Zone.xml>`__

    This function works for both grided source, including RMLEs. 

    :param file:
        An XML file in the format used for the 2014 USGS NSHMP. 

    :returns:
        A dictionary with two entries: 'Settings' and 'Nodes'. 
        'Settings" is a dictionary, 'Nodes' is a list of dictionaries with
        length equal to the total number of RLME nodes. The Node dictionaries
        include elements for 'lat', 'lon', 'dep', 'rate', etc. Note that there
        are minor differences in these values depending on whether or not the
        gridded seismicity is for RLMEs or not. 

    """

    tree = etree.parse(file)
    root = tree.getroot()
    settings = root.findall('Settings')[0]
    nodes = root.findall('Nodes')[0].getchildren()

    #----------------------------------------------
    # Settings; includes 
    #    list of DefaultMfds and 
    #    dictionary of SourceProperties
    #----------------------------------------------
    setdict = {}
    for c in settings.getchildren():

        if c.tag == 'DefaultMfds':
            mfdlist = []

            for m in c.getchildren():
                mdict = {}
                for k,v in m.attrib.iteritems():
                    mdict.update({k: v})
                mfdlist.append(mdict)
            setdict.update({'DefaultMfds':mfdlist})

        if c.tag == 'SourceProperties':
            pdict = {}
            for k,v in c.attrib.iteritems():
                pdict.update({k: v})
            setdict.update({'SourceProperties':pdict})

    #----------------------------------------------
    # Nodes; list of dictionaries
    #----------------------------------------------
    nlist = []
    for n in nodes:
        ndict = {}
        for k,v in n.attrib.iteritems():
            ndict.update({k: v})
        loc = n.text.split(',')
        ndict['lat'] = float(loc[1])
        ndict['lon'] = float(loc[0])
        ndict['dep'] = float(loc[2])
        nlist.append(ndict)

    return {'Settings':setdict,'Nodes':nlist}
