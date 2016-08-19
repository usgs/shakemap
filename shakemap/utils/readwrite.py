#!/usr/bin/env python

import numpy as np
from lxml import etree


def read_nshmp_fault_xml(file):
    """
    Method for reading the XML format used by the NSHMP for faults. 
    `[example] <https://github.com/emthompson-usgs/nshmp-model-cous-2014/blob/master/Western%20US/Fault/Geologic%20Model%20Full%20Rupture.xml>`__

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
