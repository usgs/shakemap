#!/usr/bin/env python

#stdlib modules
import StringIO
import sys

#third party imports
import numpy as np
from openquake.hazardlib.geo import mesh

class Fault(object):
    """
    Class to handle fault files of various types and output fault data in various ways.
    """
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.reference = ''
        self.faultfile = ''

    def readFromGMTFormat(self,faultfile):
        """
        Read fault file format as defined in ShakeMap Software Guide.  
        Input: 
        Fault file in GMT psxy format, where
        
         * Fault vertices are space separated lat,lon,depth triplets on a single line.
         * Fault segments are separated by lines containing ">"
         * Fault segments must be closed.
         * Fault segments must be all clockwise or all counter-clockwise.

        Raises Exception when above conditions are not met.
        """
        if isinstance(faultfile,str) or isinstance(faultfile,unicode):
            self.faultfile = faultfile
            faultlines = open(faultfile,'rt').readlines()
        else:
            self.faultfile = 'File-like object'
            faultlines = faultfile.readlines()
        reference = ''
        for line in faultlines:
            sline = line.strip()
            if sline.startswith('#'):
                reference += sline
                continue
            if sline.startswith('>'):
                if len(self.x): #start of new line segment
                    self.x.append(np.nan)
                    self.y.append(np.nan)
                    self.z.append(np.nan)
                    continue
                else: #start of file
                    continue
            parts = sline.split()
            if len(parts) < 3:
                raise Exception('Finite fault file %s has no depth values.' % self.faultfile)
            self.y.append(float(parts[0]))
            self.x.append(float(parts[1]))
            self.z.append(float(parts[2]))

        self.reference = reference
        if np.isnan(self.x[-1]):
            self.x = self.x[0:-1]
            self.y = self.y[0:-1]
            self.z = self.z[0:-1]
        self._validate()
        faultfile.close()

    def getReference(self):
        """
        Return whatever reference information was contained in fault file.
        """
        return self.reference
        
    def getNumSegments(self):
        """
        Return a count of the number of fault segments.
        """
        return len(np.where(np.isnan(self.x))[0]) + 1

    def getFaultAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.  Fault segments are separated by 
        numpy.NaN values.
        """
        return (np.array(self.x),np.array(self.y),np.array(self.z))

    def getFaultAsMesh(self):
        """
        Return fault segments as a OQ-Hazardlib Mesh object. 
        
        https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py
        """
        fault = mesh.Mesh(self.x,self.y,self.z)
        return fault

    def _validate(self):
        #TODO - implement ccw algorithm...
        #http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        if len(self.x) != len(self.y) or len(self.x) != len(self.z):
            raise Exception("Fault coordinates don't match")
        inan = np.where(np.isnan(self.x))[0]
        if not np.isnan(self.x[inan[-1]]):
            inan = list(inan).append(len(self.x))
        istart = 0
        for i in range(0,len(inan)):
            iend = inan[i]-1
            x1 = self.x[istart]
            x2 = self.x[iend]
            y1 = self.y[istart]
            y2 = self.y[iend]
            z1 = self.z[istart]
            z2 = self.z[iend]
            if x1 != x2 or y1 != y2 or z1 != z2:
                raise Exception('Unclose segments exist in fault file.')
            istart = inan[i]+1
        
def _test():
    #TODO - actually test against expected output
    samplefault = """#Hartzell, S. (pers. comm., 2011)
    >
    30.685 103.333 0
    31.566 104.277 0
    31.744 104.058 15
    30.856 103.115 15
    30.685 103.333 0
    >
    30.762 103.237 0
    31.643 104.181 0
    31.781 104.006 21
    30.897 103.064 21
    30.762 103.237 0
    >
    31.610 104.232 0
    32.815 105.562 0
    32.901 105.460 16
    31.694 104.122 16
    31.610 104.232 0
    >"""
    newfaultlist = []
    for line in samplefault.split('\n'):
        newfaultlist.append(line.strip())
    newfault = '\n'.join(newfaultlist)
    faultfile = StringIO.StringIO(newfault)
    f = Fault()
    f.readFromGMTFormat(faultfile)
    (x,y,z) = f.getFaultAsArrays()
    
    print 'Fault has %i segments' % f.getNumSegments()
    print 'Fault points:'
    for i in range(0,len(x)):
        print y[i],x[i],z[i]
    
if __name__ == '__main__':
    if len(sys.argv) > 1:
        faultfile = sys.argv[1]
        fault = Fault()
        fault.readFromGMTFormat(faultfile)
        print 'Fault has %i segments' % fault.getNumSegments()
        print 'Fault points:'
        (x,y,z) = fault.getFaultAsArrays()
        for i in range(0,len(x)):
            print y[i],x[i],z[i]
    else:
        _test()

