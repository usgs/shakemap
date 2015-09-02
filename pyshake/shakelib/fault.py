#!/usr/bin/env python

#stdlib modules
import StringIO
import sys

#third party imports
import numpy as np
from openquake.hazardlib.geo import mesh
from ecef import latlon2ecef
from ecef import Vector

#CONSTANTS
#what is the maximum ratio of distance out of the plane defined by 3 points a 4th point can be before
#being considered non-co-planar?
OFFPLANE_TOLERANCE = 0.05

class FaultException(Exception):
    """
    Class to represent errors in the Fault class.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    

class Fault(object):
    """
    Class to handle fault files of various types and output fault data in various ways.
    """
    def __init__(self,lon,lat,depth,reference):
        """
        Constructor for Fault.
        :param lon:
            float longitude of hypocenter.
        :param lat:
            float latitude of hypocenter.
        :param depth:
            float depth of hypocenter.
        :param reference:
            string citeable reference for Fault.
        """
        self.lon = lon
        self.lat = lat
        self.depth = depth
        self.reference = reference
        self._validate()
        self.setQuadrilaterals()
        
    @classmethod
    def readFaultFile(cls,faultfile):
        """
        Read fault file format as defined in ShakeMap Software Guide.  
        :param faultfile: 
            Path to fault file OR file-like object in GMT psxy format, where
             * Fault vertices are space separated lat,lon,depth triplets on a single line.
             * Fault segments are separated by lines containing ">"
             * Fault segments must be closed.
             * Fault segments must be all clockwise or all counter-clockwise.
        :returns:
           Fault object.
        :raises Exception:
            when any of above conditions are not met.
        """
        x = []
        y = []
        z = []
        if isinstance(faultfile,str) or isinstance(faultfile,unicode):
            faultfile = open(faultfile,'rt')
            faultlines = faultfile.readlines()
        else:
            faultlines = faultfile.readlines()
        reference = ''
        for line in faultlines:
            sline = line.strip()
            if sline.startswith('#'):
                reference += sline
                continue
            if sline.startswith('>'):
                if len(x): #start of new line segment
                    x.append(np.nan)
                    y.append(np.nan)
                    z.append(np.nan)
                    continue
                else: #start of file
                    continue
            parts = sline.split()
            if len(parts) < 3:
                raise Exception('Finite fault file %s has no depth values.' % faultfile)
            y.append(float(parts[0]))
            x.append(float(parts[1]))
            z.append(float(parts[2]))
        faultfile.close()
        if np.isnan(x[-1]):
            x = x[0:-1]
            y = y[0:-1]
            z = z[0:-1]

        return cls(x,y,z,reference)

    def getQuadrilaterals(self):
        """
        Return a list of quadrilaterals.
        :returns:
            each quad is a tuple of four Point objects
            (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py)
        """
        return self.Quadrilaterals

    def getQuadWidth(self,p0,p1,p3):
        """
        Return width of an individual planar trapezoid, where the p0-p1 distance represents the long side.
        :param p0: 
            ECEF x,y,z point representing the first vertex of a quadrilateral.
        :param p1: 
            ECEF x,y,z point representing the second vertex of a quadrilateral.
        :param p3: 
            ECEF x,y,z point representing the fourth vertex of a quadrilateral.
        :returns:
           width of planar trapezoid.
        """
        AB = p0-p1
        AC = p0-p3
        numer = np.cross(np.cross(AB,AC),AB)
        denom = numer.norm()
        width = numer/denom
        return width

    def getStrike(self):
        """
        Return representative strike angle from first vertex of first quad to second vertex of last quad.
        :returns:
          float strike angle
        """
        P0 = self.Quadrilaterals[0][0]
        P1 = self.Quadrilaterals[-1][1]
        return P0.azimuth(P1)

    def getTopOfRupture(self):
        """
        Determine shallowest vertex of entire fault.
        :returns:
            float shallowest depth of all vertices.
        """
        mindep = 9999999
        for quad in self.Quadrilaterals:
            P0,P1,P2,P3 = quad
            depths = np.array([P0.depth,P1.depth,P2.depth,P3.depth])
            if np.min(depths) < mindep:
                mindep = np.min(depths)
        return mindep
    
    def getDip(self):
        """
        Return average dip of all quadrilaterals in the fault.
        :returns:
           Average dip in degrees.
        """
        dipsum = 0.0
        for quad in self.Quadrilaterals:
            P0,P1,P2,P3 = quad
            d1 = P1.depth * -1
            d2 = P2.depth * -1
            dx = (P2-P1).mag()
            dz = d1 - d2
            dip = np.degrees(np.tan(dx/dz))
            dipsum += dip
        dip = dipsum/len(self.Quadrilaterals)
        return dip
    
    def getWidth(self):
        """
        Return the average fault width (km) for all quadrilaterals defined for the fault.
        :returns:
            Average width in km of all fault quadrilaterals.
        """
        wsum = 0.0
        for quad in self.Quadrilaterals:
            P0,P1,P2,P3 = quad
            p0 = ecef.Vector.fromPoint(P0)
            p1 = ecef.Vector.fromPoint(P1)
            p3 = ecef.Vector.fromPoint(P3)
            wsum += self.getQuadWidth(p0,p1,p3)
        mwidth = (wsum/len*self.Quadrilaterals)/1000.0
        return mwidth
    
    def getTrapMeanLength(self,p0,p1,p2,p3):
        """
        Return the sqrt of the area of a quadrilateral (used for QA of fault plane).
        :param p0: 
            ECEF x,y,z point representing the first vertex of a quadrilateral.
        :param p1: 
            ECEF x,y,z point representing the second vertex of a quadrilateral.
        :param p2: 
            ECEF x,y,z point representing the third vertex of a quadrilateral.
        :param p3: 
            ECEF x,y,z point representing the fourth vertex of a quadrilateral.
        :returns:
            square root of trapezoid area.
        """
        #area of a trapezoid: A = (a+b)/2 * h (https://en.wikipedia.org/wiki/Trapezoid)
        h = self.getQuadWidth(p0,p1,p3)
        a = (p1-p0).mag
        b = (p2-p3).mag
        A = ((a+b)/2.0)*h
        length = np.sqrt(A)
        return length
    
    def getDistanceToPlane(self,planepoints,otherpoint):
        """
        Calculate a point's distance to a plane.  Used to figure out if a quadrilateral points are all co-planar.
        :param planepoints:
            List of three points (Vector objects) defining a plane.
        :param otherpoint:
            4th Vector to compare to points defining the plane
        :returns:
            Distance (in meters) from otherpoint to plane.
        """
        #from https://en.wikipedia.org/wiki/Plane_(geometry)#Describing_a_plane_through_three_points
        p0,p1,p2 = planepoints
        x1,y1,z1 = p0.getArray()
        x2,y2,z2 = p1.getArray()
        x3,y3,z3 = p2.getArray()
        D = np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
        d = -1
        at = np.array([[1,y1,z1],[1,y2,z2],[1,y3,z3]])
        bt = np.array([[x1,1,z1],[x2,1,z2],[x3,1,z3]])
        ct = np.array([[x1,y1,1],[x2,y2,1],[x3,y3,1]])
        a = np.dot((d/D),at)
        b = np.dot((d/D),bt)
        c = np.dot((d/D),ct)

        numer = np.abs(a*otherpoint.x + b*otherpoint.y + c*otherpoint.z + d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer/denom
        return dist

    
    def validateQuad(self,P0,P1,P2,P3):
        """
        Validate and fix* a given quadrilateral (*currently "fix" means check third vertex for co-planarity
        with other three points, and force it to be co-planar if it's not wildly out of the plane.()
        :param P0:
            First vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P1:
            Second vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P2:
            Third vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :param P3:
            Fourth vertex https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/point.py
        :returns:
           Tuple of (potentially) modified vertices.
        :raises FaultException:
           * if top and bottom edges are not parallel to surface
           * if dip angle is not dipping to the right relative to strike (defined by first two vertices)
           * if all 4 points are not reasonably co-planar (P2 is more than 5% of mean length of trapezoid out of plane)
        """
        #TODO: Someday fix the rule about dip angle being clockwise and 0-90 degrees
        #In theory, you could flip the quadrilateral by 180 degrees and it would be ok.
        
        #Are the top and bottom edges both parallel to the surface?
        topDepthsEqual = P0.depth == P1.depth
        bottomDepthsEqual = P2.depth == P3.depth
        if not topDepthsEqual or not bottomDepthsEqual:
            raise FaultException('Top and bottom edges of fault quadrilateral must be parallel to the surface')
        #Is top edge defined by first two vertices?
        #Is dip angle clockwise and btw 0-90 degrees?
        if P1.depth < P2.depth:
            raise FaultException('Top edge of a quadrilateral must be defined by the first two vertices')
        #Are all 4 points (reasonably) co-planar?
        #Translate vertices to ECEF
        p0 = Vector.fromPoint(P0.lat,P0.lon,P0.depth)
        p1 = Vector.fromPoint(P1.lat,P1.lon,P1.depth)
        p2 = Vector.fromPoint(P2.lat,P2.lon,P2.depth)
        p3 = Vector.fromPoint(P3.lat,P3.lon,P3.depth)
        #Calculate normalized vector along top edge
        v0 = (P1-P0).norm()
        #Calculate distance btw p3 and p2
        d = (p3-p2).mag()
        #get the new P2 value
        v1 = v0*d
        newp2 = p3 + v1
        planepoints = [p0,p1,p2]
        dnormal = pself.getDistanceToPlane(planepoints,p2)
        geometricMean = self.getTrapMeanLength(p0,p1,newp2,p3)
        if dnormal/geometricMean > OFFPLANE_TOLERANCE:
            raise FaultException('Points in quadrilateral are not co-planar')
        newP0 = p0.toPoint()
        newP1 = p1.toPoint()
        newP2 = newp2.toPoint()
        newP3 = p3.toPoint()
        return (newP0,newP1,newP2,newP3)
        
    
    def setQuadrilaterals(self):
        """
        Create internal list of N quadrilaterals.
        """
        #Fault QA rules
        #1) Fault must consist of 1 or more quadrilaterals, where each quad top/bottom edges are
        #   parallel to the surface
        #2) The strike angle of each quadrilateral is defined by the first two vertices of that quad
        #3) The dip angle is defined by segments 2 and 3, or 1 and 4.  This angle must be clockwise with respect to the 
        #   strike angle, and between 0 and 90 degrees.
        #4) The top edge of each quad must be defined by the first two vertices of that quad.
        #5) 4 points of quadrilateral must be co-planar
        lon = np.array(lon)
        lat = np.array(lat)
        depth = np.array(depth)
        inan = np.isnan(lon)
        numnans = len(lon[inan])
        numsegments = numnans + 1
        #requirements:
        # 1) Coordinate arrays must be same length
        # 2) Polygons must be quadrilaterals
        # 3) Quads must be closed
        # 4) Quads must be planar
        if len(lon) != len(lat) != len(depth):
            raise IndexError('Length of input lon,lat,depth arrays must be equal')
        
        istart = 0
        endpoints = list(np.where(np.isnan(lon))[0])
        endpoints.append(len(lon))
        self.Quadrilaterals = []
        for iend in endpoints:
            lonseg = lon[istart:iend]
            latseg = lat[istart:iend]
            depthseg = depth[istart:iend]
            #each segment can have many contiguous quadrilaterals defined in it
            #separations (nans) between segments mean that segments are not contiguous.
            npoints = len(lonseg)
            nquads = ((npoints - 5)/2) + 1
            ioff = 0
            for i in range(0,nquads):
                endidx = ioff-1 #we have the closing vertex that we're not interested in here
                topLeft = point.Point(lonseg[ioff],latseg[ioff],depthseg[ioff])
                topRight = point.Point(lonseg[ioff+1],latseg[ioff+1],depthseg[ioff])
                bottomRight = point.Point(lonseg[endidx-2],latseg[endidx-2],depthseg[endidx-2])
                bottomLeft = point.Point(lonseg[endidx-1],latseg[endidx-1],depthseg[endidx-1])
                surface = self.validateQuad(topLeft,topRight,bottomRight,bottomLeft)
                self.Quadrilaterals.append(surface)
                ioff += 1
            istart = iend+1
            
        
    def getReference(self):
        """
        Return whatever reference information was contained in fault file.
        :returns:
           string citeable reference
        """
        return self.reference
        
    def getNumSegments(self):
        """
        Return a count of the number of fault segments.
        :returns:
            number of fault segments
        """
        return len(np.where(np.isnan(self.x))[0]) + 1

    def getFaultAsArrays(self):
        """
        Return a 3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.  Fault segments are separated by 
        numpy.NaN values.
        :returns:
            3-tuple of numpy arrays indicating X,Y,Z (lon,lat,depth) coordinates.
        """
        return (np.array(self.x),np.array(self.y),np.array(self.z))

    def getFaultAsMesh(self):
        """
        Return fault segments as a OQ-Hazardlib Mesh object. 
        :returns:
            Mesh (https://github.com/gem/oq-hazardlib/blob/master/openquake/hazardlib/geo/mesh.py)
        """
        fault = mesh.Mesh(self.x,self.y,self.z)
        return fault

    def _validate(self):
        """
        Ensure that all segments are closed.
        :raises Exception:
            if unclosed segments exist.
        """
        #TODO - implement ccw algorithm...
        #http://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        if len(self.lon) != len(self.lat) or len(self.lon) != len(self.depth):
            raise Exception("Fault coordinates don't match")
        inan = np.where(np.isnan(self.lon))[0]
        if not np.isnan(self.lon[inan[-1]]):
            inan = list(inan).append(len(self.lon))
        istart = 0
        for i in range(0,len(inan)):
            iend = inan[i]-1
            x1 = self.lon[istart]
            x2 = self.lon[iend]
            y1 = self.lat[istart]
            y2 = self.lat[iend]
            z1 = self.depth[istart]
            z2 = self.depth[iend]
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

