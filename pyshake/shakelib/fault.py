#!/usr/bin/env python

#stdlib modules
import StringIO
import sys
import copy

#third party imports
import numpy as np
from openquake.hazardlib.geo import mesh
from openquake.hazardlib.geo import point,utils
from ecef import latlon2ecef
from vector import Vector
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
            sequence of fault longitude vertices.
        :param lat:
            sequence of fault latitude vertices.
        :param depth:
            sequence of fault depth vertices.
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
        isFile = False
        if isinstance(faultfile,str) or isinstance(faultfile,unicode):
            isFile = True
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
            if not len(sline.strip()):
                continue
            parts = sline.split()
            if len(parts) < 3:
                raise Exception('Finite fault file %s has no depth values.' % faultfile)
            y.append(float(parts[0]))
            x.append(float(parts[1]))
            z.append(float(parts[2]))
        if isFile:
            faultfile.close()
        if np.isnan(x[-1]):
            x = x[0:-1]
            y = y[0:-1]
            z = z[0:-1]

        return cls(x,y,z,reference)

    def multiplyFaultLength(self,factor):
        for i in range(0,len(self.Quadrilaterals)):
            quad = self.Quadrilaterals[i]
            P0,P1,P2,P3 = quad
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p2 = Vector.fromPoint(P2)
            p3 = Vector.fromPoint(P3)
            dtop = p1-p0
            toplen = dtop.mag()
            topnorm = dtop.norm()
            dbottom = p2-p3
            bottomlen = dbottom.mag()
            bottomnorm = dbottom.norm()
            newtoplen = toplen*factor
            newbottomlen = bottomlen*factor
            newp1 = p0+topnorm*newtoplen
            newp2 = p3+bottomnorm*newbottomlen
            newP1 = newp1.toPoint()
            newP2 = newp2.toPoint()
            newP1.depth = P0.depth
            newP2.depth = P3.depth
            self.Quadrilaterals[i][1] = newP1
            self.Quadrilaterals[i][2] = newP2
    
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
        t1 = (AB.cross(AC).cross(AB)).norm()
        width = t1.dot(AC)
        
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

            dist = P0.distance(P3)
            vert_dist = P3.depth - P0.depth
            dip = np.degrees(np.arcsin(vert_dist / dist))
            
            d1 = P1.depth * -1
            d2 = P2.depth * -1
            dx = P2.distance(P1)
            dz = d1 - d2
            dip = np.degrees(np.arctan2(dx,dz))
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
            p0 = Vector.fromPoint(P0)
            p1 = Vector.fromPoint(P1)
            p3 = Vector.fromPoint(P3)
            wsum += self.getQuadWidth(p0,p1,p3)
        mwidth = (wsum/len(self.Quadrilaterals))/1000.0
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
        a = (p1-p0).mag()
        b = (p2-p3).mag()
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
        D = np.linalg.det(np.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]]))
        d = -1
        at = np.linalg.det(np.array([[1,y1,z1],[1,y2,z2],[1,y3,z3]]))
        bt = np.linalg.det(np.array([[x1,1,z1],[x2,1,z2],[x3,1,z3]]))
        ct = np.linalg.det(np.array([[x1,y1,1],[x2,y2,1],[x3,y3,1]]))
        a = (-d/D) *at
        b = (-d/D) *bt
        c = (-d/D) *ct

        numer = np.abs(a*otherpoint.x + b*otherpoint.y + c*otherpoint.z + d)
        denom = np.sqrt(a**2 + b**2 + c**2)
        dist = numer/denom
        return dist

    def isPointToRight(self,P0,P1,P2):
        east = np.min([P0.longitude,P1.longitude,P2.longitude])
        west = np.max([P0.longitude,P1.longitude,P2.longitude])
        south = np.min([P0.latitude,P1.latitude,P2.latitude])
        north = np.max([P0.latitude,P1.latitude,P2.latitude])
        proj = utils.get_orthographic_projection(west,east,north,south)
        p0x,p0y = proj(P0.longitude,P0.latitude)
        p1x,p1y = proj(P1.longitude,P1.latitude)
        p2x,p2y = proj(P2.longitude,P2.latitude)
        dx = p1x-p0x
        dy = p1y-p0y
        theta = np.arctan2(dx,dy)
        R = np.array([[np.cos(theta),-np.sin(theta)],
                      [np.sin(theta),np.cos(theta)]])
        xy = np.array([[p2x],[p2y]])
        oxy = np.array([[p0x],[p0y]])
        xp,yp = np.dot(R,xy)
        ox,oy = np.dot(R,oxy)
        if ox > xp:
            return False
        return True

    def reverseQuad(self,P0,P1,P2,P3):
        newP0 = copy.deepcopy(P1)
        newP1 = copy.deepcopy(P0)
        newP2 = copy.deepcopy(P3)
        newP3 = copy.deepcopy(P2)
        if not self.isPointToRight(newP0,newP1,newP2):
            raise FaultException('Third vertex of quadrilateral must be to the right of the second vertex')
        return (newP0,newP1,newP2,newP3)
    
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
        if P1.depth > P2.depth:
            raise FaultException('Top edge of a quadrilateral must be defined by the first two vertices')
        #Is dip angle clockwise and btw 0-90 degrees?
        if not self.isPointToRight(P0,P1,P2):
            P0,P1,P2,P3 = self.reverseQuad(P0,P1,P2,P3)
            print 'Reversing quad where dip not between 0 and 90 degrees.'
        #Are all 4 points (reasonably) co-planar?
        #Translate vertices to ECEF
        p0 = Vector.fromPoint(P0)
        p1 = Vector.fromPoint(P1)
        p2 = Vector.fromPoint(P2)
        p3 = Vector.fromPoint(P3)
        #Calculate normalized vector along top edge
        v0 = (p1-p0).norm()
        #Calculate distance btw p3 and p2
        d = (p3-p2).mag()
        #get the new P2 value
        v1 = v0*d
        newp2 = p3 + v1
        planepoints = [p0,p1,p2]
        dnormal = self.getDistanceToPlane(planepoints,p2)
        geometricMean = self.getTrapMeanLength(p0,p1,newp2,p3)
        if dnormal/geometricMean > OFFPLANE_TOLERANCE:
            raise FaultException('Points in quadrilateral are not co-planar')
        newP0 = p0.toPoint()
        newP1 = p1.toPoint()
        newP2 = newp2.toPoint()
        newP3 = p3.toPoint()
        return [newP0,newP1,newP2,newP3]
        
    
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
        self.lon = np.array(self.lon)
        self.lat = np.array(self.lat)
        self.depth = np.array(self.depth)
        inan = np.isnan(self.lon)
        numnans = len(self.lon[inan])
        numsegments = numnans + 1
        #requirements:
        # 1) Coordinate arrays must be same length
        # 2) Polygons must be quadrilaterals
        # 3) Quads must be closed
        # 4) Quads must be planar
        if len(self.lon) != len(self.lat) != len(self.depth):
            raise IndexError('Length of input lon,lat,depth arrays must be equal')
        
        istart = 0
        endpoints = list(np.where(np.isnan(self.lon))[0])
        endpoints.append(len(self.lon))
        self.Quadrilaterals = []
        for iend in endpoints:
            lonseg = self.lon[istart:iend][0:-1] #remove closing points
            latseg = self.lat[istart:iend][0:-1]
            depthseg = self.depth[istart:iend][0:-1]
            #each segment can have many contiguous quadrilaterals defined in it
            #separations (nans) between segments mean that segments are not contiguous.
            npoints = len(lonseg)
            nquads = ((npoints - 4)/2) + 1
            startidx = 0
            endidx = -1
            for i in range(0,nquads):
                topLeft = point.Point(lonseg[startidx],latseg[startidx],depthseg[startidx])
                topRight = point.Point(lonseg[startidx+1],latseg[startidx+1],depthseg[startidx+1])
                bottomRight = point.Point(lonseg[endidx-1],latseg[endidx-1],depthseg[endidx-1])
                bottomLeft = point.Point(lonseg[endidx],latseg[endidx],depthseg[endidx])
                surface = self.validateQuad(topLeft,topRight,bottomRight,bottomLeft)
                self.Quadrilaterals.append(surface)
                startidx += 1
                endidx -= 1
            istart = iend+1

    def plot(self,ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            if 'xlim3d' not in ax.properties.keys():
                raise FaultException('Non-3d axes object passed to plot() method.')
        for quad in self.Quadrilaterals:
            P0,P1,P2,P3 = quad
            ax.plot([P0.longitude],[P0.latitude],[-P0.depth],'B.')
            ax.text([P0.longitude],[P0.latitude],[-P0.depth],'P0')
            ax.plot([P1.longitude],[P1.latitude],[-P1.depth],'b.')
            ax.text([P1.longitude],[P1.latitude],[-P1.depth],'P1')
            ax.plot([P2.longitude],[P2.latitude],[-P2.depth],'b.')
            ax.text([P2.longitude],[P2.latitude],[-P2.depth],'P2')
            ax.plot([P3.longitude],[P3.latitude],[-P3.depth],'b.')
            ax.text([P3.longitude],[P3.latitude],[-P3.depth],'P3')
            
        
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
        if not len(inan):
            return
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

def _test_northridge():
    #this should fail!
    fault_text = """
    # Source: Wald, D. J., T. H. Heaton, and K. W. Hudnut (1996). The Slip History of the 1994 Northridge, California, Earthquake Determined from Strong-Motion, Teleseismic, GPS, and Leveling Data, Bull. Seism. Soc. Am. 86, S49-S70.
    34.315 -118.421 5.000
    34.401 -118.587 5.000
    34.261 -118.693 20.427
    34.175 -118.527 20.427
    34.315 -118.421 5.000
    """
    cbuf = StringIO.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
    quad = fault.getQuadrilaterals()[0]
    topdist = quad[0].distance(quad[1])
    fault.multiplyFaultLength(2.0)    
    quad = fault.getQuadrilaterals()[0]
    topdist2 = quad[0].distance(quad[1])
    x = 1
    
def _test_correct():
    #this fault should parse correctly
    fault_text = """#SOURCE: Barka, A., H. S. Akyz, E. Altunel, G. Sunal, Z. Akir, A. Dikbas, B. Yerli, R. Armijo, B. Meyer, J. B. d. Chabalier, T. Rockwell, J. R. Dolan, R. Hartleb, T. Dawson, S. Christofferson, A. Tucker, T. Fumal, R. Langridge, H. Stenner, W. Lettis, J. Bachhuber, and W. Page (2002). The Surface Rupture and Slip Distribution of the 17 August 1999 Izmit Earthquake (M 7.4), North Anatolian Fault, Bull. Seism. Soc. Am. 92, 43-60.
    40.70985 29.33760 0
    40.72733 29.51528 0
    40.72933 29.51528 20
    40.71185 29.33760 20
    40.70985 29.33760 0
    >
    40.70513 29.61152 0
    40.74903 29.87519 0
    40.75103 29.87519 20
    40.70713 29.61152 20
    40.70513 29.61152 0
    >
    40.72582 29.88662 0
    40.72336 30.11126 0
    40.73432 30.19265 0
    40.73632 30.19265 20
    40.72536 30.11126 20
    40.72782 29.88662 20
    40.72582 29.88662 0
    >
    40.71210 30.30494 0
    40.71081 30.46540 0
    40.70739 30.56511 0
    40.70939 30.56511 20
    40.71281 30.46540 20
    40.71410 30.30494 20
    40.71210 30.30494 0
    >
    40.71621 30.57658 0
    40.70068 30.63731 0
    40.70268 30.63731 20
    40.71821 30.57658 20
    40.71621 30.57658 0
    >
    40.69947 30.72900 0
    40.79654 30.93655 0
    40.79854 30.93655 20
    40.70147 30.72900 20
    40.69947 30.72900 0
    >
    40.80199 30.94688 0
    40.84501 31.01799 0
    40.84701 31.01799 20
    40.80399 30.94688 20
    40.80199 30.94688 0"""

    cbuf = StringIO.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)

def _test_incorrect():
    fault_text = """# Source: Ji, C., D. V. Helmberger, D. J. Wald, and K.-F. Ma (2003). Slip history and dynamic implications of the 1999 Chi-Chi, Taiwan, earthquake, J. Geophys. Res. 108, 2412, doi:10.1029/2002JB001764.
    24.27980 120.72300	0 
    24.05000 121.00000	17
    24.07190 121.09300	17
    24.33120 121.04300	17
    24.27980 120.72300	0 
    >   
    24.27980 120.72300	0
    23.70000 120.68000	0
    23.60400 120.97200	17
    24.05000 121.00000	17
    24.27980 120.72300	0
    >
    23.60400 120.97200	17 
    23.70000 120.68000	0 
    23.58850 120.58600	0
    23.40240 120.78900	17
    23.60400 120.97200	17"""

    cbuf = StringIO.StringIO(fault_text)
    fault = Fault.readFaultFile(cbuf)
    

if __name__ == '__main__':
    _test_northridge()
    # _test_correct()
    # _test_incorrect()

