#!/usr/bin/env python

#stdlib imports
import os.path
import sys

#third party
import numpy as np

#hack the path so that I can debug these functions if I need to
homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
shakedir = os.path.abspath(os.path.join(homedir,'..'))
sys.path.insert(0,shakedir) #put this at the front of the system path, ignoring any installed mapio stuff

from shakemap.grind.distance import *
from shakemap.grind.vector import Vector

if __name__ == '__main__':
    flat = [28.427,27.986,27.469,27.956]
    flon = [84.614,86.179,85.880,84.461]
    fdep = [20.0,20.0,13.0,13.0]
    lat0,lat1,lat2,lat3 = flat
    lon0,lon1,lon2,lon3 = flon
    dep0,dep1,dep2,dep3 = fdep
    
    P0 = Vector.fromPoint(point.Point(lon0,lat0,dep0))
    P1 = Vector.fromPoint(point.Point(lon1,lat1,dep1))
    P2 = Vector.fromPoint(point.Point(lon2,lat2,dep2))
    P3 = Vector.fromPoint(point.Point(lon3,lat3,dep3))

    lons = np.arange(81.5,87.5,0.0083)
    lats = np.arange(25.5,31.0,0.0083)
    lonmat,latmat = np.meshgrid(lons,lats)
    depmat = np.zeros_like(lonmat)

    x,y,z = latlon2ecef(lonmat,latmat,depmat)
    nr,nc = x.shape
    x.shape = (nr*nc,1)
    y.shape = (nr*nc,1)
    z.shape = (nr*nc,1)
    points = np.hstack((x,y,z))
    t1 = datetime.now()
    dist = calcRuptureDistance(P0,P1,P2,P3,points)
    t2 = datetime.now()
    dt = t2-t1
    sec = dt.seconds + float(dt.microseconds)/1e6
    npoints = nr*nc
    print('RRupt: %.4f seconds for %i point grid' % (sec,npoints))
    dist.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.imshow(dist,interpolation='nearest')
        plt.colorbar()
        plt.savefig('dist.png')
    else:
        print(dist)

    #test this against Bruce's C program
    nr2,nc2 = points.shape
    f = open('points.bin','wb')
    for i in range(0,nr2):
        xt,yt,zt = points[i,:]
        xb = struct.pack('d',xt)
        yb = struct.pack('d',yt)
        zb = struct.pack('d',zt)
        f.write(xb)
        f.write(yb)
        f.write(zb)
    f.close()

    #now write out the quads
    f = open('fault.txt','wt')
    f.write('5\n')
    f.write('%.11f %.11f %.11f\n' % tuple(P0.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P1.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P2.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P3.getArray()))
    f.write('%.11f %.11f %.11f\n' % tuple(P0.getArray()))
    f.close()

    cmd = '/Users/mhearne/src/python/ShakeMap/src/contour/dist_rrupt fault.txt < points.bin > outpoints.bin'
    res,stdout,stderr = getCommandOutput(cmd)
    f = open('outpoints.bin','rb')
    dlist = []
    buf = f.read(8)
    while buf:
        d = struct.unpack('d',buf)
        dlist.append(d[0])
        buf = f.read(8)
    f.close()
    d2 = np.array(dlist)
    d2.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(d2,interpolation='nearest')
        plt.colorbar()
        plt.savefig('dist_c.png')
    else:
        print(dist)
    
    P0,P1 = getTopEdge(flat,flon,fdep)
    rxdist = calcRxDistance(P0,P1,points)
    rxdist.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(rxdist,interpolation='nearest')
        plt.colorbar()
        plt.savefig('distrx.png')
    else:
        print(rxdist)

    #test against Bruce's C implementation
    f = open('points.txt','wt')
    for i in range(0,nr2):
        xt,yt = (points[i,0],points[i,1])
        f.write('%.11f %.11f\n' %(xt,yt))
    f.close()
    cmd = '/Users/mhearne/src/python/ShakeMap/src/contour/dist_rx %.11f %.11f %.11f %.11f < points.txt > outpoints.txt' % (P0.x,P0.y,P1.x,P1.y)
    res,stdout,stderr = getCommandOutput(cmd)
    rxdist2 = []
    for line in open('outpoints.txt','rt').readlines():
        rxdist2.append(float(line.strip()))
    rxdist2 = np.array(rxdist2)
    rxdist2.shape = (nr,nc)
    if nr > 1 and nc > 1:
        plt.close()
        plt.imshow(rxdist2,interpolation='nearest')
        plt.colorbar()
        plt.savefig('distrx_c.png')
    else:
        print(rxdist2)
