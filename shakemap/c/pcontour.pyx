#cython: language_level=3
import numpy as np
cimport numpy as np
cimport cython

cdef extern from "contour.h":
    cdef struct cpoint:
        double x
        double y
        int ti
        int ei
        int gei
        int order
        cpoint *next

    cdef struct cseg:
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        int ci
        int nesting
        int neg_depth
        int pos_depth
        cseg *neg_holes
        cseg *neg_next
        cseg *pos_holes
        cseg *pos_next
        cseg *tmp_hole
        cseg *next
        cpoint *first
        cpoint *last

    cdef struct cinfo:
        double cv
        cseg *closed
        cseg *open

    cdef struct cresult:
        size_t ncarray
        cinfo *contour0
        cinfo *carray

    cresult contour_grid(double *grid, size_t gnx, size_t gny, double gdx,
                         double gdy, double ul_x, double ul_y,
                         double *contour_levels, size_t ncont,
                         int outstyle, int verb)


def pcontour(np.ndarray[double, ndim=2, mode='c']grid, dx, dy, ul_x, ul_y,
             np.ndarray[double, ndim=1, mode='c']contour_levels, outstyle,
             verb=0, fmt=0):
    cdef Py_ssize_t ny = grid.shape[0]
    cdef Py_ssize_t nx = grid.shape[1]
    cdef np.ndarray[double, ndim=1] fg = grid.flatten(order='C')
    cdef cresult cres
    cdef cinfo *crlist
    cdef double cv
    cdef cseg *cp
    cdef cseg *hp
    cdef cpoint *pt
    cdef int i

    gjson = {"type": "FeatureCollection",
             "features": []}

    cres = contour_grid(<double *>fg.data, nx, ny, dx, dy, ul_x, ul_y,
                        <double *>contour_levels.data,
                        np.size(contour_levels),
                        outstyle, verb)

    if outstyle == 3 or outstyle == 4:
        for i in range(-1, <int>cres.ncarray, 1):
            if i == -1:
                crlist = cres.contour0
            else:
                crlist = &(cres.carray[i])
            cv = crlist[0].cv
            if fmt == 0:
                feature = {"type": "Feature",
                           "properties": {"value": cv},
                           "geometry": {"type": "MultiPolygon",
                                        "coordinates": []} }
            else:
                feature = {"type": "Feature",
                           "properties": {"AREA": 0,
                                          "PERIMETER": 0,
                                          "PGAPOL_": i + 1,
                                          "PGAPOL_ID": i + 1,
                                          "GRID_CODE": 0,
                                          "PARAMVALUE": cv},
                           "geometry": {"type": "MultiPolygon",
                                        "coordinates": []} }
            cp = crlist[0].closed
            while cp != NULL:
                clist = []
                plist = []
                pt = cp[0].first
                while pt != NULL:
                    plist.append([round(pt.x, 6), round(pt.y, 6)])
                    pt = pt[0].next
                clist.append(plist)

                hp = cp[0].neg_holes
                while hp != NULL:
                    plist = []
                    pt = hp[0].first
                    while pt != NULL:
                        plist.append([round(pt.x, 6), round(pt.y, 6)])
                        pt = pt[0].next
                    clist.append(plist)
                    hp = hp[0].neg_next

                hp = cp[0].pos_holes
                while hp != NULL:
                    plist = []
                    pt = hp[0].first
                    while pt != NULL:
                        plist.append([round(pt.x, 6), round(pt.y, 6)])
                        pt = pt[0].next
                    clist.append(plist)
                    hp = hp[0].pos_next

                feature['geometry']['coordinates'].append(clist)
                cp = cp[0].next
            if len(feature['geometry']['coordinates']) > 0:
                gjson['features'].append(feature)
    else:
        for i in range(-1, <int>cres.ncarray, 1):
            if i == -1:
                # crlist = cres.contour0
                continue
            else:
                crlist = &(cres.carray[i])
            cv = crlist[0].cv
            feature = {"type": "Feature",
                    "properties": {
                        "value": cv
                        },
                    "geometry": {
                        "type": "MultiLineString",
                        "coordinates": []
                        }
                    }
            clist = []
            cp = crlist[0].closed
            while cp != NULL:
                plist = []
                pt = cp[0].first
                while pt != NULL:
                    plist.append([round(pt.x, 6), round(pt.y, 6)])
                    pt = pt[0].next
                clist.append(plist)
                cp = cp[0].next
            cp = crlist[0].open
            while cp != NULL:
                plist = []
                pt = cp[0].first
                while pt != NULL:
                    plist.append([round(pt.x, 6), round(pt.y, 6)])
                    pt = pt[0].next
                clist.append(plist)
                cp = cp[0].next
            feature['geometry']['coordinates'] = clist
            if len(feature['geometry']['coordinates']) > 0:
                gjson['features'].append(feature)

    return gjson
