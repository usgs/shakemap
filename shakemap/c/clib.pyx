#cython: language_level=3
import numpy as np
cimport cython
from cython.parallel import prange
from libc.math cimport (sqrt,
                        cos,
                        sin,
                        asin,
                        exp)

@cython.boundscheck(False)
@cython.wraparound(False)
def make_sigma_matrix(double[:, ::1]corr12, double[:]sdsta, double[:]sdarr):
    cdef Py_ssize_t ny = corr12.shape[0]
    cdef Py_ssize_t nx = corr12.shape[1]

    cdef double *c12p
    cdef double sdval
    cdef double tmp
    cdef Py_ssize_t x, y

    for y in prange(ny, nogil=True, schedule=dynamic):
        c12p = &corr12[y, 0]
        sdval = sdarr[y]
        for x in range(nx):
            # Putting these operations all on one line seems to
            # allow the compiler to do things that result in the
            # output matrix being very slightly asymmetric.
            tmp = sdsta[x] * sdval
            c12p[x] = c12p[x] * tmp
    return


@cython.boundscheck(False)
@cython.wraparound(False)
def geodetic_distance_fast(double[::1]lons1, double[::1]lats1,
                           double[::1]lons2, double[::1]lats2,
                           double[:, ::1]result):
    cdef double EARTH_RADIUS = 6371.
    cdef Py_ssize_t nx = lons1.shape[0]
    cdef Py_ssize_t ny = lons2.shape[0]

    cdef double lon2, lat2
    cdef double *res
    cdef Py_ssize_t x, y

    if &lons1[0] == &lons2[0] and &lats1[0] == &lats2[0]:
        for y in prange(ny, nogil=True, schedule='guided'):
            lon2 = lons2[y]
            lat2 = lats2[y]
            for x in range(y+1):
                result[y, x] = result[x, y] = (
                    EARTH_RADIUS *
                    sqrt(((lons1[x] - lon2) *
                        cos(0.5 * (lats1[x] + lat2)))**2 +
                        (lats1[x] - lat2)**2))
    else:
        for y in prange(ny, nogil=True, schedule=dynamic):
            res = &result[y, 0]
            lon2 = lons2[y]
            lat2 = lats2[y]
            for x in range(nx):
                res[x] = (
                    EARTH_RADIUS *
                    sqrt(((lons1[x] - lon2) *
                        cos(0.5 * (lats1[x] + lat2)))**2 +
                        (lats1[x] - lat2)**2))
    return


@cython.boundscheck(False)
@cython.wraparound(False)
def geodetic_distance_haversine(double[::1]lons1, double[::1]lats1,
                                double[::1]lons2, double[::1]lats2,
                                double[:, ::1]result):
    cdef double EARTH_RADIUS = 6371.
    cdef Py_ssize_t nx = lons1.shape[0]
    cdef Py_ssize_t ny = lons2.shape[0]

    cdef Py_ssize_t x, y
    cdef double diameter = 2.0 * EARTH_RADIUS

    if &lons1[0] == &lons2[0] and &lats1[0] == &lats2[0]:
        for y in prange(ny, nogil=True, schedule='guided'):
            for x in range(y+1):
                result[y, x] = result[x, y] = (
                    diameter * asin(sqrt(
                        sin((lats1[x] - lats2[y]) / 2.0)**2 +
                        cos(lats1[x]) * cos(lats2[y]) *
                        sin((lons1[x] - lons2[y]) / 2.0)**2)))
    else:
        for y in prange(ny, nogil=True, schedule=dynamic):
            for x in range(nx):
                result[y, x] = (
                    diameter * asin(sqrt(
                        sin((lats1[x] - lats2[y]) / 2.0)**2 +
                        cos(lats1[x]) * cos(lats2[y]) *
                        sin((lons1[x] - lons2[y]) / 2.0)**2)))
    return


@cython.boundscheck(False)
@cython.wraparound(False)
def eval_lb_correlation(double[:, ::1]b1, double[:, ::1]b2, double[:, ::1]b3,
                        long[:, ::1]ix1, long[:, ::1]ix2, double[:, ::1]h):
    cdef Py_ssize_t nx = ix1.shape[1]
    cdef Py_ssize_t ny = ix1.shape[0]

    cdef Py_ssize_t x, y, i, j
    cdef double hval
    cdef long *ix1p
    cdef long *ix2p
    cdef double *hp
    cdef double afact = -3.0 / 20.0
    cdef double bfact = -3.0 / 70.0

    for y in prange(ny, nogil=True, schedule=dynamic):
        hp = &h[y, 0]
        ix1p = &ix1[y, 0]
        ix2p = &ix2[y, 0]
        for x in range(nx):
            hval = hp[x]
            i = ix1p[x]
            j = ix2p[x]
            hp[x] = (b1[i, j] * exp(hval * afact) +
                     b2[i, j] * exp(hval * bfact))
            if hval == 0:
                hp[x] += b3[i, j]

    return h


@cython.boundscheck(False)
@cython.wraparound(False)
def make_sd_array(double[:, ::1]sdgrid, double[:, ::1]pout_sd2, long iy,
                  double[:, ::1]rcmatrix, double[:, ::1]sigma12):
    cdef Py_ssize_t nx = rcmatrix.shape[1]
    cdef Py_ssize_t ny = rcmatrix.shape[0]

    cdef double tmp
    cdef double *sdg = &sdgrid[iy, 0]
    cdef double *pop = &pout_sd2[iy, 0]
    cdef double *rcp
    cdef double *sgp
    cdef Py_ssize_t x, y

    for y in prange(ny, nogil=True):
        rcp = &rcmatrix[y, 0]
        sgp = &sigma12[y, 0]
        tmp = 0
        for x in range(nx):
            tmp = tmp + rcp[x] * sgp[x]
        sdg[y] = pop[y] - tmp
        if sdg[y] < 0:
            sdg[y] = 0
        # sdg[y] = sqrt(sdg[y])
    return
