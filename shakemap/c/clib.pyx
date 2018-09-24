import numpy as np
cimport cython
from cython.parallel import prange
from libc.math cimport (sqrt,
                        cos,
                        exp)

cdef double EARTH_RADIUS = 6371.


@cython.boundscheck(False)
@cython.wraparound(False)
def make_sigma_matrix(double[:, ::1]corr12, double[:, ::1]corr_adj12,
                      double[:]sdsta, double[:]sdarr):
    cdef Py_ssize_t ny = corr12.shape[0]
    cdef Py_ssize_t nx = corr12.shape[1]

    cdef Py_ssize_t x, y

    for y in prange(ny, nogil=True, schedule=dynamic):
        for x in range(nx):
            corr12[y, x] = corr12[y, x] * corr_adj12[y, x] * \
                           sdsta[y] * sdarr[x]
    return


@cython.boundscheck(False)
@cython.wraparound(False)
def geodetic_distance_fast_c(double[::1]lons1, double[::1]lats1,
                             double[::1]lons2, double[::1]lats2,
                             double[:, ::1]result):
    cdef Py_ssize_t nx = lons1.shape[0]
    cdef Py_ssize_t ny = lons2.shape[0]

    cdef Py_ssize_t x, y

    if &lons1[0] == &lons2[0] and &lats1[0] == &lats2[0]:
        for y in prange(ny, nogil=True, schedule='guided'):
            for x in range(y+1):
                result[y, x] = result[x, y] = (
                    EARTH_RADIUS *
                    sqrt(((lons1[x] - lons2[y]) *
                        cos(0.5 * (lats1[x] + lats2[y])))**2 +
                        (lats1[x] - lats2[y])**2))
    else:
        for y in prange(ny, nogil=True, schedule=dynamic):
            for x in range(nx):
                result[y, x] = (
                    EARTH_RADIUS *
                    sqrt(((lons1[x] - lons2[y]) *
                        cos(0.5 * (lats1[x] + lats2[y])))**2 +
                        (lats1[x] - lats2[y])**2))
    return


@cython.boundscheck(False)
@cython.wraparound(False)
def eval_lb_correlation(double[:, ::1]b1, double[:, ::1]b2, double[:, ::1]b3,
                        long[:, ::1]ix1, long[:, ::1]ix2, double[:, ::1]h):
    cdef Py_ssize_t nx = ix1.shape[1]
    cdef Py_ssize_t ny = ix1.shape[0]

    cdef Py_ssize_t x, y, i, j
    cdef double hval
    cdef double afact = -3.0 / 20.0
    cdef double bfact = -3.0 / 70.0

    for y in prange(ny, nogil=True, schedule=dynamic):
        for x in range(nx):
            hval = h[y, x]
            i = ix1[y, x]
            j = ix2[y, x]
            h[y, x] = (b1[i, j] * exp(hval * afact) +
                       b2[i, j] * exp(hval * bfact))
            if hval == 0:
                h[y, x] += b3[i, j]

    return h


@cython.boundscheck(False)
@cython.wraparound(False)
def make_sd_array(double[:, ::1]sdgrid, double[:, ::1]pout_sd2, long iy,
                  double[:, ::1]rcmatrix, double[::1, :]sigma12):
    cdef Py_ssize_t nx = rcmatrix.shape[1]
    cdef Py_ssize_t ny = rcmatrix.shape[0]

    tmp = np.zeros(ny, dtype=np.double)
    cdef double[::1] tmp_view = tmp

    cdef Py_ssize_t x, y

    for y in prange(ny, nogil=True):
        for x in range(nx):
            tmp_view[y] += rcmatrix[y, x] * sigma12[y, x]
        sdgrid[iy, y] = pout_sd2[iy, y] - tmp_view[y]
        if sdgrid[iy, y] < 0:
            sdgrid[iy, y] = 0
        sdgrid[iy, y] = sqrt(sdgrid[iy, y])
    return
