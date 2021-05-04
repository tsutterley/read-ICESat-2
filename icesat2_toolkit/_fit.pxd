#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

# Inline distance functions

from libc.math cimport fabs, INFINITY

# check if distance is within window
cdef inline double in_window(double d, double w):
    return d if d < w else INFINITY

# calculate windowed manhattan (rectilinear) distances
cdef inline dist(double[:] u, double[:] v, double[:] window):
    cdef int n = u.shape[0]
    cdef double d = 0.0
    for j in range(n):
        d += in_window(fabs(u[j] - v[j]), window[j])
    return d
