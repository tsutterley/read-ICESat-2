#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True

# PURPOSE: create distance metric for windowed classifier
def windowed_manhattan(u, v, window=None):
    """
    Create a windowed manhattan distance metric

    Arguments
    ---------
    u: Input array
    v: Input array for distance

    Keyword arguments
    -----------------
    window: distance window for reducing neighbors
    """
    # calculate manhattan (rectilinear) distances
    cdef double d = 0.0
    d = dist(u, v, window)
    return d
