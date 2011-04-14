cimport cython

import numpy as np
cimport numpy as np

cdef inline int sign(double a):
    if a < 0.0: return -1
    elif a > 0.0: return 1
    else: return 0

def same_sign_or_zero(double a0, double a1, double a2, double a3, double lval=0.0):
    cdef double arr[4]
    arr[0] = a0-lval; arr[1] = a1-lval; arr[2] = a2-lval; arr[3] = a3-lval
    cdef int i
    cdef int sign0, signi

    for i in range(4):
        signi = sign(arr[i])
        if signi != 0:
            sign0 = signi

    for i in range(4):
        signi = sign(arr[i])
        if signi and sign0 != signi:
            return 0

    return 1

def same_sign(double a0, double a1, double a2, double a3):
    cdef double arr[4]
    arr[0] = a0; arr[1] = a1; arr[2] = a2; arr[3] = a3
    cdef int sign0 = sign(arr[0])
    cdef int i
    for i in range(4):
        signi = sign(arr[i])
        if signi != sign0:
            return 0
    return 1

cpdef int clip(int i, int nx):
    while i < 0:
        i += nx
    while i >= nx:
        i -= nx
    return i

def neighbors8(int i, int j, size_t nx, size_t ny):
    return (
            (clip(i-1, nx), clip(j-1, ny)),
            (clip(i+1, nx), clip(j-1, ny)),
            (clip(i+1, nx), clip(j+1, ny)),
            (clip(i-1, nx), clip(j+1, ny)),

            (clip(i+1, nx), j),
            (i, clip(j+1, ny)),
            (clip(i-1, nx), j),
            (i, clip(j-1, ny)))


def neighbors4(int i, int j, size_t nx, size_t ny):
    # return (((i+1)%nx, j),
            # (i, (j+1)%ny),
            # ((i-1)%nx, j),
            # (i, (j-1)%ny))
    return ((clip(i+1, nx), j),
            (i, clip(j+1, ny)),
            (clip(i-1, nx), j),
            (i, clip(j-1, ny)))


@cython.boundscheck(False)
@cython.wraparound(False)
def flood_fill(np.ndarray[int, ndim=2, mode='c'] arr, int i, int j, int border_color, int fill_color):
    cdef size_t nx, ny
    nx = arr.shape[0]; ny = arr.shape[1]
    i %= nx; j %= ny
    if arr[i, j] == border_color:
        return
    edge = [(i,j)]
    arr[i,j] = fill_color
    while edge:
        newedge = []
        for (i, j) in edge:
            for (s, t) in neighbors4(i, j, nx, ny):
                if arr[s,t] not in (border_color, fill_color):
                    arr[s,t] = fill_color
                    newedge.append((s,t))
        edge = newedge

@cython.boundscheck(False)
@cython.wraparound(False)
def _get_corner_vals(np.ndarray[double, ndim=2, mode='c'] f, int i, int j):

    # #  y --------->
    # #  
    # #  f00 ---- f01  x
    # #  |         |   |
    # #  |         |   |
    # #  f10 ---- f11  V

    cdef int nx = f.shape[0], ny = f.shape[1]
    cdef int istep = clip(i+1, nx)
    cdef int jstep = clip(j+1, ny)

    return (f[i, j], f[i, jstep], f[istep, j], f[istep, jstep])
