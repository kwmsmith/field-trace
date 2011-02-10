cdef extern from "gsl_interp2d.h":

    ctypedef struct interp2d_t:
        pass

    interp2d_t *interp2d_make_periodic( double x0max, double x1max,
            double *arr, size_t nrows, size_t ncols)
    double interp2d_eval(interp2d_t *interp2d, double x, double y)
    void interp2d_free(interp2d_t *interp2d)

import numpy as np
cimport numpy as np

cdef class Interp2DPeriodic:

    cdef interp2d_t *interp2d

    def __cinit__(self,
                  double x0max, double x1max,
                  np.ndarray[double, ndim=2] arr):
        self.interp2d = interp2d_make_periodic(x0max, x1max, <double*>arr.data, arr.shape[0], arr.shape[1])
        if not self.interp2d:
            raise RuntimeError("failed to make an interp2d object")

    cpdef double eval(self, double x, double y):
        return interp2d_eval(self.interp2d, x, y)

    def __dealloc__(self):
        if self.interp2d:
            interp2d_free(self.interp2d)
