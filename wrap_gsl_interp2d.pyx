cdef extern from "gsl/gsl_interp.h":
    ctypedef struct gsl_interp_type:
        pass

cdef extern from "gsl/gsl_spline.h":
    gsl_interp_type *gsl_interp_cspline
    gsl_interp_type *gsl_interp_cspline_periodic

cdef extern from "gsl_interp2d.h":

    ctypedef struct interp2d_t:
        pass

    interp2d_t *interp2d_alloc(gsl_interp_type *T, size_t nrows, size_t ncols)
    int interp2d_init(interp2d_t *interp2d, double *x0, double *x1, double *ya, size_t nrows, size_t ncols)
    double interp2d_eval(interp2d_t *interp2d, double x, double y)
    void interp2d_free(interp2d_t *interp2d)

import numpy as np
cimport numpy as np

cdef class Interp2DPeriodic:

    cdef interp2d_t *interp2d
    cdef size_t nrows, ncols

    def __cinit__(self,
                  np.ndarray[double, ndim=1] x0,
                  np.ndarray[double, ndim=1] x1,
                  np.ndarray[double, ndim=2] arr):
        import sys
        self.nrows = x0.shape[0]
        self.ncols = x1.shape[0]
        self.interp2d = interp2d_alloc(gsl_interp_cspline_periodic,
                                       self.nrows, self.ncols)
        if not self.interp2d:
            raise RuntimeError("error allocating interp2d_t struct")

        if interp2d_init(self.interp2d, <double*>x0.data, <double*>x1.data, <double*>arr.data, self.nrows, self.ncols):
            # interp2d_free(self.interp2d)
            raise RuntimeError("error initializing interp2d_t struct")

    cpdef double eval(self, double x, double y):
        return interp2d_eval(self.interp2d, x, y)

    def __dealloc__(self):
        if self.interp2d:
            interp2d_free(self.interp2d)
