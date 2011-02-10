cdef extern from "gsl_interp2d.h":

    ctypedef struct interp2d_t:
        pass

    interp2d_t *interp2d_make_periodic( double x0max, double x1max,
            double *arr, size_t nrows, size_t ncols)
    double interp2d_eval(interp2d_t *interp2d, double x, double y)
    void interp2d_free(interp2d_t *interp2d)

cimport numpy as np

cdef class Interp2DPeriodic:

    cdef interp2d_t *interp2d

    cpdef double eval(self, double x, double y)
