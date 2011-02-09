cdef extern from "gsl/gsl_interp.h":
    ctypedef struct gsl_interp_type:
        pass

cdef extern from "gsl/gsl_spline.h":
    gsl_interp_type *gsl_interp_cspline
    gsl_interp_type *gsl_interp_cspline_periodic

cdef extern from "gsl_interp2d.h":

    ctypedef struct interp2d_t:
        pass

    interp2d_t *interp2d_alloc(gsl_interp_type *row_type,
                               gsl_interp_type *col_type, size_t nrows, size_t ncols)
    int interp2d_init(interp2d_t *interp2d, double *x0, double max_x0,
                      double *x1, double max_x1, double *ya)
    double interp2d_eval(interp2d_t *interp2d, double x, double y)
    void interp2d_free(interp2d_t *interp2d)

import numpy as np
cimport numpy as np

cdef class Interp2DPeriodic:

    cdef interp2d_t *interp2d

    def __cinit__(self,
                  double x0max, double x1max,
                  np.ndarray[double, ndim=2] arr):

        cdef size_t nrows, ncols

        # 1D interpolation in gsl requires periodic interpolation to have the
        # first and last values of the X & Y arrays equal to each other.  The
        # 2D interpolation here is periodic in the 2nd dimension, so we need to
        # do a copy of `arr` into a larger array `arr_periodic_y`, and copy the
        # first column to the new last column.
        cdef np.ndarray[double, ndim=2] arr_periodic_y = \
                np.empty((arr.shape[0], arr.shape[1]+1), dtype=np.double)

        arr_periodic_y[:, :-1] = arr
        arr_periodic_y[:,-1] = arr[:,0]

        # Set-up the coordinate indices for the 0th and 1st dimensions -- x0
        # and x1 -- minding the comment above about gsl's requirements for
        # periodicity.
        cdef np.ndarray[double, ndim=1] x0, x1
        x0 = np.linspace(0.0, x0max, arr_periodic_y.shape[0], endpoint=False)
        x1 = np.arange(arr_periodic_y.shape[1], dtype=np.double)
        x1_idx_max = np.max(x1)
        x1 = (x1max / x1_idx_max) * x1

        nrows = x0.shape[0]
        ncols = x1.shape[0]

        self.interp2d = interp2d_alloc(gsl_interp_cspline_periodic,
                                       gsl_interp_cspline,
                                       nrows, ncols)
        # interp2d_alloc returns NULL on error.
        if not self.interp2d:
            raise RuntimeError("error allocating interp2d_t struct")

        # interp2d_init returns nonzeror on error.
        if interp2d_init(self.interp2d, <double*>x0.data, x0max,
                <double*>x1.data, x1max, <double*>arr_periodic_y.data):
            interp2d_free(self.interp2d)
            raise RuntimeError("error initializing interp2d_t struct")

    cpdef double eval(self, double x, double y):
        return interp2d_eval(self.interp2d, x, y)

    def __dealloc__(self):
        if self.interp2d:
            interp2d_free(self.interp2d)
