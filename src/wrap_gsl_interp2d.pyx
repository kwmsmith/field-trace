cdef class Interp2DPeriodic:

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
