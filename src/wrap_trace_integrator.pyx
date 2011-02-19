import sys

cdef class TraceIntegrator:

    def __cinit__(self,
            int ntracers, double scale,
            double x0max, double x1max,
            Interp2DPeriodic v0, Interp2DPeriodic v1,
            double eps_abs, double eps_rel):
        
        self.params = \
                make_params(x0max, x1max,
                        scale, <size_t>ntracers,
                        v0.interp2d, v1.interp2d)
        if not self.params:
            raise RuntimeError("failure to make parameter struct")

        # FIXME: we hardcode the integrator here -- how about making a configurable option?
        self.integrator = \
                make_trace_integrator(gsl_odeiv_step_rk8pd,
                        eps_abs, eps_rel,
                        self.params)
        if not self.integrator:
            if self.params:
                params_free(self.params)
            raise RuntimeError("failure to make parameter struct")


    def evolve(self, double t, double t1, double h, np.ndarray[double, ndim=1] y):
        
        retval = trace_integrator_evolve(self.integrator, &t, t1, &h, <double*>y.data)
        if retval:
            raise RuntimeError("error in evole")

        return (t, h)


    def __dealloc__(self):
        if self.params:
            params_free(self.params)

        if self.integrator:
            trace_integrator_free(self.integrator)

def clip(double x, double xmax):
    return clip_coord(x, xmax)
