import sys

cdef class TraceIntegrator:

    def __cinit__(self,
            int ntracers, double scale,
            double x0max, double x1max,
            Interp2DPeriodic v0, Interp2DPeriodic v1,
            double eps_abs, double eps_rel):

        sys.stderr.write("here 1\n")

        self.params = \
                make_params(x0max, x1max,
                        scale, <size_t>ntracers,
                        v0.interp2d, v1.interp2d)
        if not self.params:
            raise RuntimeError("failure to make parameter struct")

        sys.stderr.write("here 2\n")

        # FIXME: we hardcode the integrator here -- how about making a configurable option?
        self.integrator = \
                make_trace_integrator(gsl_odeiv_step_rk8pd,
                        eps_abs, eps_rel,
                        self.params)
        if not self.integrator:
            if self.params:
                params_free(self.params)
            raise RuntimeError("failure to make parameter struct")

        sys.stderr.write("here 3\n")

    def __dealloc__(self):
        sys.stderr.write("__dealloc__ 1\n")
        if self.params:
            params_free(self.params)
        sys.stderr.write("__dealloc__ 2\n")

        if self.integrator:
            trace_integrator_free(self.integrator)
        sys.stderr.write("__dealloc__ 3\n")
