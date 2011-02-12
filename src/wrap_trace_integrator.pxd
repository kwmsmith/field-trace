cdef extern from "gsl/gsl_odeiv.h":

    ctypedef struct gsl_odeiv_step_type:
        pass

    gsl_odeiv_step_type *gsl_odeiv_step_rk2
    gsl_odeiv_step_type *gsl_odeiv_step_rk4
    gsl_odeiv_step_type *gsl_odeiv_step_rkf45
    gsl_odeiv_step_type *gsl_odeiv_step_rkck
    gsl_odeiv_step_type *gsl_odeiv_step_rk8pd

cdef extern from "gsl_interp2d.h":

    ctypedef struct interp2d_t:
        pass

cdef extern from "trace_integrator.h":

    ctypedef struct integrator_t:
        pass

    ctypedef struct params_t:
        pass

    params_t *make_params(
            double x0max, double x1max,
            double scale, size_t ntracers,
            interp2d_t *interp2d_v0, interp2d_t *interp2d_v1)

    void params_free(params_t *params)

    integrator_t *make_trace_integrator(
            gsl_odeiv_step_type *T,
            double eps_abs, double eps_rel,
            params_t *params)

    void trace_integrator_free(integrator_t *i)

    int trace_integrator_evolve(integrator_t *i, double *t, double t1, double *h, double *y)


from wrap_gsl_interp2d cimport Interp2DPeriodic
cimport numpy as np

cdef class TraceIntegrator:

    cdef params_t *params
    cdef integrator_t *integrator
