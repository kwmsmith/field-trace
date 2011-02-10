#ifndef __WRAP_GSL_ODE_H__
#define __WRAP_GSL_ODE_H__

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "gsl_interp2d.h"

typedef struct {
    double x0max, x1max;
    double scale;
    size_t ntracers;
    interp2d_t *interp2d_v0;
    interp2d_t *interp2d_v1;
} params_t;

params_t * make_params(
        const double x0max, const double x1max,
        const double scale, const size_t ntracers,
        interp2d_t *interp2d_v0, interp2d_t *interp2d_v1);

void params_free(params_t *params);

typedef struct {
    gsl_odeiv_step    *s;
    gsl_odeiv_control *c;
    gsl_odeiv_evolve  *e;

    gsl_odeiv_system sys;

} integrator_t;

integrator_t *make_trace_integrator(
        const gsl_odeiv_step_type *T,
        const double eps_abs, const double eps_rel,
        const params_t *params);

void trace_integrator_free(integrator_t *i);

int trace_integrator_evolve(integrator_t *i, double *t, double t1, double *h, double y[]);

#endif
