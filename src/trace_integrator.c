#include "trace_integrator.h"

#include <assert.h>
#include <stdio.h>

    static double
clip_coord(double x, const double xmax)
{
    while(x < 0.0) {
        x += xmax;
    }
    while(x >= xmax) {
        x -= xmax;
    }
    return x;
}

    params_t *
make_params(const double x0max,
        const double x1max,
        const double scale,
        const size_t ntracers,
        interp2d_t *interp2d_v0,
        interp2d_t *interp2d_v1)
{
    params_t *params = (params_t*)malloc(sizeof(params_t));


    /* FIXME: error checking of the above allocation */

    params->x0max = x0max;
    params->x1max = x1max;
    params->scale = scale;
    params->ntracers = ntracers;
    params->interp2d_v0 = interp2d_v0;
    params->interp2d_v1 = interp2d_v1;

    return params;
}

    void
params_free(params_t *params)
{
    if(params)
        free(params);
}

/*
 * t is unused.
 * y[2*ntracers] -- array of (x0, x1) positions of tracer points.
 * f[2*ntracers] -- will be set to an array of (v0, v1) tracer point velocities.
 * v_params -- a params struct.
 */
    static int
RHS(double t, const double y[], double f[],
        void *v_params)
{
    /* unpack the v_params struct */
    params_t *params = (params_t *)v_params;
    size_t ntracers = params->ntracers;
    double scale = params->scale;
    interp2d_t *interp2d_v0 = params->interp2d_v0;
    interp2d_t *interp2d_v1 = params->interp2d_v1;

    double x0, x1;
    int i;

    for(i=0; i<2*ntracers; i+=2) {
        x0 = clip_coord(y[i], params->x0max);
        x1 = clip_coord(y[i+1], params->x1max);

        f[i]   = scale * interp2d_eval(interp2d_v0, x0, x1);
        f[i+1] = scale * interp2d_eval(interp2d_v1, x0, x1);
    }

    return GSL_SUCCESS;
}

    integrator_t *
make_trace_integrator(const gsl_odeiv_step_type *T,
        const double eps_abs, const double eps_rel,
        const params_t *params)
{

    size_t ndim = 2 * params->ntracers;

    integrator_t *integrator = (integrator_t*)malloc(sizeof(integrator));

    /* FIXME: error checking of the above allocations */

    integrator->s = gsl_odeiv_step_alloc(T, ndim);
    integrator->c = gsl_odeiv_control_y_new(eps_abs, eps_rel);
    integrator->e = gsl_odeiv_evolve_alloc(ndim);

    integrator->sys.function = RHS;
    integrator->sys.jacobian = NULL;
    integrator->sys.dimension = ndim;
    integrator->sys.params = (void*)params;
    
    return integrator;
}

    void
trace_integrator_free(integrator_t *i)
{
    if(i) {
        if(i->e){
            gsl_odeiv_evolve_free(i->e);
            i->e = NULL;
        }
        if(i->c) {
            gsl_odeiv_control_free(i->c);
            i->c = NULL;
        }
        if(i->s) {
            gsl_odeiv_step_free(i->s);
            i->s = NULL;
        }
        free(i);
    }
}

    int
trace_integrator_evolve(integrator_t *i, double *t, double t1, double *h, double y[])
{
    /* XXX: should we do a reset?? */

    int status = GSL_SUCCESS;

    while (*t < t1)
    {
        status = 
            gsl_odeiv_evolve_apply(
                    i->e, i->c, 
                    i->s, &i->sys,
                    t, t1,
                    h, y);

        if (status != GSL_SUCCESS)
            break;
    }
    return status;
}
