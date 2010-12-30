#include "gsl_interp2d.h"

interp2d_t *
interp2d_alloc(const gsl_interp_type *T, size_t nrows, size_t ncols)
{
    int i;
    gsl_spline *tmp_spline = NULL;

    /* the interp2d struct */
    interp2d_t *interp2d = (interp2d_t *) malloc(sizeof(interp2d_t));
    if( NULL == interp2d) {
        return interp2d;
    }

    interp2d->nrows = nrows;
    interp2d->ncols = ncols;

    /* row splines array and each element */
    interp2d->row_splines = (gsl_spline **) malloc(nrows * sizeof(gsl_spline*));
    if(NULL == interp2d->row_splines) {
        goto fail_rowsplines_arr;
    } else {
        for(i=0; i<nrows; i++)
            interp2d->row_splines[i] = NULL;
    }

    for(i=0; i<nrows; i++) {
        if(NULL == (interp2d->row_splines[i] = gsl_spline_alloc(T, ncols))) {
            goto fail_rowsplines;
    }

    /* row_spline_accels array and each element */
    interp2d->row_spline_accels = (gsl_interp_accel **) malloc(nrows *sizeof(gsl_interp_accel*));
    if(NULL == interp2d->row_spline_accels) {
        goto fail_rowsplines;
    } else {
        for(i=0; i<nrows; i++) {
            interp2d->row_spline_accels[i] = NULL;
        }
    }
    for(i=0; i<nrows; i++) {
        if(NULL == (interp2d->row_spline_accels[i] = gsl_interp_accel_alloc())) {
            goto fail_accels;
        }
    }

    /* col_spline and col_spline_accel */
    interp2d->col_spline = gsl_spline_alloc(T, gsl_spline_min_size(interp2d->row_splines[0]));
    if(NULL == interp2d->col_spline) {
        goto fail_accels;
    }

    if(NULL == (interp2d->col_spline_accel = gsl_interp_accel_alloc())) {
        goto fail_col_spline;
    }

    goto success;

fail_col_spline:
    free(interp2d->col_spline);
    interp2d->col_spline = NULL;
fail_accels:
    for(i=0; i<nrows; i++) {
        if(interp2d->row_spline_accels[i]) {
            gsl_interp_accel_free(interp2d->row_spline_accels[i]);
            interp2d->row_spline_accels[i] = NULL;
        }
    }
    free(interp2d->row_spline_accels);
    interp2d->row_spline_accels = NULL;
fail_rowsplines:
    for(i=0; i<nrows; i++) {
        if(interp2d->row_splines[i]) {
            gsl_spline_free(interp2d->row_splines[i]);
            interp2d->row_splines[i] = NULL;
        }
    }
    free(interp2d->row_splines);
    interp2d->row_splines = NULL;
fail_rowsplines_arr:
    free(interp2d);
    interp2d = NULL;
success:
    return interp2d;
}

int
interp2d_init(interp2d_t *interp2d, const double *xa, const double *ya, size_t nrows, size_t ncols)
{
}

double
interp2d_eval(const interp2d_t *interp2d, double x, double y)
{
}

void
interp2d_free(interp2d_t *interp2d)
{
    int i;
    size_t nrows = interp2d->nrows;
    size_t ncols = interp2d->ncols;

    if(interp2d) {
        if(interp2d->col_spline_accel) {
            free(interp2d->col_spline_accel);
            interp2d->col_spline_accel = NULL;
        }
        if(interp2d->col_spline) {
            free(interp2d->col_spline);
            interp2d->col_spline = NULL;
        }
        if(interp2d->row_spline_accels) {
            for(i=0; i<nrows; i++) {
                if(interp2d->row_spline_accels[i]) {
                    gsl_interp_accel_free(interp2d->row_spline_accels[i]);
                    interp2d->row_spline_accels[i] = NULL;
                }
            }
            free(interp2d->row_spline_accels);
            interp2d->row_spline_accels = NULL;
        }
        if(interp2d->row_splines) {
            for(i=0; i<nrows; i++) {
                if(interp2d->row_splines[i]) {
                    gsl_spline_free(interp2d->row_splines[i]);
                    interp2d->row_splines[i] = NULL;
                }
            }
            free(interp2d->row_splines);
            interp2d->row_splines = NULL;
        }
        free(interp2d);
        interp2d = NULL;
    }
}
