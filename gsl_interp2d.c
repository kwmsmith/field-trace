#include "gsl_interp2d.h"

#include <stdio.h>
#include <math.h>

#include <assert.h>

#define COL_SPLINE_LEN 5

    interp2d_t *
interp2d_alloc(const gsl_interp_type *row_type, const gsl_interp_type *col_type, size_t nrows, size_t ncols)
{
    int i;

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
        if(NULL == (interp2d->row_splines[i] = gsl_spline_alloc(row_type, ncols))) {
            goto fail_rowsplines;
        }
    }

    /* row_spline_accel */
    if(NULL == (interp2d->row_spline_accel = gsl_interp_accel_alloc())) {
        goto fail_accel;
    }

    assert(COL_SPLINE_LEN < interp2d->nrows);

    /* col_spline and col_spline_accel */
    interp2d->col_spline_len = COL_SPLINE_LEN;
    interp2d->col_spline = gsl_spline_alloc(col_type, interp2d->col_spline_len);
    if(NULL == interp2d->col_spline) {
        goto fail_col_spline;
    }

    if(NULL == (interp2d->col_spline_accel = gsl_interp_accel_alloc())) {
        goto fail_col_spline;
    }

    interp2d->global_x0 = (double*)malloc(interp2d->nrows*sizeof(double));
    interp2d->col_spline_x = (double*)malloc(interp2d->col_spline_len*sizeof(double));
    interp2d->col_spline_y = (double*)malloc(interp2d->col_spline_len*sizeof(double));

    goto success;

fail_col_spline:
    free(interp2d->col_spline);
    interp2d->col_spline = NULL;
fail_accel:
    if(interp2d->row_spline_accel)
        gsl_interp_accel_free(interp2d->row_spline_accel);
    interp2d->row_spline_accel = NULL;
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

/* NOTE: x0[0..nrows-1] is an array of the coordinate values along the 0th dimension.
 *       x1[0..ncols-1] is an array of the coordinate values along the 1st dimension.
 *       ya[0..nrows, 0..ncols] is a 2D array of Y values, of shape NROWS x NCOLS.
 */
    int
interp2d_init(
        interp2d_t *interp2d,
        const double *x0,
        const double max_x0,
        const double *x1,
        const double max_x1,
        const double *ya)
{
    int retval;
    size_t i;
    interp2d->max_x0 = max_x0;
    interp2d->max_x1 = max_x1;
    for(i=0; i<interp2d->nrows; i++) {
        interp2d->global_x0[i] = x0[i];
        if((retval = gsl_spline_init(interp2d->row_splines[i], x1, &(ya[interp2d->ncols*i]), interp2d->ncols))) {
            return retval;
        }
    }
    return 0;
}

    double
interp2d_eval(const interp2d_t *interp2d, double x, double y)
{
    double retval;
    long x0_idx, lower_idx, upper_idx, i, idx;

    /* NOTE: should the interp2d->nrows be instead `interp2d->nrows-1` ??
     */
    x0_idx = gsl_interp_bsearch(interp2d->global_x0, x, 0, interp2d->nrows);

    /* for now we only cover odd-length splines.
     */
    assert(interp2d->col_spline_len % 2 == 1);

    lower_idx = x0_idx - (long)interp2d->col_spline_len / 2;
    upper_idx = x0_idx + (long)interp2d->col_spline_len / 2;

    /* initialize col_spline_x
     */
    for(i=0, idx=lower_idx; idx<=upper_idx; i++, idx++) {
        if(idx < 0) {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx + interp2d->nrows] - interp2d->max_x0;
        } else if(idx >= interp2d->nrows) {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx - interp2d->nrows] + interp2d->max_x0;
        } else {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx];
        }
    }

    /* set col_spline_y array.
     */
    lower_idx = x0_idx - (long)interp2d->col_spline_len / 2 + interp2d->nrows;
    upper_idx = x0_idx + (long)interp2d->col_spline_len / 2 + interp2d->nrows;
    for(i=0, idx=lower_idx; idx<=upper_idx; i++, idx++) {
        interp2d->col_spline_y[i] = 
            gsl_spline_eval(
                    interp2d->row_splines[idx % interp2d->nrows],
                    y, interp2d->row_spline_accel);
    }

    /* FIXME: no error checking here */
    gsl_spline_init(
            interp2d->col_spline, 
            interp2d->col_spline_x, 
            interp2d->col_spline_y, 
            interp2d->col_spline_len);

    gsl_interp_accel_reset(interp2d->col_spline_accel);
    retval = gsl_spline_eval(interp2d->col_spline, x, interp2d->col_spline_accel);
    return retval;
}

    void
interp2d_free(interp2d_t *interp2d)
{
    int i;
    size_t nrows = interp2d->nrows;

    if(interp2d) {
        if(interp2d->global_x0) {
            free(interp2d->global_x0);
            interp2d->global_x0 = NULL;
        }
        if(interp2d->col_spline_x) {
            free(interp2d->col_spline_x);
            interp2d->col_spline_x = NULL;
        }
        if(interp2d->col_spline_y) {
            free(interp2d->col_spline_y);
            interp2d->col_spline_y = NULL;
        }
        if(interp2d->col_spline_accel) {
            free(interp2d->col_spline_accel);
            interp2d->col_spline_accel = NULL;
        }
        if(interp2d->col_spline) {
            free(interp2d->col_spline);
            interp2d->col_spline = NULL;
        }
        if(interp2d->row_spline_accel) {
            gsl_interp_accel_free(interp2d->row_spline_accel);
            interp2d->row_spline_accel = NULL;
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
