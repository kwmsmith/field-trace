#include "gsl_interp2d.h"

#include <stdio.h>
#include <math.h>
#include <assert.h>

#define COL_SPLINE_LEN 5

#define IDX(i,j,ncols) ((j)+(i)*(ncols))

    interp2d_t *
interp2d_make_periodic(const double x0max, const double x1max, const double *arr, const size_t nrows, const size_t ncols)
{
    interp2d_t *interp2d = NULL;
    int i,j;
    size_t ncols_periodic = ncols+1;
    double *x0 = NULL, *x1 = NULL;
    double *arr_periodic_y = (double*)malloc(nrows * ncols_periodic * sizeof(double));
    if(NULL == arr_periodic_y) {
        return NULL;
    }

    /* Initialize arr_periodic_y array.  
     * 1D interpolation in gsl requires periodic interpolation to have the
     * first and last values of the X & Y arrays equal to each other.  The 2D
     * interpolation here is periodic in the 2nd dimension, so we need to do a
     * copy of `arr` into a larger array `arr_periodic_y`, and copy the first
     * column to the new last column.  
     */
    for(i=0; i<nrows; i++) {
        for(j=0; j<ncols; j++) {
            arr_periodic_y[IDX(i,j,ncols_periodic)] = arr[IDX(i,j,ncols)];
        }
        /* set the last element in the column dimension. */
        arr_periodic_y[IDX(i,ncols_periodic-1,ncols_periodic)] = arr_periodic_y[IDX(i,0,ncols_periodic)];
    }

    /* Set-up the coordinate indices for the 0th and 1st dimensions -- x0 and
     * x1 -- minding the comment above about gsl's requirements for
     * periodicity.
     */
    x0 = (double*)malloc(nrows * sizeof(double));
    x1 = (double*)malloc(ncols_periodic * sizeof(double));
    if(NULL == x0 || NULL == x1)
        goto fail_coordinates;

    /* x0 = np.linspace(0.0, x0max, arr_periodic_y.shape[0], endpoint=False) */
    for(i=0; i<nrows; i++) {
        x0[i] = x0max * i / nrows;
    }

    /* x1 = np.arange(arr_periodic_y.shape[1], dtype=np.double)
     * x1_idx_max = np.max(x1)
     * x1 = (x1max / x1_idx_max) * x1
     */
    for(i=0; i<ncols_periodic; i++) {
        x1[i] = (x1max / (ncols_periodic-1)) * i;
    }

    interp2d = interp2d_alloc(gsl_interp_cspline_periodic, gsl_interp_cspline, nrows, ncols_periodic);
    if(NULL == interp2d) {
        goto fail_alloc_interp2d;
    }

    if(interp2d_init(interp2d, x0, x0max, x1, x1max, arr_periodic_y)) {
        goto fail_init_interp2d;
    }

    goto success;

fail_init_interp2d:
    if(interp2d) {
        interp2d_free(interp2d);
        interp2d = NULL;
    }
success:
fail_alloc_interp2d:
fail_coordinates:
    if(x0) {
        free(x0);
    }
    if(x1) {
        free(x1);
    }
    if(arr_periodic_y) {
        free(arr_periodic_y);
    }
    return interp2d;
}


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

    /* NOTE: should the interp2d->nrows be instead `interp2d->nrows-1` ??  */
    x0_idx = gsl_interp_bsearch(interp2d->global_x0, x, 0, interp2d->nrows);

    /* for now we only cover odd-length splines.  */
    assert(interp2d->col_spline_len % 2 == 1);

    lower_idx = x0_idx - (long)interp2d->col_spline_len / 2;
    upper_idx = x0_idx + (long)interp2d->col_spline_len / 2;

    /* initialize col_spline_x */
    for(i=0, idx=lower_idx; idx<=upper_idx; i++, idx++) {
        if(idx < 0) {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx + interp2d->nrows] - interp2d->max_x0;
        } else if(idx >= interp2d->nrows) {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx - interp2d->nrows] + interp2d->max_x0;
        } else {
            interp2d->col_spline_x[i] = interp2d->global_x0[idx];
        }
    }

    /* set col_spline_y array.  */
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
