#include <stdlib.h>

#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"

/*----------------------------------------------------------------------------
 * Layout for reference:
 *
 * Note that the first index specifies the rows, and the second index the
 * columns.  We use standard C array layout, or 'row major' storage; second
 * index varies fastest.
 *
 * COLS ---------------->
 * (0,0) ---------------> (0,NCOLS-1)  ROWS
 *   |                         |         |
 *   |                         |         |
 *   |                         |         |
 *   v                         v         v
 * (NROWS-1,0) --------> (NROWS-1,NCOLS-1)
 *
 */

typedef struct {
    size_t ncols;
    size_t nrows;

    /* row_splines is an `nrows` length array of splines. Each spline has a
     * size of `ncols`.
     */
    gsl_spline **row_splines;
    gsl_interp_accel **row_spline_accels;

    size_t col_spline_len;
    gsl_spline *col_spline;
    gsl_interp_accel *col_spline_accel;
    double *col_spline_x;
} interp2d_t;

interp2d_t *
interp2d_alloc(const gsl_interp_type *T, size_t nrows, size_t ncols);

int
interp2d_init(interp2d_t *interp2d, 
        const double *x0, const double *x1, 
        const double *ya, size_t nrows, size_t ncols);

double
interp2d_eval(const interp2d_t *interp2d, double x, double y);

void
interp2d_free(interp2d_t *interp2d);
