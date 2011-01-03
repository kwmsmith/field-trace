/* |+ general interpolation object +| */
/* typedef struct { */
  /* gsl_interp * interp; */
  /* double  * x; */
  /* double  * y; */
  /* size_t  size; */
/* } gsl_spline; */

typedef struct {
    size_t ncols;
    size_t nrows;

    /* row_splines is an `nrows` length array of splines. Each spline has a
     * size of `ncols`.
     */
    gsl_spline *row_splines[];
    gsl_interp_accel *row_spline_accels[];

    size_t col_spline_len;
    gsl_spline *col_spline;
    gsl_interp_accel *col_spline_accel;
} interp2d_t;

/* LAYOUT:
 *
 * COLS ------------------------->
 * (0,0) -----------------> (0,NCOLS-1)        ROWS
 *   |                           |               |
 *   |                           |               |
 *   |                           |               |
 *   v                           v               |
 * (NROWS-1,0) ----------> (NROWS-1,NCOLS-1)     v
 *
 */

interp2d_t *
interp2d_alloc(const gsl_interp_type *T, size_t nrows, size_t ncols);

int
interp2d_init(interp2d_t *interp2d, const double *xa, const double *ya, size_t nrows, size_t ncols);

double
interp2d_eval(const interp2d_t *interp2d, double x, double y);

void
interp2d_free(interp2d_t *interp2d);

/* gsl_spline * */
/* gsl_spline_alloc(const gsl_interp_type * T, size_t size); */

/* int */
/* gsl_spline_init(gsl_spline * spline, const double xa[], const double ya[], size_t size); */

/*
 * xa and ya are pointers to 2D arrays of dimensions [nrows] X [ncols].
 */

/* const char * gsl_spline_name(const gsl_spline * spline); */
/* unsigned int gsl_spline_min_size(const gsl_spline * spline); */

/* int */
/* gsl_spline_eval_e(const gsl_spline * spline, double x, */
                  /* gsl_interp_accel * a, double * y); */

/* double */
/* gsl_spline_eval(const gsl_spline * spline, double x, gsl_interp_accel * a); */


/* int */
/* gsl_spline_eval_deriv_e(const gsl_spline * spline, */
                        /* double x, */
                        /* gsl_interp_accel * a, */
                        /* double * y); */

/* double */
/* gsl_spline_eval_deriv(const gsl_spline * spline, */
                      /* double x, */
                      /* gsl_interp_accel * a); */

/* int */
/* gsl_spline_eval_deriv2_e(const gsl_spline * spline, */
                         /* double x, */
                         /* gsl_interp_accel * a, */
                         /* double * y); */

/* double */
/* gsl_spline_eval_deriv2(const gsl_spline * spline, */
                       /* double x, */
                       /* gsl_interp_accel * a); */

/* int */
/* gsl_spline_eval_integ_e(const gsl_spline * spline, */
                        /* double a, double b, */
                        /* gsl_interp_accel * acc, */
                        /* double * y); */

/* double */
/* gsl_spline_eval_integ(const gsl_spline * spline, */
                      /* double a, double b, */
                      /* gsl_interp_accel * acc); */

/* void */
/* gsl_spline_free(gsl_spline * spline); */

