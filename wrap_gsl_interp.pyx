cdef extern from "gsl/gsl_interp.h":
    ctypedef struct gsl_interp_accel:
        pass
    ctypedef struct gsl_interp_type:
        pass
    gsl_interp_accel *gsl_interp_accel_alloc()
    void gsl_interp_accel_free(gsl_interp_accel * a)

cdef extern from "gsl/gsl_spline.h":
    ctypedef struct gsl_spline:
        pass
    gsl_interp_type *gsl_interp_cspline
    gsl_spline *gsl_spline_alloc(gsl_interp_type *, size_t)
    int gsl_spline_init(gsl_spline * spline, double xa[], double ya[], size_t size)
    double gsl_spline_eval(gsl_spline * spline, double x, gsl_interp_accel * a)
    void gsl_spline_free(gsl_spline * spline)

import numpy as np
cimport numpy as np


cdef class Interp:
    cdef gsl_spline *spline
    cdef gsl_interp_accel *acc
    cdef size_t size

    def __cinit__(self, np.ndarray[double, ndim=1] xa, np.ndarray[double, ndim=1] ya):
        self.size = len(xa)
        if len(xa) != len(ya):
            raise ValueError("len(xa) != len(ya)")
        self.acc = gsl_interp_accel_alloc()
        self.spline = gsl_spline_alloc(gsl_interp_cspline, self.size)
        gsl_spline_init(self.spline, <double*>xa.data, <double*>ya.data, self.size)

    cpdef double eval(self, double xi):
        return gsl_spline_eval(self.spline, xi, self.acc)

    def __dealloc__(self):
        gsl_spline_free(self.spline)
        gsl_interp_accel_free(self.acc)

cdef class Interp2D:
    cdef int nx, ny, npoints
    cdef list interp_cols

    def __init__(self, np.ndarray[double, ndim=2] a, npoints=3):
        self.npoints = npoints
        self.nx = a.shape[0]
        self.ny = a.shape[1]
        Yidxs = np.arange(self.ny, dtype=np.double)
        self.interp_cols = []
        for row in a:
            self.interp_cols.append(Interp(Yidxs, row))

    cpdef double eval(self, double xloc, double yloc):
        # can be sped up.  xis can be just a few elements long...
        cdef int idx
        cdef int ctr = int(round(xloc))
        idxs = range(ctr-self.npoints//2, ctr+self.npoints//2+1)
        cdef np.ndarray[double, ndim=1] xis = np.array(idxs, dtype=np.double)
        cdef np.ndarray[double, ndim=1] yis = np.empty(self.npoints, dtype=np.double)
        for idx in range(len(idxs)):
            yis[idx] = self.interp_cols[idxs[idx]].eval(yloc)
        interp_row = Interp(xis, yis)
        return interp_row.eval(xloc)

'''
/* general interpolation object */
typedef struct {
  gsl_interp * interp;
  double  * x;
  double  * y;
  size_t  size;
} gsl_spline;

gsl_spline *
gsl_spline_alloc(const gsl_interp_type * T, size_t size);
     
int
gsl_spline_init(gsl_spline * spline, const double xa[], const double ya[], size_t size);

const char * gsl_spline_name(const gsl_spline * spline);
unsigned int gsl_spline_min_size(const gsl_spline * spline);


int
gsl_spline_eval_e(const gsl_spline * spline, double x,
                  gsl_interp_accel * a, double * y);

double
gsl_spline_eval(const gsl_spline * spline, double x, gsl_interp_accel * a);

int
gsl_spline_eval_deriv_e(const gsl_spline * spline,
                        double x,
                        gsl_interp_accel * a,
                        double * y);

double
gsl_spline_eval_deriv(const gsl_spline * spline,
                      double x,
                      gsl_interp_accel * a);

int
gsl_spline_eval_deriv2_e(const gsl_spline * spline,
                         double x,
                         gsl_interp_accel * a,
                         double * y);

double
gsl_spline_eval_deriv2(const gsl_spline * spline,
                       double x,
                       gsl_interp_accel * a);

int
gsl_spline_eval_integ_e(const gsl_spline * spline,
                        double a, double b,
                        gsl_interp_accel * acc,
                        double * y);

double
gsl_spline_eval_integ(const gsl_spline * spline,
                      double a, double b,
                      gsl_interp_accel * acc);

void
gsl_spline_free(gsl_spline * spline);

__END_DECLS

#endif /* __GSL_INTERP_H__ */
'''

'''
/* evaluation accelerator */
typedef struct {
  size_t  cache;        /* cache of index   */
  size_t  miss_count;   /* keep statistics  */
  size_t  hit_count;
}
gsl_interp_accel;


/* interpolation object type */
typedef struct {
  const char * name;
  unsigned int min_size;
  void *  (*alloc) (size_t size);
  int     (*init)    (void *, const double xa[], const double ya[], size_t size);
  int     (*eval)    (const void *, const double xa[], const double ya[], size_t size, double x, gsl_interp_accel *, double * y);
  int     (*eval_deriv)  (const void *, const double xa[], const double ya[], size_t size, double x, gsl_interp_accel *, double * y_p);
  int     (*eval_deriv2) (const void *, const double xa[], const double ya[], size_t size, double x, gsl_interp_accel *, double * y_pp);
  int     (*eval_integ)  (const void *, const double xa[], const double ya[], size_t size, gsl_interp_accel *, double a, double b, double * result);
  void    (*free)         (void *);

} gsl_interp_type;


/* general interpolation object */
typedef struct {
  const gsl_interp_type * type;
  double  xmin;
  double  xmax;
  size_t  size;
  void * state;
} gsl_interp;


/* available types */
GSL_VAR const gsl_interp_type * gsl_interp_linear;
GSL_VAR const gsl_interp_type * gsl_interp_polynomial;
GSL_VAR const gsl_interp_type * gsl_interp_cspline;
GSL_VAR const gsl_interp_type * gsl_interp_cspline_periodic;
GSL_VAR const gsl_interp_type * gsl_interp_akima;
GSL_VAR const gsl_interp_type * gsl_interp_akima_periodic;

gsl_interp_accel *
gsl_interp_accel_alloc(void);

int
gsl_interp_accel_reset (gsl_interp_accel * a);

void
gsl_interp_accel_free(gsl_interp_accel * a);

gsl_interp *
gsl_interp_alloc(const gsl_interp_type * T, size_t n);
     
int
gsl_interp_init(gsl_interp * obj, const double xa[], const double ya[], size_t size);

const char * gsl_interp_name(const gsl_interp * interp);
unsigned int gsl_interp_min_size(const gsl_interp * interp);


int
gsl_interp_eval_e(const gsl_interp * obj,
                  const double xa[], const double ya[], double x,
                  gsl_interp_accel * a, double * y);

double
gsl_interp_eval(const gsl_interp * obj,
                const double xa[], const double ya[], double x,
                gsl_interp_accel * a);

int
gsl_interp_eval_deriv_e(const gsl_interp * obj,
                        const double xa[], const double ya[], double x,
                        gsl_interp_accel * a,
                        double * d);

double
gsl_interp_eval_deriv(const gsl_interp * obj,
                      const double xa[], const double ya[], double x,
                      gsl_interp_accel * a);

int
gsl_interp_eval_deriv2_e(const gsl_interp * obj,
                         const double xa[], const double ya[], double x,
                         gsl_interp_accel * a,
                         double * d2);

double
gsl_interp_eval_deriv2(const gsl_interp * obj,
                       const double xa[], const double ya[], double x,
                       gsl_interp_accel * a);

int
gsl_interp_eval_integ_e(const gsl_interp * obj,
                        const double xa[], const double ya[],
                        double a, double b,
                        gsl_interp_accel * acc,
                        double * result);

double
gsl_interp_eval_integ(const gsl_interp * obj,
                      const double xa[], const double ya[],
                      double a, double b,
                      gsl_interp_accel * acc);

void
gsl_interp_free(gsl_interp * interp);

INLINE_DECL size_t
gsl_interp_bsearch(const double x_array[], double x,
                   size_t index_lo, size_t index_hi);

#ifdef HAVE_INLINE

/* Perform a binary search of an array of values.
 * 
 * The parameters index_lo and index_hi provide an initial bracket,
 * and it is assumed that index_lo < index_hi. The resulting index
 * is guaranteed to be strictly less than index_hi and greater than
 * or equal to index_lo, so that the implicit bracket [index, index+1]
 * always corresponds to a region within the implicit value range of
 * the value array.
 *
 * Note that this means the relationship of 'x' to x_array[index]
 * and x_array[index+1] depends on the result region, i.e. the
 * behaviour at the boundaries may not correspond to what you
 * expect. We have the following complete specification of the
 * behaviour.
 * Suppose the input is x_array[] = { x0, x1, ..., xN }
 *    if ( x == x0 )           then  index == 0
 *    if ( x > x0 && x <= x1 ) then  index == 0, and sim. for other interior pts
 *    if ( x == xN )           then  index == N-1
 *    if ( x > xN )            then  index == N-1
 *    if ( x < x0 )            then  index == 0 
 */

INLINE_FUN size_t
gsl_interp_bsearch(const double x_array[], double x,
                   size_t index_lo, size_t index_hi)
{
  size_t ilo = index_lo;
  size_t ihi = index_hi;
  while(ihi > ilo + 1) {
    size_t i = (ihi + ilo)/2;
    if(x_array[i] > x)
      ihi = i;
    else
      ilo = i;
  }
  
  return ilo;
}
#endif

INLINE_DECL size_t 
gsl_interp_accel_find(gsl_interp_accel * a, const double x_array[], size_t size, double x);

#ifdef HAVE_INLINE
INLINE_FUN size_t
gsl_interp_accel_find(gsl_interp_accel * a, const double xa[], size_t len, double x)
{
  size_t x_index = a->cache;
 
  if(x < xa[x_index]) {
    a->miss_count++;
    a->cache = gsl_interp_bsearch(xa, x, 0, x_index);
  }
  else if(x >= xa[x_index + 1]) {
    a->miss_count++;
    a->cache = gsl_interp_bsearch(xa, x, x_index, len-1);
  }
  else {
    a->hit_count++;
  }
  
  return a->cache;
}
#endif /* HAVE_INLINE */


__END_DECLS

#endif /* __GSL_INTERP_H__ */
'''
