
#define REAL8 double
#define INTEGER int

/*
 * surfit and bispev subroutine interfaces.
 */

/*****************************************************************************
 * subroutine surfit(iopt,m,x,y,z,w,xb,xe,yb,ye,kx,ky,s,nxest,nyest,
 *          *  nmax,eps,nx,tx,ny,ty,c,fp,wrk1,lwrk1,wrk2,lwrk2,iwrk,kwrk,ier)
 *
 * real*8 xb,xe,yb,ye,s,eps,fp
 * integer iopt,m,kx,ky,nxest,nyest,nmax,nx,ny,lwrk1,lwrk2,kwrk,ier
 * real*8 x(m),y(m),z(m),w(m),tx(nmax),ty(nmax),c((nxest-kx-1)*(nyest-ky-1)),wrk1(lwrk1),wrk2(lwrk2)
 * integer iwrk(kwrk)
 *****************************************************************************/
void surfit_(INTEGER *iopt,
             INTEGER *m,
             REAL8 x[],
             REAL8 y[],
             REAL8 z[],
             REAL8 w[],
             REAL8 *xb,
             REAL8 *xe,
             REAL8 *yb,
             REAL8 *ye,
             INTEGER *kx,
             INTEGER *ky,
             REAL8 *s,
             INTEGER *nxest,
             INTEGER *nyest,
             INTEGER *nmax,
             REAL8 *eps,
             INTEGER *nx,
             REAL8 tx[],
             INTEGER *ny,
             REAL8 ty[],
             REAL8 c[],
             REAL8 *fp,
             REAL8 wrk1[],
             INTEGER *lwrk1,
             REAL8 wrk2[],
             INTEGER *lwrk2,
             INTEGER iwrk[],
             INTEGER *kwrk,
             INTEGER *ier);

/*****************************************************************************
 * subroutine bispev(tx,nx,ty,ny,c,kx,ky,x,mx,y,my,z,wrk,lwrk,
 *      * iwrk,kwrk,ier)
 *
 * integer nx,ny,kx,ky,mx,my,lwrk,kwrk,ier
 * integer iwrk(kwrk)
 * real*8 tx(nx),ty(ny),c((nx-kx-1)*(ny-ky-1)),x(mx),y(my),z(mx*my),wrk(lwrk)
 *****************************************************************************/
void bispev_(REAL8 tx[],
             INTEGER *nx,
             REAL8 ty[],
             INTEGER *ny,
             REAL8 c[],
             INTEGER *kx,
             INTEGER *ky,
             REAL8 x[],
             INTEGER *mx,
             REAL8 y[],
             INTEGER *my,
             REAL8 z[],
             REAL8 wrk[],
             INTEGER *lwrk,
             INTEGER iwrk[],
             INTEGER *kwrk,
             INTEGER *ier);
