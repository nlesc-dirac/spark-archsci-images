/* Large FFT
 *
   Copyright (C) 2013 Sarod Yatawatta <sarod@users.sf.net>  
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 $Id$
*/


#ifndef LARGEFFT_H
#define LARGEFFT_H
#ifdef __cplusplus
        extern "C" {
#endif

#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdlib.h>

#include<wcs.h>

#ifndef MIN
#define MIN(x,y) \
  ((x)<=(y)? (x): (y))
#endif

#ifndef MAX
#define MAX(x,y) \
  ((x)>=(y)? (x): (y))
#endif

/* radians to degrees constants */
#ifndef CPI_180 /* pi/180 */
#define CPI_180  0.0174532925199433
#endif

#ifndef C180_PI /* 180/pi */
#define C180_PI 57.2957795130823 
#endif

/* for gcc 4.8 and above */
#ifndef complex
#define complex _Complex
#endif

/* tolerance for floating point operations */
#ifndef TOL
#define TOL 1e-12f
#endif

/* imaging modes */
#ifndef IMG_I0 /* only make I image, no grid, weight maps */
#define IMG_I0 0
#endif
#ifndef IMG_I /* make I image and extra maps of weights etc (GR,GI) */
#define IMG_I 1
#endif
#ifndef IMG_IQUV0 /* make IQUV maps, no extra info */
#define IMG_IQUV0 2
#endif
#ifndef IMG_IQUV /* make IQUV maps, extra maps of weights etc (GR,GI) */
#define IMG_IQUV 3
#endif
#ifndef IMG_IQUVF /* make IQUV maps, extra maps of weights etc (GR,GI) for IQUV */
#define IMG_IQUVF 4
#endif
/******************** mmio.c  ***********************/
/* create binary files needed */
/* uvgridname: binary file for uv gridded data
   d: pointer for uvgridded data
   fid0: file id
   imname: binary file for image 
   im: pointer for image
   fid1: file id
   wtgridname: binary file to store weights (and average images)
   fid2: file id
   wgt: pointer to weight grid
   d and im and wgt are in row major order
   Nx,Ny : grid size
   d: Nx Ny complex float elements  + 3 (Nx Ny) for IQUV images
   im: Nx Ny  complex float elements 
   wgt: Nx Ny  complex float elements (and average images)
   avgimg: nx ny float elements (only if snapshotmode>0)
   imgmode: IMG_I0,IMG_I,IMG_IQUV0,IMG_IQUV,IMG_IQUVF
*/
extern int
open_binary_file(const char *uvgridname, int *fid0, complex float **d, const char *imname, int *fid1, complex float **im, const char *wtgridname, int *fid2, complex float **wgt, int Nx, int Ny, int nx, int ny, int imgmode,int snapshotmode);
/* sync files */
extern int
sync_binary_file(complex float *d, complex float *im, complex float *wgt, int Nx, int Ny, int imgmode);
/* free files */
extern int
close_binary_file(const char *uvgridname, int fid0, complex float *d, const char *imname, int fid1, complex float *im,  const char *wtgridname, int fid2, complex float *wgt, int Nx, int Ny, int nx, int ny, int imgmode, int snapshotmode);
/* reset file to zero (for re-use) */
extern int
reset_binary_file(complex float *d,  complex float *im,  complex float *wgt, int Nx, int Ny,int imgmode);

/******************** fftutil.c **********************/
extern void
init_fftw(int Nt);
extern void
stop_fftw(void);
/* create,destroy plans */
/* Nt: no of threads */
extern int
do_create_fftw_plan(fftwf_plan *p, complex float *d, complex float *im, int n0, int n1);
extern int
do_destroy_fftw_plan(fftwf_plan p);

/* do FFT shift of input complex array, (row major order)
   also scale data points by scalefac
   and write back in column major order,
   real and imaginary parts seperately
   also performa apodization correction  
 
   din: input data array, size n0xn1 complex float
   dout: output data array, size 2*n0xn1 float
   din: row major, dout: col major 
   scalefac: scale for correction of FFT scaling
   M:  PSWF bandwidth in uv pixels, used for calculation of apodization correction
   Nx0: unpadded image size
   Npad: padding, so n0=n1=Nx0+2*Npad

   tol: cutoff tolerance for apodization correction

   output, the first Nx0*Nx0 values starting at dout[0] is the true image (real)
   and Nx0*Nx0 values starting at dout[2*n0*n1] is  the imaginary image
 */
extern int
do_data_fftshift(complex float *din, float *dout, float scalefac, int n0, int n1,int M, int Nx0, int Npad, float tol);

/* FFTshift the weights, remove zeropadding */
extern int
do_weight_fftshift(float *din, float *dout, int n0, int n1, int Nx0, int Npad);

/* FFTshift the gridded data, remove zeropadding
   write real and imag parts separately 
   din => dout
   scratch (n0xn1) is used for temp storage */
extern int
do_grid_fftshift(complex float *din, float *dout,  float scalefac, float *scratch, int n0, int n1, int Nx0, int Npad);



/* IFFTshift to undo the FFTshift */
/* N even
    [0,N/2-1] -> [N/2,N-1]  (N/2 values) and [N/2,N-1] -> [0,N/2-1] (N/2 values)
   N odd
    [0,N/2-1] -> [N/2+1,N-1]  (N/2 values) and [N/2,N-1] -> [0,N/2]  (N/2+1 values)
*/
extern int
ifftshift_index(int N, int x);

/* create,destroy plans */
/* FFT din => dout */
extern int
do_create_fftw_plan_2d(fftwf_plan *p, complex float *din, complex float *dout, int n0, int n1, int Nt);

/* in-place fftshift */
/* in: NxN array, row major order
   N0: location of zero pixel
   [N0,N-1] -> [0,N-1-N0] and [0,N0-1] -> [N-N0,N-1]
*/
extern int
do_fftshift_inplace(complex float *in, int N, int N0);

/* in-place ifftshift */
/* in: NxN array, row major order
   N0: location of zero pixel
   [0,N-1-N0] -> [N0,N-1]  and [N-N0,N-1] -> [0,N0-1]  
*/
extern int
do_ifftshift_inplace(complex float *in, int N, int N0);
/******************** myblasf.c **********************/
/* y = x */
/* read x values spaced by Nx (so x size> N*Nx) */
/* write to y values spaced by Ny  (so y size > N*Ny) */
extern void
my_fcopy(int N, float *x, int Nx, float *y, int Ny);

/* blas ccopy */
/* y = x */
/* read x values spaced by Nx (so x size> N*Nx) */
/* write to y values spaced by Ny  (so y size > N*Ny) */
extern void
my_ccopy(int N, complex float *x, int Nx, complex float *y, int Ny);

/* blas scale */
/* x = a. x */
extern void
my_fscal(int N, float a, float *x);
/* blas rotation */
/* y <= rot(x) */
extern void
my_frot(int N, float *x, float *y, float cosval, float sinval);
/* blas scale */
/* x = a. x */
extern void
my_csscal(int N, float a, complex float *x);
/* ||x|| 2 norm */
extern float
my_fnrm2(int N, float *x);
/* sum||x||_1 */
float 
my_fasum(int N, float *x);
/* BLAS y = a.x + y */
extern void
my_faxpy(int N, float *x, int incx, float a, float *y, int incy);

/* BLAS x^T*y */
float
my_fdot(int N, float *x, float *y);

/* LAPACK routine */
/* eigenvalue/eigenvector of a symmetric tridiagonal matrix */
/* N: full matrix size
   D: main diagonal Nx1
   E: subdiagonal Nx1, only first N-1 values 
   W: output eigenvalues Nx1
   Z: output eigenvectors NxN 
   WORK: work size LWORK
   IWORK: work size LIWORK
*/
extern int
my_dstemr(int N, double *D, double *E, double *W, double *Z, double *WORK, int LWORK, int *IWORK, int LIWORK);

/* eigenvalue/eigenvector of a real symmetric */
/* N: full matrix size
   A: lower triangle of matrix in column major order
   W: output eigenvalues Nx1
   Z: output eigenvectors NxN 
   WORK: work size LWORK
   IWORK: work size LIWORK
*/
extern int
my_dsyevr(int N, double *A, double *W, double *Z, double *WORK, int LWORK, int *IWORK, int LIWORK);


/* BLAS SGEMV  y = alpha*op(A)*x+ beta*y : op 'T' or 'N' */
extern void
my_sgemv(char trans, int M, int N, float alpha, float *A, int lda, float *x, int incx, float beta, float *y, int incy);

/******************** myblas2.c **********************/
/* BLAS wrappers */
/* blas dcopy */
/* y = x */
/* read x values spaced by Nx (so x size> N*Nx) */
/* write to y values spaced by Ny  (so y size > N*Ny) */
extern void
my_dcopy(int N, double *x, int Nx, double *y, int Ny);

/* blas scale */
/* x = a. x */
extern void
my_dscal(int N, double a, double *x);

/* x^T*y */
extern double
my_ddot(int N, double *x, double *y);

/* ||x||_2 */
extern double
my_dnrm2(int N, double *x);

/* sum||x||_1 */
extern double
my_dasum(int N, double *x);

/* BLAS y = a.x + y */
extern void
my_daxpy(int N, double *x, double a, double *y);

/* BLAS y = a.x + y */
extern void
my_daxpys(int N, double *x, int incx, double a, double *y, int incy);

/* max |x|  index (start from 1...)*/
extern int
my_idamax(int N, double *x, int incx);
extern int
my_icamax(int N, complex float *x, int incx);

/* BLAS DGEMM C = alpha*op(A)*op(B)+ beta*C */
extern void
my_dgemm(char transa, char transb, int M, int N, int K, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc);

/* BLAS DGEMV  y = alpha*op(A)*x+ beta*y : op 'T' or 'N' */
extern void
my_dgemv(char trans, int M, int N, double alpha, double *A, int lda, double *x, int incx,  double beta, double *y, int incy);

/* following routines used in LAPACK solvers */
/* cholesky factorization: real symmetric */
extern int
my_dpotrf(char uplo, int N, double *A, int lda);

/* solve Ax=b using cholesky factorization */
extern int
my_dpotrs(char uplo, int N, int nrhs, double *A, int lda, double *b, int ldb);

/* solve Ax=b using QR factorization */
extern int
my_dgels(char TRANS, int M, int N, int NRHS, double *A, int LDA, double *B, int LDB, double *WORK, int LWORK);

/* A=U S VT, so V needs NOT to be transposed */
extern int
my_dgesvd(char JOBU, char JOBVT, int M, int N, double *A, int LDA, double *S,
   double *U, int LDU, double *VT, int LDVT, double *WORK, int LWORK);

/******************** pswf.c **********************/
/* binary search */
 /*
 * xarr: array of sorted values, make sure x is within the range 
 * x: value to search
 * i_start: starting index of array to search 
 * i_end: end index of array to search 
 *
 * return value: index k, such that xarr[k]<= x< xarr[k+1]
 * for errors: return negative values
 */
extern int 
bin_search(float *xarr,float x,int i_start,int i_end);

/* 
   c>0 real, bandwidth 0<c<=(pi/2)*(1+0.5)
   for VLA memo, c=pi*m/2, m: convolution support in uv pixels
   x  array of Nx1 values in [0,1]  (PSWF is symmetric, so only half is calculated)
   y output array of Nx1 PSWF values, normalized such that y[0]=1.0
*/
extern int
pswf_nn(float c, int N, float *x, float *y);

/* given arrays x,y for PSWF(x) values, and given z,
   find the interpolated value */
/* x : (N+2)x1 with guard values
   y : size Nx1 
   x: grid, y=PSWF(x) */
extern float
pswf_eval(int N,float *x, float *y, float z);

/* 
  ugrid: array (N+2)x1 of uv pixel values with guard
  wkernel: NxN x Nw kernel values
  k: lower bound w plane, such that w[k] <= |w| <= w[k+1]
*/
extern complex float
conv_eval(int N, float *ugrid, complex float *wkernel, float *wparr, float z, int Nw, int k, float w,float deltau,float deltav);


/* 
  ugrid: array (N+2)x1 of uv pixel values with guard
  kernel: NxN kernel values
  deltau, deltav: (u-u_i), (v-v_i) of this pixel
*/
float
conv_plane_eval(int N, float *ugrid, float *kernel, float deltau, float deltav);


/* 
  create convolution kernel to be written as a FITS file 
  buff: image buffer: size NxN
  N: image size: NxN
  M: conv. kernel bandwidth in uvpixels 
*/
extern int
create_conv_kernel(float *buff, int N, int M);
/******************** hankel.c **********************/
/* setup matrix C and multiplication vector m
   for Hankel transform
   N: no. of data points <= 3000
   C: NxN matrix, column major order
   mr: Nx1 pre-multiplication vector for data
   mv: Nx1 post-multiplication vector for transformed data
   R: support in time, data points in [0,R]
   transformed data points are in [0,1/(2 pi R)]
*/
extern int
setup_hankel_tf(int N, float *C, float *mr, float *mv, float R);


/* perform Hankel transform of data vector x to get y
   C: NxN matrix
   mr,mv: Nx1 vectors
   x: Nx1 data
   y: Nx1 output
   y=(C*(x.*mr)).*mv;
*/
extern int
do_hankel_tf(int N, float *C, float *mr, float *mv, float *x, float *y);

/******************** wkernel.c **********************/
extern int
generate_wconv_kernels(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, float *wkernel,  float tol, int *support);

/******************** wkernel2.c **********************/
extern int
generate_2d_wconv_kernels(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel,  float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot, int Nt);
extern int
generate_2d_wconv_kernels_over(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel,  float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot, double r, int Nt);
extern int
generate_2d_wconv_kernels_over_gpu(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel,  float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot, double r, int Nt);

/******************** lmwcs.c **********************/
/* 
  create a default WCS structure to calculate l,m coords
  deltalm: pixel width (deg)
  int N0: reference pixel 
  ra0,dec0: reference coordinates (deg)
  only if gid==0, enable memory allocation for WCS
*/
extern int
generate_def_wcs(struct wcsprm *wcs, double deltalm, double N0, double ra0, double dec0,int gid);


/******************** imagetf.c *********************/
/* transform image pixels
  l' = l +a(sqrt(1-l^2-m^2)-1)
  m' = m +b(sqrt(1-l^2-m^2)-1)
  din is in l',m' grid, with pixel values
  dout is l,m grid, with no pixel value
  update dout pixel values using din pixel values and right coords
  both din,dout size Nx x Ny

  delta_l: pixel size (rad)
  ra0,dec0: phase center (rad)
  Nt: no. threads
*/
extern int
do_image_lm_transform(float *din,float *dout,int Nx,int Ny,double a,double b,double delta_l,double ra0, double dec0, int Nt);

/****************************** clmfit_nocuda.c ****************************/
/* fit Gaussian PSF
   psf: PSF image pixels (peak normalized to 1)
   lgrid,mgrid: grid values
   psf,lgrid,mgrid:  N x 1 vectors 
   
   fit a function exp(-l^2/bmaj^2-m^2/bmin^2)
   where 
   l=lgrid*cos(bpa)+mgrid*sin(bpa)
   m=-lgrid*sin(bpa)+mgrid*cos(bpa)
*/
int
fit_gaussian_psf(double *psf,double *lgrid,double *mgrid, int N, double *bmaj, double *bmin, double *bpa);

/********* constants - from levmar ******************/
#define CLM_INIT_MU       1E-03
#define CLM_STOP_THRESH   1E-17
#define CLM_DIFF_DELTA    1E-06
#define CLM_EPSILON       1E-12
#define CLM_ONE_THIRD     0.3333333334 /* 1.0/3.0 */
#define CLM_OPTS_SZ       5 /* max(4, 5) */
#define CLM_INFO_SZ       10

/****************************** pcgsolve.c ****************************/
/* a function to solve the matrix equation by conjugate gradient 
 method  Ax=y */
extern void
solve_conjugate(float *A, float *x, float *y, int N);

/* Non-negative least squares routine
  Author:   Yuancheng Luo
  E-mail:   yluo1@mail.umd.edu
*/
/* A * x =b, NNLS solution
   A : size m x n
   nSys: Number of systems to solve = 1
   isTransposed: A is transposed in memory = 0
*/
extern void
nnlsOMPSysMKLUpdates(float *A, float *b, float *x, int isTransposed, int maxNNLSIters, int maxLSIters, int nSys, int m, int n, float TOL_TERMINATION);

#ifdef __cplusplus
     } /* extern "C" */
#endif
#endif /* LARGEFFT_H */
