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

#include "largefft.h"
#include <string.h>


/* blas fcopy */
/* y = x */
/* read x values spaced by Nx (so x size> N*Nx) */
/* write to y values spaced by Ny  (so y size > N*Ny) */
void
my_fcopy(int N, float *x, int Nx, float *y, int Ny) {
  extern void scopy_(int *N, float *x, int *incx, float *y, int *incy);
  /* use memcpy if Nx=Ny=1 */
  if (Nx==1&&Ny==1) {
   memcpy((void*)y,(void*)x,sizeof(float)*(size_t)N);
  } else {
   scopy_(&N,x,&Nx,y,&Ny);
  }
}

/* blas ccopy */
/* y = x */
/* read x values spaced by Nx (so x size> N*Nx) */
/* write to y values spaced by Ny  (so y size > N*Ny) */
void
my_ccopy(int N, complex float *x, int Nx, complex float *y, int Ny) {
  extern void ccopy_(int *N, complex float *x, int *incx, complex float *y, int *incy);
  /* use memcpy if Nx=Ny=1 */
  if (Nx==1&&Ny==1) {
   memcpy((void*)y,(void*)x,sizeof(complex float)*(size_t)N);
  } else {
   ccopy_(&N,x,&Nx,y,&Ny);
  }
}
/* blas scale */
/* x = a. x */
void
my_fscal(int N, float a, float *x) {
  extern void sscal_(int *N, float *alpha, float *x, int *incx);
  int i=1;
  sscal_(&N,&a,x,&i);
}
/* blas rotation */
/* y <= rot(x) */
void
my_frot(int N, float *x, float *y, float cosval, float sinval) {
  extern void srot_(int *N, float *x, int *incx, float *y, int *incy, float *csval, float *snval);
  int i=1;
  srot_(&N,x,&i,y,&i,&cosval,&sinval);
}
/* blas scale */
/* x = a. x */
void
my_csscal(int N, float a, complex float *x) {
  extern void csscal_(int *N, float *alpha, complex float *x, int *incx);
  int i=1;
  csscal_(&N,&a,x,&i);
}

/* ||x|| 2 norm */
float
my_fnrm2(int N, float *x) {
   extern float snrm2_(int *N, float *x, int *incx);
   int i=1;
   return snrm2_(&N,x,&i);
}

/* sum||x||_1 */
float 
my_fasum(int N, float *x) {
  extern float sasum_(int *N, float *x, int *incx);
  int i=1;
  return(sasum_(&N,x,&i));
}

/* BLAS y = a.x + y */
void
my_faxpy(int N, float *x, int incx, float a, float *y, int incy) {
    extern void saxpy_(int *N, float *alpha, float *x, int *incx, float *y, int *incy);
    saxpy_(&N,&a,x,&incx,y,&incy);
}

/* BLAS x^T*y */
float
my_fdot(int N, float *x, float *y) {
  extern float sdot_(int *N, float *x, int *incx, float *y, int *incy);
  int i=1;
  return(sdot_(&N,x,&i,y,&i));
}


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
int
my_dstemr(int N, double *D, double *E, double *W, double *Z, double *WORK, int LWORK, int *IWORK, int LIWORK) {
 extern void dstemr_(char *JOBZ, char *RANGE, int *N, double *D, double *E, double *VL, double *VU, int *IL, int *IU, int *M, double *W, double *Z,  int *LDZ,
        int *NZC,  int *ISUPPZ, int *TRYRAC,  double *WORK,  int *LWORK, int *IWORK, int *LIWORK, int *INFO);
  int info;
  char V='V';
  char RANGE='I';
  double VL=0.0; double VU=1e6; /* not used */
  int IL=1; int IU=N; /* find all eigenvalues/vectors */
  int M;
  int LDZ=N;
  int NZC=N;
  int *ISUPPZ;
  int TRYRAC=1;

  if ((ISUPPZ=(int*)calloc((size_t)(2*N),sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  dstemr_(&V, &RANGE, &N, D, E, &VL, &VU, &IL, &IU, &M, W, Z,  &LDZ,
        &NZC,  ISUPPZ, &TRYRAC,  WORK,  &LWORK, IWORK, &LIWORK, &info);
  free(ISUPPZ);
  return info;

}

/* LAPACK routine */
/* eigenvalue/eigenvector of a real symmetric */
/* N: full matrix size
   A: lower triangle of matrix in column major order
   W: output eigenvalues Nx1
   Z: output eigenvectors NxN 
   WORK: work size LWORK
   IWORK: work size LIWORK
*/
int
my_dsyevr(int N, double *A, double *W, double *Z, double *WORK, int LWORK, int *IWORK, int LIWORK) {
 extern void dsyevr_(char *JOBZ, char *RANGE, char *UPLO, int *N, double *A, int *LDA, double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, double *Z,  int *LDZ, 
        int *ISUPPZ,  double *WORK,  int *LWORK, int *IWORK, int *LIWORK, int *INFO);
  int info;
  char V='V';
  char RANGE='I';
  char UL='U';
  double VL=0.0; double VU=1e6; /* not used */
  int IL=1; int IU=N; /* find all eigenvalues/vectors */
  int M;
  int LDZ=N;
  int *ISUPPZ;
  double eps=1e-9;

  if ((ISUPPZ=(int*)calloc((size_t)(2*N),sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  dsyevr_(&V, &RANGE, &UL, &N, A, &N, &VL, &VU, &IL, &IU, &eps, &M, W, Z,  &LDZ,
         ISUPPZ,  WORK,  &LWORK, IWORK, &LIWORK, &info);
  free(ISUPPZ);
  return info;
}


/* BLAS SGEMV  y = alpha*op(A)*x+ beta*y : op 'T' or 'N' */
void
my_sgemv(char trans, int M, int N, float alpha, float *A, int lda, float *x, int incx, float beta, float *y, int incy) {
 extern void sgemv_(char *TRANS, int *M, int *N, float *ALPHA, float *A, int *LDA, float *X, int *INCX, float *BETA, float *Y, int *INCY);
 sgemv_(&trans, &M, &N, &alpha, A, &lda, x, &incx, &beta, y, &incy);
}

