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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "largefft.h"

/* reference:
  H Xiao, V Rokhlin and N Yarvin
  Prolate spheroidal wavefunctions, quadrature and
interpolation
  Inverse Problems 17 (2001) 805–838
*/

/* Legendre polynomial x in [-1,1] */
/* M : how many terms 0..M 
    P: size M+1 x 1
*/
static int 
Lfn(int M, double x, double *P){
 int ci;
 P[0]=1.0/sqrt(2.0);
 P[1]=sqrt(1.5)*x;
 for (ci=1; ci<M; ci++) {
  P[ci+1]=sqrt((double)((2*(ci+1)+1)*(2*(ci+1)-1)))*x*P[ci]/(double)(ci+1)-sqrt((double)((2*(ci+1)+1)/(double)(2*(ci+1)-3)))*(double)(ci)*P[ci-1]/(double)(ci+1);
 }
 return 0;
}


/* c>0 real, bandwidth 0<c<=(pi/2)*(1+0.5)
   x in [-1,1]
   for VLA memo, c=pi*m/2, m: convolution support
   beta: expansion coefficients, size m+1 x 1
*/
static int
setup_pswf_nn(int m, double c, double* beta) {
 int k; 
 
 double *A,*Z,*WORK,*W;
 int *IWORK;
 int info,work_sz,iwork_sz;
 int M=m+1;
 if ((A=(double*)calloc((size_t)M*M,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((Z=(double*)calloc((size_t)M*M,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((W=(double*)calloc((size_t)M,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((WORK=(double*)calloc((size_t)2,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((IWORK=(int*)calloc((size_t)2,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }

 /* workspace query */
 info=my_dsyevr(M, A, W, Z, WORK, -1, IWORK, -1);
 if (info) {
  fprintf(stderr,"%s: %d: LAPACK error\n",__FILE__,__LINE__);
  exit(1);
 }
 work_sz=(int)WORK[0];
 iwork_sz=(int)IWORK[0];
 //printf("info=%d, work_sz=%d iwork_sz=%d\n",info,work_sz,iwork_sz);
 free(WORK);
 if ((WORK=(double*)calloc((size_t)work_sz,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 free(IWORK);
 if ((IWORK=(int*)calloc((size_t)iwork_sz,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 double c2=c*c;
 /* setup main/super diagonals */
 for (k=0; k<M-2; k++) {
   /* A(k,k) */
   A[k*M+k]=(double)(k*(k+1))+(double)(2*k*(k+1)-1)/((double)((2*k+3)*(2*k-1)))*c2;
   /* A(k,k+2) */
   A[(k+2)*M+k]=(double)((k+2)*(k+1))/((double)((2*k+3))*sqrt((double)((2*k+1)*(2*k+5))))*c2;
   //printf("%d %lf %lf\n",k,A[k*M+k],A[(k+2)*M+k]);
 }
 for (k=M-2; k<M; k++) {
   /* A(k,k) */
   A[k*M+k]=(double)(k*(k+1))+(double)(2*k*(k+1)-1)/((double)((2*k+3)*(2*k-1)))*c2;
 }
 info=my_dsyevr(M, A, W, Z, WORK, work_sz, IWORK, iwork_sz);
 if (info) {
  fprintf(stderr,"%s: %d: LAPACK error\n",__FILE__,__LINE__);
  exit(1);
 }
 /* find min eigenvalue */
 double min_e=1e9;
 int min_idx=0;
 for (k=0; k<M; k++) {
  if (min_e>W[k]) {
   min_e=W[k];
   min_idx=k;
   //printf("eigenvalue %d =%lf\n",k,W[k]);
  }
 }
 //printf("min eigenvalue %d=%lf\n",min_idx,min_e);
 /* copy values for max eigenvalue to beta */
 for (k=0; k<M; k++) {
  //printf("%d : %lf\n",k,Z[min_idx*M+k]);
  beta[k]=Z[min_idx*M+k];
 }
 free(A);
 free(Z);
 free(W);
 free(WORK);
 free(IWORK);
 return 0;
}




/* 
   c>0 real, bandwidth 0<c<=(pi/2)*(1+0.5)
   for VLA memo, c=pi*m/2, m: convolution support in uv pixels
   x  array of Nx1 values in [0,1]  (PSWF is symmetric, so only half is calculated)
   y output array of Nx1 PSWF values, normalized such that y[0]=1.0
*/
int
pswf_nn(float c, int N, float *x, float *y) {

 /* expansion coefficients for P_{m-14},P_{m-12}..P_{m}....P_{m+12},P_{m+14}
   depending on m, some of them are zero as P_{n} exists for n>=0 */
 /* now evaluate */
 double *beta,*Px;
 int M=32; /* M=2*N+30, N=No. of pswf modes = 1 */
 if ((beta=(double*)calloc((size_t)(M+1),sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((Px=(double*)calloc((size_t)(M+1),sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 setup_pswf_nn(M, (double)c, beta);
 int ci;
 ci=0;
 Lfn(M, x[ci], Px);
 y[ci]=my_ddot(M+1,beta,Px);
 /* scale y such that y[0]=1.0 */
 for (ci=1; ci<N; ci++) {
  Lfn(M, x[ci], Px);
  y[ci]=my_ddot(M+1,beta,Px)/y[0];
 }
 y[0]=1.0f;
 free(beta);
 free(Px);
 return 0;
}


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
int 
bin_search1(float *__restrict xarr,float x,int i_start,int i_end) {
 int i;
 //trivial case
 if (i_start==i_end) {
    if (xarr[i_start]==x)
        return i_start;
    else {
      fprintf(stderr,"%s: %d: ERROR: bin search error 1\n",__FILE__,__LINE__);
      return -1;
    }
  }
 //trivial case with length 2 array
  if (i_end==i_start+1) {
    if (x>=xarr[i_start] && x<xarr[i_end]) {
        return i_start;
    } else {
      fprintf(stderr,"%s: %d: ERROR: bin search error 2\n",__FILE__,__LINE__);
      return -2;
    }
  }

  //compare the mid point
  i=(int)((i_start+i_end)/2);
  if (x>=xarr[i] && x<xarr[i+1]) {
     return i;
  } else {
    //go to lower half of the upper half of the array
    if (x<xarr[i])
      return bin_search1(xarr,x,i_start,i);
    else
      return bin_search1(xarr,x,i,i_end);
  }
  //will not reach here
  fprintf(stderr,"%s: %d: ERROR: bin search error 3\n",__FILE__,__LINE__);
  return -3;
}

int
bin_search(float *__restrict xarr,float x,int i_start,int i_end) {
  int low = i_start;
  int high = i_end-1;
  while (low <= high) {
      int mid = (low + high) / 2;
      if ((x>=xarr[mid] && x<xarr[mid+1])) return mid;
      else if (xarr[mid] < x) low = mid + 1;
      else high = mid - 1;
  }
  return -1;
}

/* given arrays x,y for PSWF(x) values, and given z,
   find the interpolated value */
/* x: (N+2)x1 with guard values
   y : size Nx1 
   x: grid, y=PSWF(x) */
float
pswf_eval(int N,float *x, float *y, float z) {
  
  /* find correct grid value */
  int k=bin_search(x,z,0,N+2);
  /* now x[k] <= z < x[k+1] */
  /* x[0]  x[1]  x[2] ..... x[N]     x[N+1] 
      XX   y[0]  y[1] .....   y[N-1]   XX
  */
  if (k==0) {
   return y[k];
  }
  if (k==N+1) {
   return y[k-2];
  } 
  k=k-1;
  /* use bilinear interpolation */
  return y[k]+(z-x[k+1])*(y[k+1]-y[k])/(x[k+2]-x[k+1]);
}


/* 
  ugrid: array (N+2)x1 of uv pixel values with guard
  wkernel: NxN x Nw kernel values
  wparr: (Nw+2)x1 array of w values, with guard
  z: (abs(w))
  k: lower bound w plane, such that w[k] <= |w| <= w[k+1]
  w: w value of this pixel
  deltau, deltav: (u-u_i), (v-v_i) of this pixel
*/
complex float
conv_eval(int N, float *ugrid, complex float *wkernel, float *wparr, float z, int Nw, int k, float w, float deltau, float deltav) {
  complex float ppxy;
  /* find correct grid value */
  int ui=bin_search(ugrid,deltau,0,N+2);
  int vi=bin_search(ugrid,deltav,0,N+2);
  /* now pixel is in [ui,ui+1] and [vi,vi+1],
    but since 0 is a guard value, shift down pixels */
  ui--;
  vi--;
  if (ui<0 || vi<0 || ui>=N-1 || vi>=N-1) { 
//    printf("%f : %d in [%f,%f] %f: %d in [%f,%f]\n",deltau,ui,ugrid[ui+1],ugrid[ui+2],deltav,vi,ugrid[vi+1],ugrid[vi+2]);
    return 0.0f+_Complex_I*0.0f; 
  }

  complex float c11=wkernel[k*N*N+ui*N+vi];
  complex float c21=wkernel[k*N*N+(ui+1)*N+vi];
  float uratio=(deltau-ugrid[ui+1])/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  //complex float c1=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  complex float c1=c11+uratio*(c21-c11);
  c11=wkernel[k*N*N+ui*N+vi+1];
  c21=wkernel[k*N*N+(ui+1)*N+vi+1];
  //complex float c2=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  complex float c2=c11+uratio*(c21-c11);

  float vratio=(deltav-ugrid[vi+1])/(TOL+ugrid[vi+2]-ugrid[vi+1]);
  //complex float c3=c1+(deltav-ugrid[vi+1])*(c2-c1)/(TOL+ugrid[vi+2]-ugrid[vi+1]);
  complex float c3=c1+vratio*(c2-c1);

  c11=wkernel[(k+1)*N*N+ui*N+vi];
  c21=wkernel[(k+1)*N*N+(ui+1)*N+vi];
  //c1=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  c1=c11+uratio*(c21-c11);
  c11=wkernel[(k+1)*N*N+ui*N+vi+1];
  c21=wkernel[(k+1)*N*N+(ui+1)*N+vi+1];
  //c2=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  c2=c11+uratio*(c21-c11);

  //complex float c4=c1+(deltav-ugrid[vi+1])*(c2-c1)/(TOL+ugrid[vi+2]-ugrid[vi+1]);
  complex float c4=c1+vratio*(c2-c1);

  ppxy=c3+(c4-c3)*(z-wparr[k+1])/(TOL+wparr[k+2]-wparr[k+1]);
  if (w<0.0f) {
   return conjf(ppxy);
  }
  return ppxy;
}



/* 
  ugrid: array (N+2)x1 of uv pixel values with guard
  kernel: NxN kernel values
  deltau, deltav: (u-u_i), (v-v_i) of this pixel
*/
float
conv_plane_eval(int N, float *__restrict ugrid, float *__restrict kernel, float deltau, float deltav) {
  /* find correct grid value */
  int ui=bin_search(ugrid,deltau,0,N+2);
  int vi=bin_search(ugrid,deltav,0,N+2);

  /* now pixel is in [ui,ui+1] and [vi,vi+1],
    but since 0 is a guard value, shift down pixels */
  ui--;
  vi--;
  if (ui<0 || vi<0 || ui>=N-1 || vi>=N-1) { 
    return 0.0f; 
  }
  //printf("%f : %d in [%f,%f] %f: %d in [%f,%f]\n",deltau,ui,ugrid[ui+1],ugrid[ui+2],deltav,vi,ugrid[vi+1],ugrid[vi+2]);

  float c11=kernel[ui*N+vi];
  float c21=kernel[(ui+1)*N+vi];
  float uratio=(deltau-ugrid[ui+1])/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  //float c1=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  float c1=c11+uratio*(c21-c11);
  //printf("c11=%f c21=%f c1=%f ratio=%f\n",c11,c21,c1,uratio);
  c11=kernel[ui*N+vi+1];
  c21=kernel[(ui+1)*N+vi+1];
  //float c2=c11+(deltau-ugrid[ui+1])*(c21-c11)/(TOL+ugrid[ui+2]-ugrid[ui+1]);
  float c2=c11+uratio*(c21-c11);
  //printf("c11=%f c21=%f c2=%f\n",c11,c21,c1);

  float c3=c1+(deltav-ugrid[vi+1])*(c2-c1)/(TOL+ugrid[vi+2]-ugrid[vi+1]);

  return c3;
}





/* 
  create convolution kernel to be written as a FITS file 
  buff: image buffer: size NxN
  N: image size: NxN
  M: conv. kernel bandwidth in uvpixels 
*/
int
create_conv_kernel(float *buff, int N, int M) {
  int ci,cj;

  float pswf_c=M_PI*0.5f*(float)M;
  int N0=N/2; /* 0 pixel */
  int Np=N/2+1;
  float pswf_delta=1.0f/((float)(N/2));
  float *pswf_x,*pswf_y;
  if ((pswf_x=(float*)calloc((size_t)Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_y=(float*)calloc((size_t)Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  pswf_x[0]=0.0f;
  for (ci=1; ci<Np; ci++) {
      pswf_x[ci]=pswf_x[ci-1]+pswf_delta;
  }
  pswf_nn(pswf_c,Np,pswf_x,pswf_y);
  /* multiplied by the sinc() of pixel width */
  for (ci=1; ci<Np; ci++) {
     pswf_y[ci] *=sinf(pswf_x[ci])/pswf_x[ci]; 
  }
  /* copy this to 0..N-1 positions, split at the middle */
  cj=0;
  for (ci=N0; ci<N; ci++) {
     buff[ci]=pswf_y[cj]; 
     cj++;
  }
  for (ci=0; ci<N0; ci++) {
     buff[ci]=pswf_y[Np-ci-1];
  }
  /* now copy this column to the remaining N-1 columns of the buffer */
  for (ci=1; ci<N; ci++) {
   memcpy(&buff[ci*N],&buff[0],sizeof(float)*N);
  }

  free(pswf_x);
  free(pswf_y);

  /* create a backup of 1st column */
  float *y0,*d0;
  if ((y0=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  memcpy(y0,&buff[0],sizeof(float)*N);
  if ((d0=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* scale each row */
  for (ci=0; ci<N; ci++) {
    my_fcopy(N,&buff[ci],N,d0,1);
    for (cj=0; cj<N; cj++) {
     d0[cj]*=y0[cj];
    }
    my_fcopy(N,d0,1,&buff[ci],N);
  }

  free(y0);
  free(d0);
  return 0;
}
