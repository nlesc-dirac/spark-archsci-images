/*
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

#include "gridder.h"
#include <cufft.h>
#include <cublas_v2.h>

__global__ void 
fftshift_1D(cuFloatComplex *u_d, int N) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if(i < N) {
    //float a = powf(-1.0f,i&1);
    float a = float(1-2*((i)&1));
    u_d[i].x *= a;
    u_d[i].y *= a;
    }
}
#define IDX2R(i,j,N) (((i)*(N))+(j))
__global__ void 
fftshift_2D(cuFloatComplex *data, int N, int N0) {
    int i = threadIdx.y + blockDim.y * blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if (i < N && j < N) {
    //float a = powf(-1.0f, (i+j)&1);
    float a = float(1-2*((N0-i+N0-j)&1));
    data[IDX2R(i,j,N)].x *= a;
    data[IDX2R(i,j,N)].y *= a;
    }
}

__global__ void
kernel_lmpswf(cuFloatComplex *din, double *lmgrid, double *denom, float *pswfxy, float w, int Np, int Nz0, int Npad) {
    int i = threadIdx.y + blockDim.y * blockIdx.y;
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if (i<Np && j<Np) {
       float phaseterm=float(lmgrid[i*Np+j])*(w);
//printf("i=%d j=%d Np=%d Nz0=%d Npad=%d w=%f phase=%f\n",i,j,Np,Nz0,Npad,w,phaseterm);
       float cosp,sinp;
       sincosf(phaseterm,&sinp,&cosp);
       float invdenom=fdividef(pswfxy[(i+Npad)*Nz0+Npad+j],(float)denom[i*Np+j]);
//printf("i=%d j=%d Np=%d Nz0=%d Npad=%d\n",i,j,Np,Nz0,Npad);
       din[(i+Npad)*Nz0+Npad+j].x=cosp*invdenom;
       din[(i+Npad)*Nz0+Npad+j].y=sinp*invdenom;
//__syncthreads();
    }
}

extern "C" {

static void
checkCudaError(cudaError_t err, const char *file, int line) {
    if(!err)
        return;
    fprintf(stderr,"GPU (CUDA): %s %s %d\n", cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
}


/*  pswfx: Nz0xNz0
    lmgrid, denom : NpdxNpd
    wkernel: complex float -> cast as float, only copy NpxNp values
*/
void
evaluate_wplane_fft(int Nz0, int Np, int Npd, int Npad, int Npad1, float *pswfxy, double *lmgrid, double *denom, float *wkernel, float *wparr, int Nw, float *peakval) {
  cudaError_t err;
  cufftResult cffterr;
  cublasStatus_t cbstatus;
  cublasHandle_t cbhandle;

//printf("Nz0=%d Np=%d Npd=%d Npad=%d Npad1=%d\n",Nz0,Np,Npd,Npad,Npad1);
  float *dpswfxy;
  double *dlmgrid,*ddenom;
  cuFloatComplex *ddin;
  cufftHandle plan0;
  cudaSetDevice(0);
  cbstatus=cublasCreate(&cbhandle);
  if (cbstatus!=CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr,"%s: %d: CUBLAS create fail\n",__FILE__,__LINE__);
    exit(1);
  }


  cffterr=cufftPlan2d(&plan0, Nz0, Nz0, CUFFT_C2C);
  if (cffterr!=CUFFT_SUCCESS) {
    fprintf(stderr,"%s: %d: CUFFT error\n",__FILE__,__LINE__);
    exit(1);
  }

  err=cudaMalloc((void**)&dpswfxy, sizeof(float)*(Nz0*Nz0));
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMalloc((void**)&dlmgrid, sizeof(double)*(Npd*Npd));
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMalloc((void**)&ddenom, sizeof(double)*(Npd*Npd));
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMalloc((void**)&ddin, sizeof(cuFloatComplex)*(Nz0*Nz0));
  checkCudaError(err,__FILE__,__LINE__);

  err=cudaMemcpy(dpswfxy,pswfxy,sizeof(float)*Nz0*Nz0,cudaMemcpyHostToDevice);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMemcpy(dlmgrid,lmgrid,sizeof(double)*Npd*Npd,cudaMemcpyHostToDevice);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMemcpy(ddenom,denom,sizeof(double)*Npd*Npd,cudaMemcpyHostToDevice);
  checkCudaError(err,__FILE__,__LINE__);

  /* scale by -2pi */
  double alpha=-2.0*M_PI;
  cbstatus=cublasDscal(cbhandle, Npd*Npd, &alpha, dlmgrid, 1);

  dim3 threadsPerBlock(16, 16);
  dim3 numBlocks((Nz0+threadsPerBlock.x-1)/threadsPerBlock.x,
               (Nz0+threadsPerBlock.y-1)/threadsPerBlock.y);
  dim3 numBlocks1((Npd+threadsPerBlock.x-1)/threadsPerBlock.x,
               (Npd+threadsPerBlock.y-1)/threadsPerBlock.y);

  cudaMemset(ddin, 0, sizeof(cuFloatComplex)*Nz0*Nz0);

  err=cudaGetLastError(); /* reset errors */
  for (int nw=0; nw<Nw; nw++) {
   kernel_lmpswf<<<numBlocks1,threadsPerBlock>>>(ddin, dlmgrid, ddenom, dpswfxy, wparr[nw], Npd, Nz0, Npad);
   cudaDeviceSynchronize();
   err=cudaGetLastError();
   checkCudaError(err,__FILE__,__LINE__);


   fftshift_2D<<<numBlocks,threadsPerBlock>>>(ddin, Nz0, Npad+Npd/2);
   cudaDeviceSynchronize();
   err=cudaGetLastError();
   checkCudaError(err,__FILE__,__LINE__);

   cffterr=cufftExecC2C(plan0, (cufftComplex*)ddin, (cufftComplex *)ddin, CUFFT_INVERSE);
   if (cffterr!=CUFFT_SUCCESS) {
     fprintf(stderr,"%s: %d: CUFFT error\n",__FILE__,__LINE__);
     exit(1);
   }

   fftshift_2D<<<numBlocks,threadsPerBlock>>>(ddin, Nz0, Npad+Npd/2);
   cudaDeviceSynchronize();
   err=cudaGetLastError();
   checkCudaError(err,__FILE__,__LINE__);

   for (int nrow=Npad+Npad1; nrow<Npad+Npad1+Np; nrow++) {
    /* size is 2 times, because we copy float */
    err=cudaMemcpy(&wkernel[nw*Np*Np*2+(nrow-Npad-Npad1)*Np*2],&ddin[nrow*Nz0+Npad+Npad1],sizeof(float)*Np*2,cudaMemcpyDeviceToHost);
    checkCudaError(err,__FILE__,__LINE__);
   }
   cudaDeviceSynchronize();

   cudaMemset(ddin, 0, sizeof(cuFloatComplex)*Nz0*Nz0);
  }

  /* also calculate kernel for w=0, to normalize */
  kernel_lmpswf<<<numBlocks1,threadsPerBlock>>>(ddin, dlmgrid, ddenom, dpswfxy, 0.0f, Npd, Nz0, Npad);
  err=cudaGetLastError();
  checkCudaError(err,__FILE__,__LINE__);


  fftshift_2D<<<numBlocks,threadsPerBlock>>>(ddin, Nz0, Npad+Npd/2);
  err=cudaGetLastError();
  checkCudaError(err,__FILE__,__LINE__);

  cffterr=cufftExecC2C(plan0, (cufftComplex*)ddin, (cufftComplex *)ddin, CUFFT_INVERSE);
  if (cffterr!=CUFFT_SUCCESS) {
     fprintf(stderr,"%s: %d: CUFFT error\n",__FILE__,__LINE__);
     exit(1);
  }

  cuFloatComplex w0;
  err=cudaMemcpy(&w0,&ddin[(Nz0/2)*Nz0+Nz0/2],sizeof(cuFloatComplex),cudaMemcpyDeviceToHost);
  checkCudaError(err,__FILE__,__LINE__);

  *peakval=cuCabsf(w0);

  cufftDestroy(plan0);
  cudaFree(dpswfxy);
  cudaFree(dlmgrid);
  cudaFree(ddenom);
  cudaFree(ddin);

  cbstatus=cublasDestroy(cbhandle);
  if (cbstatus!=CUBLAS_STATUS_SUCCESS) {
    fprintf(stderr,"%s: %d: CUBLAS create fail\n",__FILE__,__LINE__);
    exit(1);
  }

}

}
