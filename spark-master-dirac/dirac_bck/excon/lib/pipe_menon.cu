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
#include <glib.h>

//declare texture reference
texture<float,cudaTextureType2D,cudaReadModeElementType> texreference;

/* sum up all N elements of vector input 
 and save (per block) in output (size > number of blocks) */
__global__ void
plus_reduce_multi(float *input, int N, int blockDim_2, float *output) {
 // Each block loads its elements into shared memory
 extern __shared__ float x[];
 int tid = threadIdx.x;
 int i = blockIdx.x*blockDim.x + threadIdx.x;
 x[tid] = (i<N) ? input[i] : 0.0f; // last block may pad with 0’s
 __syncthreads();
 // Build summation tree over elements, handling case where B is not a power of two.
  int nTotalThreads = blockDim_2; // Total number of threads, rounded up to the next power of two
  while(nTotalThreads > 1) {
   int halfPoint = (nTotalThreads >> 1); // divide by two
    if (tid < halfPoint) {
     int thread2 = tid + halfPoint;
     if (thread2 < blockDim.x) { // Skipping the fictitious threads blockDim.x ... blockDim_2-1
      x[tid] = x[tid]+x[thread2];
     }
    }
    __syncthreads();
    nTotalThreads = halfPoint; // Reducing the binary tree size by two
 }

 /* add back to total */
 if( tid == 0 ) {
  output[blockIdx.x]=x[tid];
 }
}


/* sum up all N elements of vector input 
 NOTE: only 1 block should be used */
__global__ void
plus_reduce(float *input, int N, int blockDim_2, float *total) {
 // Each block loads its elements into shared memory
 extern __shared__ float x[];
 int tid = threadIdx.x;
 int i = blockIdx.x*blockDim.x + threadIdx.x;
 x[tid] = (i<N) ? input[i] : 0.0f; // last block may pad with 0’s
 __syncthreads();
 // Build summation tree over elements, handling case where B is not a power of two.
  int nTotalThreads = blockDim_2; // Total number of threads, rounded up to the next power of two
  while(nTotalThreads > 1) {
   int halfPoint = (nTotalThreads >> 1); // divide by two
    if (tid < halfPoint) {
     int thread2 = tid + halfPoint;
     if (thread2 < blockDim.x) { // Skipping the fictitious threads blockDim.x ... blockDim_2-1
      x[tid] = x[tid]+x[thread2];
     }
    }
    __syncthreads();
    nTotalThreads = halfPoint; // Reducing the binary tree size by two
 }

 /* add back to total */
 if( tid == 0 ) {
  *total=*total+x[tid];
 }
}


__global__ void
kernel_ncpweight(float uf,float vf,float *wtd, float uvscale) {
 __shared__ float a[6];
 __shared__ float b[6];
 __shared__ float c[6];
 __shared__ float x;
 int tid = threadIdx.x;
 if (tid==0) {
  a[0] =0.2589f;
  b[0] =109.4f;
  c[0] =13.09f;
  a[1] =0.6783f;
  b[1] =88.86f;
  c[1] =35.7f;
  a[2] =0.0868f;
  b[2] =212.1f;
  c[2] =10.6f;
  a[3] =-0.5993f;
  b[3] =300.0f;
  c[3] =84.17f;
  a[4] =1.476e+05f;
  b[4] =-4327.0f;
  c[4] =1391.0f;
  a[5] =-5.714f;
  b[5] =13.01f;
  c[5] =185.3f;

  x=sqrtf(uf*uf+vf*vf)/uvscale; /* scale by inverse scale get x in [0,800] */
 }
 __syncthreads();
 if (tid<6) {
  if (x<25.0f||x>900.0f) {
   wtd[tid]=0.0f;
  } else if (x<65.0f) {
   wtd[tid]=2.0517f/(1.0f+expf(-(x-40.0f)*0.333333333f))/6.0f;
  } else if (x>800.f) {
   float x2=(x-800.0f);
   wtd[tid]=0.1832f*expf(-x2*x2*0.001f)/6.0f;
  } else {
   float t=(x-b[tid])/c[tid];
   wtd[tid]=expf(-t*t)*a[tid];
  }
 }
}

__global__ void
kernel_pmconvolution(int N, float *ud, float *vd, float *wtd, float uf, float vf, float *ed) {
  // Each thread saves error into shared memory
  extern __shared__ float ek[];
  int ui=threadIdx.x + blockIdx.x * blockDim.x;
  int tid=threadIdx.x;
  ek[tid]=0.0f;
  if (ui<N) {
   float x=-ud[ui]+uf+0.5f;
   float y=vd[ui]+vf+0.5f;
   /* x,y in [0,1], take absolute value of kernel */
   float tt=fabsf(tex2D(texreference,x,y));
   ek[tid]=wtd[ui]*tt;
//printf("uf,vf %f,%f ud,vd %f,%f x=%f y=%f tt=%f e=%f\n",uf,vf,ud[ui],vd[ui],x,y,tt,ek[tid]);
  }
  __syncthreads();
  // Build summation tree over elements, assuming blockDim.x is power of 2.
  for(int s=blockDim.x/2; s>0; s=s/2) {
    if(tid < s) ek[tid] += ek[tid + s];
   __syncthreads();
  }

  /* copy back the sum to proper location in ed */
  if(tid==0) {
   ed[blockIdx.x]=ek[0];
  }

}

extern "C" {

static void
checkCudaError(cudaError_t err, const char *file, int line)
{
    if(!err)
        return;
    fprintf(stderr,"GPU (CUDA): %s %s %d\n", cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
}

/* need power of 2 for tree reduction to work */
static int
NearestPowerOf2 (int n){
  if (!n) return n;  //(0 == 2^0)

  int x = 1;
  while(x < n) {
      x <<= 1;
  }
  return x;
}

/* u,v,wt: Nx1 arrays */
static float 
cudakernel_pmconvolution(int card, int N,float *u,float *v,float *wt,float uf,float vf, float uvscale, int convmode) {
  cudaError_t err;

  float *ud,*vd,*wtd,*ed;
  float *totald;
  cudaSetDevice(card);
  int threadsPerBlock=128;
  int BlocksPerGrid=(N+threadsPerBlock-1)/threadsPerBlock;
  if (BlocksPerGrid==0) { /* catch situation when N=1 */
   BlocksPerGrid=1;
  } 

  err=cudaMalloc((void**)&ud, sizeof(float)*(N));
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMalloc((void**)&vd, sizeof(float)*(N));
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMalloc((void**)&wtd, sizeof(float)*((N>6?N:6))); /* note: make sure length>6, because we use it for temp storage */
  checkCudaError(err,__FILE__,__LINE__);
  /* to store sum of each block */
  err=cudaMalloc((void**)&ed, sizeof(float)*(BlocksPerGrid));
  checkCudaError(err,__FILE__,__LINE__);

  err=cudaMemcpy(ud, (void*)u, sizeof(float)*(N), cudaMemcpyHostToDevice);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMemcpy(vd, (void*)v, sizeof(float)*(N), cudaMemcpyHostToDevice);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaMemcpy(wtd, (void*)wt, sizeof(float)*(N), cudaMemcpyHostToDevice); 
  checkCudaError(err,__FILE__,__LINE__);

  //printf("card=%d N=%d threads=%d blocks=%d\n",card,N,threadsPerBlock,BlocksPerGrid);
  kernel_pmconvolution<<<BlocksPerGrid,threadsPerBlock,sizeof(float)*threadsPerBlock>>>(N,ud,vd,wtd,uf,vf,ed);
  cudaDeviceSynchronize();

  err = cudaGetLastError();
  checkCudaError(err,__FILE__,__LINE__);

  err=cudaMalloc((void**)&totald, sizeof(float));
  checkCudaError(err,__FILE__,__LINE__);
  cudaMemset(totald, 0, sizeof(float));

  /* summation over ed */
  if (BlocksPerGrid<threadsPerBlock) {
    /* one kernel launch is enough */
    plus_reduce<<< 1, BlocksPerGrid, sizeof(float)*BlocksPerGrid>>>(ed, BlocksPerGrid, NearestPowerOf2(BlocksPerGrid), totald);
    cudaDeviceSynchronize();
  } else {
    /* multiple kernel launches */
    int L=(BlocksPerGrid+threadsPerBlock-1)/threadsPerBlock;
    /* reuse wtd as temp storage */
    plus_reduce_multi<<< L, threadsPerBlock, sizeof(float)*threadsPerBlock>>>(ed, BlocksPerGrid, NearestPowerOf2(threadsPerBlock), wtd);
    cudaDeviceSynchronize();
    plus_reduce<<< 1, L, sizeof(float)*L>>>(wtd, L, NearestPowerOf2(L), totald);
    cudaDeviceSynchronize();
  }
  err = cudaGetLastError();
  checkCudaError(err,__FILE__,__LINE__);
  float total;
  err=cudaMemcpy(&total,totald,sizeof(float),cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  /* now total has (W_k \odot C_k) */
  /* calculate ncp_weight/(W_k \odot C_k): 6 threads for 6 order poly */
  /* reuse wtd as temp storage */
  float ncpwt;
  if (convmode==CONV_MODE_NCP) {
  cudaMemset(totald, 0, sizeof(float));
  kernel_ncpweight<<<1,6>>>(uf,vf,wtd,uvscale);
  cudaDeviceSynchronize();
  plus_reduce<<< 1, 6, sizeof(float)*6>>>(wtd, 6, NearestPowerOf2(6), totald);
  cudaDeviceSynchronize();
  err=cudaMemcpy(&ncpwt,totald,sizeof(float),cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  } else {
   ncpwt=1.0f; /* uniform weight */
  }

  err=cudaFree(ud);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaFree(vd);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaFree(wtd);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaFree(ed);
  checkCudaError(err,__FILE__,__LINE__);
  err=cudaFree(totald);
  checkCudaError(err,__FILE__,__LINE__);
//printf("ncpwt %f conv %f\n",ncpwt,total);
  return ncpwt/(total+1e-12f);
}

/* function to set up a GPU, should be called only once */
static void
attach_gpu_to_thread(int card, float *wkernel, int Np, void **carrayp) {
 cudaError_t err;

 cudaChannelFormatDesc channel;
 cudaArray* carray;
 cudaSetDevice(card);
 
 //create channel to describe data type
 channel=cudaCreateChannelDesc<float>();

 //allocate device memory for cuda array
 err=cudaMallocArray(&carray,&channel,Np,Np);
 checkCudaError(err,__FILE__,__LINE__);

 err=cudaMemcpyToArray(carray,0,0,wkernel,sizeof(float)*Np*Np,cudaMemcpyHostToDevice);
 checkCudaError(err,__FILE__,__LINE__);


 //all coordinate axes are in [0,1]
 texreference.normalized=true;
 //set texture filter mode property
 //use cudaFilterModePoint or cudaFilterModeLinear
 texreference.filterMode=cudaFilterModeLinear;
 //set texture address mode property
 //use cudaAddressModeClamp or cudaAddressModeWrap
 texreference.addressMode[0]=cudaAddressModeClamp;
 texreference.addressMode[1]=cudaAddressModeClamp;

 //bind texture reference with cuda array
 err=cudaBindTextureToArray(texreference,carray);
 checkCudaError(err,__FILE__,__LINE__);

 *carrayp=carray;
}

static void
detach_gpu_from_thread(int card, void *carray) {
 cudaSetDevice(card);
 cudaArray *cudaarr=(cudaArray*)carray;
 //unbind texture reference to free resource
 cudaUnbindTexture(texreference);
 cudaFreeArray(cudaarr);
}


/* slave thread function */
static void *
pipeline_pm_slave_code(void *data)
{
 slave_pmtdata *td=(slave_pmtdata*)data;
 gbpmgdata *gd=(gbpmgdata*)(td->pline->data);
 int tid=td->tid;
 unsigned long int i;
 uvlist *ll;
 int card; /* which GPU */
#ifndef ONE_GPU
 card=tid%2;
#endif
#ifdef ONE_GPU
 card=0;
#endif
 while(1) {
  sync_barrier(&(td->pline->gate1)); /* stop at gate 1*/
  if(td->pline->terminate) break; /* if flag is set, break loop */
  sync_barrier(&(td->pline->gate2)); /* stop at gate 2 */
  if (gd->status[tid]==PT_DO_WORK_GRID) {
/************************* work *********************/
    //printf("thread %d from [%ld,%ld] %ld rows\n",tid,gd->startrow[tid],gd->startrow[tid]+gd->Nrows[tid]-1,gd->Nrows[tid]);
     float *ubuff,*vbuff,*wtbuff;
     int bufsz=1024;
     int bfilled=0;
     if ((ubuff=(float*)malloc(sizeof(float)*(size_t)bufsz))==0) {
        fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
        exit(1);
     }
     if ((vbuff=(float*)malloc(sizeof(float)*(size_t)bufsz))==0) {
        fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
        exit(1);
     }
     if ((wtbuff=(float*)malloc(sizeof(float)*(size_t)bufsz))==0) {
        fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
        exit(1);
     }


     for (i=gd->startrow[tid]; i<gd->startrow[tid]+gd->Nrows[tid]; i++) {
       if (!gd->darr[i].flag) {
//printf("i=%ld u,v %f,%f wold=%f\n",i,gd->darr[i].u,gd->darr[i].v,gd->wold[i]);
         /* calculate right bucket index */
          float tempu = gd->darr[i].u;
          float tempv = gd->darr[i].v;
         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*gd->uvscale;
          float vf=tempv*gd->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)gd->Nx); /*  width */
          float vi=(vf*(float)gd->Ny); /*  width */
          float poffX=-1.0f;
          float poffY=-1.0f;
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*gd->Nx+poffX);
          int y=(int)round(-vi+0.5f*gd->Ny+poffY);
         if (x>=0 && x<gd->Nx && y>=0 && y<gd->Ny) {
           int xlow=(x-gd->M>=0?x-gd->M:0);
           int xhigh=(x+gd->M>gd->Nx?gd->Nx:x+gd->M);
           int ylow=(y-gd->M>=0?y-gd->M:0);
           int yhigh=(y+gd->M>gd->Ny?gd->Ny:y+gd->M);
           bfilled=0;
           for (int bx=xlow; bx<xhigh; bx++ )
           for (int by=ylow; by<yhigh; by++ ) {
//printf("U,V %f,%f (u,v) (%f,%f) (x,y) (%d,%d) -> (bx,by) (%d,%d) [xlow,xhigh] (%d,%d) [ylow,yhigh] (%d,%d)\n",tempu,tempv,uf,vf,x,y,bx,by,xlow,xhigh,ylow,yhigh);
            unsigned long int bi=bx*gd->B+by;
            pthread_mutex_lock(gd->writelock_hash);
            ll=(uvlist*)g_hash_table_lookup(gd->ht,&bi);
            pthread_mutex_unlock(gd->writelock_hash);
            if (ll) { /* found neighbour pixel list */
              /* realloc memory if needed */
              if (bfilled+ll->P>bufsz) {
                if((ubuff=(float*)realloc((void*)ubuff,(size_t)(bfilled+ll->P)*sizeof(float)))==0){
                 fprintf(stderr, "%s: %d: no free memory\n", __FILE__,__LINE__);
                 exit(1);
                }
                if((vbuff=(float*)realloc((void*)vbuff,(size_t)(bfilled+ll->P)*sizeof(float)))==0){
                 fprintf(stderr, "%s: %d: no free memory\n", __FILE__,__LINE__);
                 exit(1);
                }
                if((wtbuff=(float*)realloc((void*)wtbuff,(size_t)(bfilled+ll->P)*sizeof(float)))==0){
                 fprintf(stderr, "%s: %d: no free memory\n", __FILE__,__LINE__);
                 exit(1);
                }

                bufsz=bfilled+ll->P;
              }
              /* copy memory */
              memcpy((void*)&ubuff[bfilled],(void*)ll->u,(size_t)(ll->P)*sizeof(float));
              memcpy((void*)&vbuff[bfilled],(void*)ll->v,(size_t)(ll->P)*sizeof(float));
              for (int cw=0; cw<ll->P; cw++) {
                wtbuff[bfilled+cw]=gd->wold[ll->id[cw]];
              }
              bfilled+=ll->P;
            }
           }
           /* now call cuda kernel */
           /* W_k+1 <= (W_k x G_k) / (W_k \odot C_k) */
           /* G_k = 1 for uniform weights, use NCP_WEIGHT function */
           float ratio=cudakernel_pmconvolution(card,bfilled,ubuff,vbuff,wtbuff,uf,vf,gd->uvscale,gd->convmode);
//printf("%ld ratio=%f\n",i,ratio);
//printf("(u,v) (%f,%f) (uf,vf) (%f,%f) (x,y) (%d,%d) [xlow,xhigh] (%d,%d) [ylow,yhigh] (%d,%d)\n",tempu,tempv,uf,vf,x,y,xlow,xhigh,ylow,yhigh);
           gd->wnew[i]=gd->wold[i]*ratio;
//printf("i=%ld old=%f new=%f\n",i,gd->wold[i],gd->wnew[i]);
         }
       }
     }
   
     free(ubuff);
     free(vbuff);
     free(wtbuff);
/************************* work *********************/
  } else if (gd->status[tid]==PT_DO_AGPU) {
  //printf("thread %d : pix %d\n",tid,gd->Np);
   /* FIXME: also copy wparr (w coords) to GPU for searching */
   attach_gpu_to_thread(card,gd->wkernel,gd->Np,&gd->carray[tid]);
  } else if (gd->status[tid]==PT_DO_DGPU) {
   detach_gpu_from_thread(card,gd->carray[tid]);
  }
 }
 return NULL;
}



/* initialize the pipeline
  and start the slaves rolling 
  create 2N slave threads */
void
init_pm_pipeline(th_pmpipeline *pline,
     void *data, int N)
{
 if ((pline->sd=(slave_pmtdata*)malloc(sizeof(slave_pmtdata)*2*N))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }

 pthread_attr_init(&(pline->attr));
 pthread_attr_setdetachstate(&(pline->attr),PTHREAD_CREATE_JOINABLE);

 init_th_barrier(&(pline->gate1),2*N+1); /* 2N+1 threads, including master */
 init_th_barrier(&(pline->gate2),2*N+1); /* 2N+1 threads, including master */
 pline->terminate=0;
 pline->data=data; /* data should have pointers to t1 and t2 */
 int ci;
 for (ci=0; ci<2*N; ci++) {
  pline->sd[ci].pline=(th_pmpipeline*)pline;
  /* link back t1, t2 to data so they could be freed */
  pline->sd[ci].tid=ci;
 }

 if ((pline->slave=(pthread_t*)malloc(sizeof(pthread_t)*2*N))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }

 for (ci=0; ci<2*N; ci++) {
  pthread_create(&(pline->slave[ci]),&(pline->attr),pipeline_pm_slave_code,(void*)&pline->sd[ci]);
 }
}

/* destroy the pipeline */
/* need to kill the slaves first */
void
destroy_pm_pipeline(th_pmpipeline *pline, int N)
{
 pline->terminate=1;
 sync_barrier(&(pline->gate1));
 int ci;
 for (ci=0; ci<2*N; ci++) {
  pthread_join(pline->slave[ci],NULL);
 }
 free(pline->slave);
 destroy_th_barrier(&(pline->gate1));
 destroy_th_barrier(&(pline->gate2));
 pthread_attr_destroy(&(pline->attr));
 free(pline->sd);
 pline->data=NULL;
}

}
