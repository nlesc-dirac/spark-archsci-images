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

//declare texture reference
texture<cuFloatComplex,cudaTextureType3D,cudaReadModeElementType> texreference;
//texture<cuFloatComplex,cudaTextureType2DLayered,cudaReadModeElementType> texreference;


/* search for val in array list of size nx1 */
__global__ static void cuda_search(int *out, float *list, float val, float lambda, int n)
{
    int start = threadIdx.x + blockIdx.x * blockDim.x;
    /* take abs(w) in lambda for wavelength conversion */
    float wl=fabs(val);
    for (int idx = start; idx < n; idx += gridDim.x * blockDim.x) {
        if (list[idx] >= wl) return;
        float next = list[idx + 1];
        if (idx == n-1 || next >= wl) {
            *out =( next == wl ? (idx + 1) : idx);
            return;
        }
     }
}


__device__ static int cuda_ifftshift_index(int N, int x) {
 if (N%2) { /* odd */
  return((x+N/2)%N);
 } else { /* even */
  return((x+N/2+1)%N);
 }
}

/* reset memory to Nx x Ny to eliminate points off grid */
__global__ void cuda_resetmemory(unsigned long int *outoff, int Mx, int Nx, int Ny) {
    int ui = threadIdx.x + blockIdx.x * blockDim.x;
    int vi = threadIdx.y + blockIdx.y * blockDim.y;
    size_t oidx = 2*(vi+ui*Mx);
    if (ui<Mx && vi<Mx) {
     outoff[oidx]=Nx*Ny;
     outoff[oidx+1]=Nx*Ny;
    }
    __syncthreads();
}


/* total support is Mx pixels, so total pixels is Mx x Mx 
 since we also grid the -ve u,v,w point, total data size is 2 Mx x Mx */
__global__ void cuda_griddata(cuFloatComplex* outputimg, cuFloatComplex *outputpsf, unsigned long int *outoff, int Mx, float u, float v, float w, float wt, cuFloatComplex xx, cuFloatComplex yy, cuFloatComplex xy, cuFloatComplex yx, float maxW, float expW, float lambda, float uvscale, float deltaU, int Nx, int Ny, int Nw, int Np, int Nz)
{

    /* note: v axis is negative */
    /* all position calculations are done using double precision */
    __shared__ float zl,zh;
    __shared__ float wh,wl,wa;
    __shared__ double uf,vf,up,vp;
    __shared__ int signofw;
    __shared__ double pix0;
    __shared__ cuFloatComplex sI;
    __shared__ double NxNp;
    if (threadIdx.x == 0) {
      pix0=-sqrt(0.8); //trueimg/zeropaddedimg?? old -1
      sI=cuCmulf(cuCaddf(xx,yy),make_cuFloatComplex(0.5f,0.0f));
      /* u,v: corresponding pixel value */
      uf=double(u)*double(uvscale); /* in [-0.5,0.5] */
      vf=double(v)*double(uvscale); /* in [-0.5,0.5] */
      /* nearest central pixel */
      up=rint(uf*double(Nx)+0.5*double(Nx))-1.0; /* in [0,Nx-1] */
      vp=rint(-vf*double(Ny)+0.5*double(Ny))-1.0; /* in [0,Ny-1] */
 
      float deltaz=fdividef(0.1f,(float)Nw); /* 0.1/Nw */
      wa=fdividef(fabsf(w),maxW);
      /* w in normalized coords [0,1], with sqrt() spacing */
      //float z=sqrtf(wa); /* to make general, use power, exponent __powf(x,y) */
      float z=powf(wa,expW);
      zl=z-deltaz;
      zl=(zl<0.0f?0.0f:zl);
      zh=z+deltaz;
      zh=(zh>1.0f?1.0f:zh);
      //wl=zl*zl;
      //wh=zh*zh;
      wl=powf(zl,1.0f/expW);
      wh=powf(zh,1.0f/expW);
      /* sign of the w term */
      signofw=signbit(w);
      NxNp=1.0/(double(Nx*Np)*double(deltaU));
    }
    __syncthreads();
    /* ui,vi: in [0,Mx-1] */
    int ui = threadIdx.x + blockIdx.x * blockDim.x;
    int vi = threadIdx.y + blockIdx.y * blockDim.y;
    /* absolute pixel value for this thread in [0,Nx-1], [0,Ny-1] */
    int x0 = int(up)-(Mx>>1)+ui;
    int y0 = int(vp)-(Mx>>1)+vi;
    int x1 = cuda_ifftshift_index(Nx,x0);
    int y1 = cuda_ifftshift_index(Ny,y0);

    /* normalized distance of this pixel from central pixel : in images pixels of [0,Nx-1] => scale to [0,Nz-1] and note that Np pixels map to [0,1]  */
    float x = float(((uf+0.5)*double(Nx)-double(x0)+pix0)*NxNp+0.5);
    float y = float(((-vf+0.5)*double(Ny)-double(y0)+pix0)*NxNp+0.5);

    size_t oidx = 2*(vi+ui*Mx);
    /* apart from checking the support, also throw out -ve pixels */
    if (ui<Mx && vi<Mx && x0>=0 && x0<Nx && y0>=0 && y0<Ny) {
     cuFloatComplex psfl=tex3D(texreference, x, y, zl);
     cuFloatComplex psfh=tex3D(texreference, x, y, zh);
     cuFloatComplex psf;
     psf.x=fdividef(psfl.x*(wh-wa)+psfh.x*(wa-wl),wh-wl);
     psf.y=fdividef(psfl.y*(wh-wa)+psfh.y*(wa-wl),wh-wl);
     /* multiply pswf with inverse sigma */
     psf=cuCmulf(make_cuFloatComplex(wt,0.0f),psf);
#ifdef ONE_GPU
     //printf("Np x delu=%f, u,v=%f,%f, scale=%f, deltau=%f, uf,vf=%f,%f up,vp=%f,%f ui,vi=(%d,%d) x0,y0=(%d,%d) x1,y1=(%d,%d) (x,y,z)=%f,%f,%f pswf=%f,%f\n",float(Np)*deltaU,u,v,uvscale,deltaU,uf,vf,up,vp,ui,vi,x0,y0,x1,y1,x,y,z,psf.x,psf.y);
     //printf("%f %f %f\n",x,y,z);
#endif
     if (signofw) { /* w is negative */
        psf=cuConjf(psf);
     }
     outputimg[oidx]=cuCmulf(sI,psf);
     outputpsf[oidx]=psf;
     outoff[oidx]=x1*Ny+y1;
    }
    __syncthreads();

    /* now handle -ve u,v,w point */
    if (threadIdx.x == 0) {
      /* flip sign of u,v */
      uf=-uf; /* in [-0.5,0.5] */
      vf=-vf; /* in [-0.5,0.5] */
      /* nearest central pixel */
      up=-rint(-(uf*double(Nx)+0.5*double(Nx)))-1.0; /* in [0,Nx-1] */
      vp=-rint(-(-vf*double(Ny)+0.5*double(Ny)))-1.0; /* in [0,Ny-1] */
    }
    __syncthreads();

    /* absolute pixel value for this thread in [0,Nx-1], [0,Ny-1] */
    x0 = int(up)-(Mx>>1)+ui;
    y0 = int(vp)-(Mx>>1)+vi;
    x1 = cuda_ifftshift_index(Nx,x0);
    y1 = cuda_ifftshift_index(Ny,y0);

    x = float(((uf+0.5)*double(Nx)-double(x0)+pix0)*NxNp+0.5);
    y = float(((-vf+0.5)*double(Ny)-double(y0)+pix0)*NxNp+0.5);

    oidx = 2*(vi+ui*Mx)+1;
    /* apart from checking the support, also throw out -ve pixels */
    if (ui<Mx && vi<Mx && x0>=0 && x0<Nx && y0>=0 && y0<Ny) {
     cuFloatComplex psfl=tex3D(texreference, x, y, zl);
     cuFloatComplex psfh=tex3D(texreference, x, y, zh);
     cuFloatComplex psf;
     psf.x=fdividef(psfl.x*(wh-wa)+psfh.x*(wa-wl),wh-wl);
     psf.y=fdividef(psfl.y*(wh-wa)+psfh.y*(wa-wl),wh-wl);
     /* multiply pswf with inverse sigma */
     psf=cuCmulf(make_cuFloatComplex(wt,0.0f),psf);
#ifdef ONE_GPU
     //printf("Np x delu=%f, u,v=%f,%f, scale=%f, deltau=%f, uf,vf=%f,%f up,vp=%f,%f ui,vi=(%d,%d) x0,y0=(%d,%d) x1,y1=(%d,%d) (x,y,z)=%f,%f,%f pswf=%f,%f\n",float(Np)*deltaU,u,v,uvscale,deltaU,uf,vf,up,vp,ui,vi,x0,y0,x1,y1,x,y,z,psf.x,psf.y);
     //printf("%f %f %f\n",x,y,z);
#endif
     //cuFloatComplex psf=tex2DLayered(texreference, x, y, z);
     if (!signofw) { /* -w is negative */
        psf=cuConjf(psf);
     }
     outputimg[oidx]=cuCmulf(cuConjf(sI),psf);
     outputpsf[oidx]=psf;
     outoff[oidx]=x1*Ny+y1;
    }
    __syncthreads();
}

/* total support is Mx pixels, so total pixels is Mx x Mx 
 since we also grid the -ve u,v,w point, total data size is 2 Mx x Mx */
__global__ void cuda_griddata_iquv(cuFloatComplex* outputimg, cuFloatComplex *outputpsf, unsigned long int *outoff, cuFloatComplex* outputimgQ, cuFloatComplex* outputimgU, cuFloatComplex* outputimgV, int Mx, float u, float v, float w, float wt, cuFloatComplex xx, cuFloatComplex yy, cuFloatComplex xy, cuFloatComplex yx, float maxW, float expW, float lambda, float uvscale, float deltaU, int Nx, int Ny, int Nw, int Np, int Nz)
{

    /* note: v axis is negative */
    __shared__ float zl,zh;
    __shared__ float wh,wl,wa;
    __shared__ double uf,vf,up,vp;
    __shared__ int signofw;
    __shared__ float pix0;
    __shared__ cuFloatComplex sI,sQ,sU,sV;
    __shared__ double NxNp;
    if (threadIdx.x == 0) {
      pix0=-sqrt(0.8); //trueimg/zeropaddedimg?? old -1
      sI=cuCmulf(cuCaddf(xx,yy),make_cuFloatComplex(0.5f,0.0f));
      /* Q =(XX-YY)/2 U=(XY+YX)/2 V=imag(YX-XY)/2 */
      sQ=cuCmulf(cuCsubf(xx,yy),make_cuFloatComplex(0.5f,0.0f));
      sU=cuCmulf(cuCaddf(xy,yx),make_cuFloatComplex(0.5f,0.0f));
      sV=cuCmulf(cuCsubf(yx,xy),make_cuFloatComplex(0.0f,0.5f));
      /* u,v: corresponding pixel value */
      uf=double(u)*double(uvscale); /* in [-0.5,0.5] */
      vf=double(v)*double(uvscale); /* in [-0.5,0.5] */
      /* nearest central pixel pixel */
      up=rint(uf*double(Nx)+0.5*double(Nx))-1.0; /* in [0,Nx-1] */
      vp=rint(-vf*double(Ny)+0.5*double(Ny))-1.0; /* in [0,Ny-1] */
 
      float deltaz=fdividef(0.1f,(float)Nw); /* 0.1/Nw */
      wa=fdividef(fabsf(w),maxW);
      /* w in normalized coords [0,1], with w^expW spacing */
      //float z=sqrtf(wa);
      float z=powf(wa,expW);
      zl=z-deltaz;
      zl=(zl<0.0f?0.0f:zl);
      zh=z+deltaz;
      zh=(zh>1.0f?1.0f:zh);
      //wl=zl*zl;
      //wh=zh*zh;
      wl=powf(zl,1.0f/expW);
      wh=powf(zh,1.0f/expW);
      /* sign of the w term */
      signofw=signbit(w);
      NxNp=1.0/(double(Nx*Np)*double(deltaU));
    }
    __syncthreads();
    /* ui,vi: in [0,Mx-1] */
    int ui = threadIdx.x + blockIdx.x * blockDim.x;
    int vi = threadIdx.y + blockIdx.y * blockDim.y;
    /* absolute pixel value for this thread in [0,Nx-1], [0,Ny-1] */
    int x0 = int(up)-(Mx>>1)+ui;
    int y0 = int(vp)-(Mx>>1)+vi;
    int x1 = cuda_ifftshift_index(Nx,x0);
    int y1 = cuda_ifftshift_index(Ny,y0);

    /* normalized distance of this pixel from central pixel : in images pixels of [0,Nx-1] => scale to [0,Nz-1] and note that Np pixels map to [0,1]  */
    float x = float(((uf+0.5)*double(Nx)-double(x0)+pix0)*NxNp+0.5);
    float y = float(((-vf+0.5)*double(Ny)-double(y0)+pix0)*NxNp+0.5);

    size_t oidx = 2*(vi+ui*Mx);
    /* apart from checking the support, also throw out -ve pixels */
    if (ui<Mx && vi<Mx && x0>=0 && x0<Nx && y0>=0 && y0<Ny) {
     cuFloatComplex psfl=tex3D(texreference, x, y, zl);
     cuFloatComplex psfh=tex3D(texreference, x, y, zh);
     cuFloatComplex psf;
     psf.x=fdividef(psfl.x*(wh-wa)+psfh.x*(wa-wl),wh-wl);
     psf.y=fdividef(psfl.y*(wh-wa)+psfh.y*(wa-wl),wh-wl);
     /* multiply pswf with inverse sigma */
     psf=cuCmulf(make_cuFloatComplex(wt,0.0f),psf);
#ifdef ONE_GPU
     //printf("Np x delu=%f, u,v=%f,%f, scale=%f, deltau=%f, uf,vf=%f,%f up,vp=%f,%f ui,vi=(%d,%d) x0,y0=(%d,%d) x1,y1=(%d,%d) (x,y,z)=%f,%f,%f pswf=%f,%f\n",float(Np)*deltaU,u,v,uvscale,deltaU,uf,vf,up,vp,ui,vi,x0,y0,x1,y1,x,y,z,psf.x,psf.y);
     //printf("%f %f %f\n",x,y,z);
#endif
     if (signofw) { /* w is negative */
        psf=cuConjf(psf);
     }
     outputpsf[oidx]=psf;
     outputimg[oidx]=cuCmulf(sI,psf);
     outoff[oidx]=x1*Ny+y1;
     outputimgQ[oidx]=cuCmulf(sQ,psf);
     outputimgU[oidx]=cuCmulf(sU,psf);
     outputimgV[oidx]=cuCmulf(sV,psf);
    }
    __syncthreads();

    /* now handle -ve u,v,w point */
    if (threadIdx.x == 0) {
      /* flip sign of u,v */
      uf=-uf; /* in [-0.5,0.5] */
      vf=-vf; /* in [-0.5,0.5] */
      /* nearest central pixel pixel */
      up=-rint(-(uf*double(Nx)+0.5f*double(Nx)))-1.0; /* in [0,Nx-1] */
      vp=-rint(-(-vf*double(Ny)+0.5f*double(Ny)))-1.0; /* in [0,Ny-1] */
    }
    __syncthreads();

    /* absolute pixel value for this thread in [0,Nx-1], [0,Ny-1] */
    x0 = int(up)-(Mx>>1)+ui;
    y0 = int(vp)-(Mx>>1)+vi;
    x1 = cuda_ifftshift_index(Nx,x0);
    y1 = cuda_ifftshift_index(Ny,y0);

    x = float(((uf+0.5)*double(Nx)-double(x0)+pix0)*NxNp+0.5);
    y = float(((-vf+0.5)*double(Ny)-double(y0)+pix0)*NxNp+0.5);

    oidx = 2*(vi+ui*Mx)+1;
    /* apart from checking the support, also throw out -ve pixels */
    if (ui<Mx && vi<Mx && x0>=0 && x0<Nx && y0>=0 && y0<Ny) {
     cuFloatComplex psfl=tex3D(texreference, x, y, zl);
     cuFloatComplex psfh=tex3D(texreference, x, y, zh);
     cuFloatComplex psf;
     psf.x=fdividef(psfl.x*(wh-wa)+psfh.x*(wa-wl),wh-wl);
     psf.y=fdividef(psfl.y*(wh-wa)+psfh.y*(wa-wl),wh-wl);
     /* multiply pswf with inverse sigma */
     psf=cuCmulf(make_cuFloatComplex(wt,0.0f),psf);
#ifdef ONE_GPU
     //printf("Np x delu=%f, u,v=%f,%f, scale=%f, deltau=%f, uf,vf=%f,%f up,vp=%f,%f ui,vi=(%d,%d) x0,y0=(%d,%d) x1,y1=(%d,%d) (x,y,z)=%f,%f,%f pswf=%f,%f\n",float(Np)*deltaU,u,v,uvscale,deltaU,uf,vf,up,vp,ui,vi,x0,y0,x1,y1,x,y,z,psf.x,psf.y);
     //printf("%f %f %f\n",x,y,z);
#endif
     //cuFloatComplex psf=tex2DLayered(texreference, x, y, z);
     if (!signofw) { /* -w is negative */
        psf=cuConjf(psf);
     }
     outputimg[oidx]=cuCmulf(cuConjf(sI),psf);
     outputpsf[oidx]=psf;
     outoff[oidx]=x1*Ny+y1;
     outputimgQ[oidx]=cuCmulf(cuConjf(sQ),psf);
     outputimgU[oidx]=cuCmulf(cuConjf(sU),psf);
     outputimgV[oidx]=cuCmulf(cuConjf(sV),psf);
    }
    __syncthreads();
}


#ifdef ONE_GPU
__global__ void 
cuda_printtexture(int Nx,int Ny,int Nw) {
 /* only one thread does work */
 if((threadIdx.x==0) && (blockIdx.x==0) ) {
  printf("Nx=%d Ny=%d Nw=%d\n",Nx,Ny,Nw);
  float delx=1.0f/float(Nx);
  float dely=1.0f/float(Ny);
  float delz=1.0f/float(Nw);
  __syncthreads();
  for (int ci=0; ci<Nx; ci++) {
   for (int cj=0; cj<Ny; cj++) {
    for (int ck=0; ck<Nw; ck++) {
     float x=float(ci)*delx;
     float y=float(cj)*dely;
     float z=float(ck)*delz;
     cuFloatComplex pswf=tex3D(texreference, y, x, 0.1f);
     //cuFloatComplex pswf=tex2DLayered(texreference, x, y, z);
     printf("%f %f %f %f %f\n",x,y,z,pswf.x,pswf.y);
     __syncthreads();
    }  
   }
  }
 }
}
#endif

extern "C" {

static void
checkCudaError(cudaError_t err, const char *file, int line)
{
    if(!err)
        return;
    fprintf(stderr,"GPU (CUDA): %s %s %d\n", cudaGetErrorString(err),file,line);
    exit(EXIT_FAILURE);
}

#ifdef ONE_GPU
static void
debug_show_texture(int card,int Nx,int Ny,int Nw) {
 printf("card=%d, Nx=%d, Ny=%d, Nw=%d\n",card,Nx,Ny,Nw);
 cudaSetDevice(card);
 printf("==================\n");
 cuda_printtexture<<<1, 1>>>(Nx,Ny,Nw);
 printf("==================\n");
}
#endif

/* function to write all buffers to output */
/* tid: 0,1,2,3 , select different ordering */
static void
write_buffers_to_output(int bfilled,int Nx,int Ny,pthread_mutex_t *writelock_img,pthread_mutex_t *writelock_psf,float *uvgrid,float *psfgrid,unsigned long int *hostoff,cuFloatComplex *hostimg,cuFloatComplex *hostpsf,int tid) {
 int ci;
 unsigned long int NN=Nx*Ny;
 if (tid%2==0) {
     pthread_mutex_lock(writelock_img);
     for (ci=0; ci<bfilled; ci++) { /* write real,imag parts separately */
      if (hostoff[ci]<NN) { /* FIXME: why need to check this? */
       uvgrid[2*hostoff[ci]]+=hostimg[ci].x;
       uvgrid[2*hostoff[ci]+1]+=hostimg[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_img);
     pthread_mutex_lock(writelock_psf);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       psfgrid[2*hostoff[ci]]+=hostpsf[ci].x;
       psfgrid[2*hostoff[ci]+1]+=hostpsf[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_psf);
 } else if (tid%2==1) {
     pthread_mutex_lock(writelock_psf);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       psfgrid[2*hostoff[ci]]+=hostpsf[ci].x;
       psfgrid[2*hostoff[ci]+1]+=hostpsf[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_psf);
     pthread_mutex_lock(writelock_img);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       uvgrid[2*hostoff[ci]]+=hostimg[ci].x;
       uvgrid[2*hostoff[ci]+1]+=hostimg[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_img);
 } 
}

/* function to write all buffers to output */
/* tid: 0,1,2,3 , select different ordering */
/* not PSF is calculated only for I, so, length/4 for PSF */
static void
write_buffers_to_output_iquv(int bfilled,int Nx,int Ny,pthread_mutex_t *writelock_img,pthread_mutex_t *writelock_psf,float *uvgrid,float *psfgrid,unsigned long int *hostoff,cuFloatComplex *hostimg,cuFloatComplex *hostpsf, 
cuFloatComplex *hostimgQ, cuFloatComplex *hostimgU,cuFloatComplex *hostimgV,
int tid) {
 int ci;
 unsigned long int NN=Nx*Ny;
 if (tid%2==0) {
     pthread_mutex_lock(writelock_img);
     for (ci=0; ci<bfilled; ci++) { /* write real,imag parts separately */
      if (hostoff[ci]<NN) { /* FIXME: why need to check this? */
       uvgrid[2*hostoff[ci]]+=hostimg[ci].x;
       uvgrid[2*hostoff[ci]+1]+=hostimg[ci].y;
       uvgrid[2*(hostoff[ci]+NN)]+=hostimgQ[ci].x;
       uvgrid[2*(hostoff[ci]+NN)+1]+=hostimgQ[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_img);
     pthread_mutex_lock(writelock_psf);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       psfgrid[2*hostoff[ci]]+=hostpsf[ci].x;
       psfgrid[2*hostoff[ci]+1]+=hostpsf[ci].y;
       /* write U,V also here to divide work evenly */
       uvgrid[2*(hostoff[ci]+2*NN)]+=hostimgU[ci].x;
       uvgrid[2*(hostoff[ci]+2*NN)+1]+=hostimgU[ci].y;
       uvgrid[2*(hostoff[ci]+3*NN)]+=hostimgV[ci].x;
       uvgrid[2*(hostoff[ci]+3*NN)+1]+=hostimgV[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_psf);
 } else if (tid%2==1) {
     pthread_mutex_lock(writelock_psf);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       psfgrid[2*hostoff[ci]]+=hostpsf[ci].x;
       psfgrid[2*hostoff[ci]+1]+=hostpsf[ci].y;
       /* write U,V also here to divide work evenly */
       uvgrid[2*(hostoff[ci]+2*NN)]+=hostimgU[ci].x;
       uvgrid[2*(hostoff[ci]+2*NN)+1]+=hostimgU[ci].y;
       uvgrid[2*(hostoff[ci]+3*NN)]+=hostimgV[ci].x;
       uvgrid[2*(hostoff[ci]+3*NN)+1]+=hostimgV[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_psf);
     pthread_mutex_lock(writelock_img);
     for (ci=0; ci<bfilled; ci++) {
      if (hostoff[ci]<NN) {
       uvgrid[2*hostoff[ci]]+=hostimg[ci].x;
       uvgrid[2*hostoff[ci]+1]+=hostimg[ci].y;
       uvgrid[2*(hostoff[ci]+NN)]+=hostimgQ[ci].x;
       uvgrid[2*(hostoff[ci]+NN)+1]+=hostimgQ[ci].y;
      }
     }
     pthread_mutex_unlock(writelock_img);
 } 
}



/* function to set up a GPU, should be called only once */
static void
attach_gpu_to_thread(int card, float *wkernel, int Np, int Nw, float *wparr, float **dwparr, void** carray) {
 cudaError_t err;

 cudaExtent volumesize;
 cudaExtent volumesizeBytes;
 cudaChannelFormatDesc channel;
 cudaMemcpy3DParms copyparms={0};
 cudaArray *cudaarr=0;

 cudaSetDevice(card);

 // copy w axis to device
 err=cudaMalloc((void**)dwparr, sizeof(float)*(Nw+2));
 checkCudaError(err,__FILE__,__LINE__);
 err=cudaMemcpy(*dwparr, (void*)wparr, sizeof(float)*(Nw+2), cudaMemcpyHostToDevice);
 checkCudaError(err,__FILE__,__LINE__);

 //set cuda array volume size NOTE: first dimension is in elements, not in
 // bytes as we use cudaMalloc3DArray and not cudaMalloc
 volumesize=make_cudaExtent(Np,Np,Nw);
 volumesizeBytes=make_cudaExtent(sizeof(cuFloatComplex)*Np,Np,Nw);

 cudaPitchedPtr d_volumeMem;
 err=cudaMalloc3D(&d_volumeMem, volumesizeBytes);
 checkCudaError(err,__FILE__,__LINE__);

 err=cudaMemcpy(d_volumeMem.ptr, (void*)wkernel, sizeof(cuFloatComplex)*Np*Np*Nw, cudaMemcpyHostToDevice);
 checkCudaError(err,__FILE__,__LINE__);

 //create channel to describe data type
 channel=cudaCreateChannelDesc<cuFloatComplex>();

 //allocate device memory for cuda array
 //err=cudaMalloc3DArray(&cudaarr,&channel,volumesize,cudaArrayLayered);
 err=cudaMalloc3DArray(&cudaarr,&channel,volumesize);
 checkCudaError(err,__FILE__,__LINE__);

 //set cuda array copy parameters
 copyparms.extent=volumesize;
 copyparms.dstArray=cudaarr;
 copyparms.kind=cudaMemcpyDeviceToDevice;
 copyparms.srcPtr=d_volumeMem;
 /* copy data */
 err=cudaMemcpy3D(&copyparms);
 checkCudaError(err,__FILE__,__LINE__);
 cudaFree(d_volumeMem.ptr);


 //all coordinate axes are in [0,1]
 texreference.normalized=true;
 //set texture filter mode property
 //use cudaFilterModePoint or cudaFilterModeLinear
 texreference.filterMode=cudaFilterModeLinear;
 //set texture address mode property
 //use cudaAddressModeClamp or cudaAddressModeWrap
 texreference.addressMode[0]=cudaAddressModeClamp;
 texreference.addressMode[1]=cudaAddressModeClamp;
 texreference.addressMode[2]=cudaAddressModeClamp;

 //bind texture reference with cuda array
 err=cudaBindTextureToArray(texreference,cudaarr,channel);
 checkCudaError(err,__FILE__,__LINE__);

 *carray=cudaarr;
}

static void
detach_gpu_from_thread(int card, float *wval, void* carray) {
 cudaSetDevice(card);
 cudaFree(wval);
 cudaArray *cudaarr=(cudaArray*)carray;
 //unbind texture reference to free resource
 cudaUnbindTexture(texreference);
 cudaFreeArray(cudaarr);
}


/* slave thread 2GPU function */
static void *
pipeline_slave_code(void *data)
{
 slave_tdata *td=(slave_tdata*)data;
 gbgdata *gd=(gbgdata*)(td->pline->data);
 int tid=td->tid;
 /* slave barrier */
 th_slave_pipeline tp;
 gb_slave_gdata *tpg;
 int ci;
 int Nt=gd->vis[tid].N;
 if ((tpg=(gb_slave_gdata*)calloc((size_t)Nt,sizeof(gb_slave_gdata)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 float *wval; /* pointer to wparr w values in GPU */
 while(1) {
  sync_barrier(&(td->pline->gate1)); /* stop at gate 1*/
  if(td->pline->terminate) break; /* if flag is set, break loop */
  sync_barrier(&(td->pline->gate2)); /* stop at gate 2 */
  if (gd->status[tid]==PT_DO_WORK_GRID) {
/************************* work *********************/
  /* update data for slave threads */
  /* wait for slave threads  to finish work */
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_WORK_GRID;
    tpg[ci].vis.flag=gd->vis[tid].flag[ci];
    tpg[ci].vis.u=gd->vis[tid].u[ci];
    tpg[ci].vis.v=gd->vis[tid].v[ci];
    tpg[ci].vis.w=gd->vis[tid].w[ci];
    tpg[ci].vis.wt=gd->vis[tid].wt[ci];
    tpg[ci].vis.xx.x=gd->vis[tid].xx[ci].x;
    tpg[ci].vis.xx.y=gd->vis[tid].xx[ci].y;
    tpg[ci].vis.yy.x=gd->vis[tid].yy[ci].x;
    tpg[ci].vis.yy.y=gd->vis[tid].yy[ci].y;
    tpg[ci].vis.xy.x=gd->vis[tid].xy[ci].x;
    tpg[ci].vis.xy.y=gd->vis[tid].xy[ci].y;
    tpg[ci].vis.yx.x=gd->vis[tid].yx[ci].x;
    tpg[ci].vis.yx.y=gd->vis[tid].yx[ci].y;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_NOTHING;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
/************************* work *********************/
  } else if (gd->status[tid]==PT_DO_AGPU) {
  printf("thread %d : w %dx%d\n",tid,gd->Np,gd->Nw);
#ifndef ONE_GPU
   attach_gpu_to_thread(tid,gd->wkernel,gd->Np,gd->Nw,gd->wparr,&wval,&gd->carray[tid]);
#endif
#ifdef ONE_GPU
   attach_gpu_to_thread(0,gd->wkernel,gd->Np,gd->Nw,gd->wparr,&wval,&gd->carray[tid]);
#endif
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_NOTHING;
#ifdef ONE_GPU
    tpg[ci].card=0;
#endif
#ifndef ONE_GPU
    tpg[ci].card=tid;
#endif
    /* also copy fixed data for gridding */
    tpg[ci].lambda=gd->lambda;
    tpg[ci].uvscale=gd->uvscale;
    tpg[ci].maxW=gd->wparr[gd->Nw]; /* remember max W value here */
    tpg[ci].expW=gd->expW; 
    tpg[ci].deltaU=gd->deltaU;
    tpg[ci].uvgrid=gd->uvgrid;
    tpg[ci].psfgrid=gd->psfgrid;
    tpg[ci].writelock_img=gd->writelock_img;
    tpg[ci].writelock_psf=gd->writelock_psf;
    tpg[ci].Nx=gd->Nx;
    tpg[ci].Ny=gd->Ny;
    tpg[ci].Nw=gd->Nw;
    tpg[ci].wparr=wval;
    tpg[ci].wpsupportX=gd->wpsupportX;
    tpg[ci].wpsupportY=gd->wpsupportY;
    tpg[ci].maxsupport=gd->maxsupport;
    tpg[ci].Np=gd->Np;
    tpg[ci].Nz=gd->Nz;
    tpg[ci].imgmode=gd->imgmode;
   }
   /* debugging: print texture Nx,Ny,Nw values */
   //if (tid==0) {debug_show_texture(0,20,20,1);}
   /* spawn slave threads */
   init_slave_pipeline(&tp,Nt,tpg);
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_AGPU;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_NOTHING;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
  } else if (gd->status[tid]==PT_DO_DGPU) {
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_DGPU;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
   sync_barrier(&(tp.gate1)); /* sync at gate 1*/
   for (ci=0; ci<Nt; ci++) {  
    tpg[ci].status=PT_DO_NOTHING;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2*/
#ifndef ONE_GPU
   detach_gpu_from_thread(tid,wval,gd->carray[tid]);
#endif
#ifdef ONE_GPU
   detach_gpu_from_thread(0,wval,gd->carray[tid]);
#endif
   /* destroy slave threads */
   destroy_slave_pipeline(&tp,Nt);
  }
 }
 free(tpg);
 return NULL;
}



/* initialize the pipeline
  and start the slaves rolling 
  Ngpu: how many slaves */
void
init_pipeline(th_pipeline *pline, int Ngpu,
     void *data)
{
 pthread_attr_init(&(pline->attr));
 pthread_attr_setdetachstate(&(pline->attr),PTHREAD_CREATE_JOINABLE);

 init_th_barrier(&(pline->gate1),Ngpu+1); /* 3 threads, including master */
 init_th_barrier(&(pline->gate2),Ngpu+1); /* 3 threads, including master */
 pline->terminate=0;
 pline->data=data; /* data should have pointers to t1 and t2 */
 pline->N=Ngpu;

 if ((pline->sd=(slave_tdata**)calloc((size_t)Ngpu,sizeof(slave_tdata*)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
 }
 int ci;
 for (ci=0; ci<Ngpu; ci++) {
  slave_tdata *t0;
  if ((t0=(slave_tdata*)malloc(sizeof(slave_tdata)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  t0->pline=pline;
  t0->tid=ci;
  pline->sd[ci]=t0;
 }
 if ((pline->slave=(pthread_t*)calloc((size_t)Ngpu,sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
    exit(1);
 }
 for (ci=0; ci<Ngpu; ci++) {
  pthread_create(&(pline->slave[ci]),&(pline->attr),pipeline_slave_code,(void*)pline->sd[ci]);
 }
}

/* destroy the pipeline */
/* need to kill the slaves first */
void
destroy_pipeline(th_pipeline *pline)
{

 pline->terminate=1;
 sync_barrier(&(pline->gate1));
 int ci;
 for (ci=0; ci<pline->N; ci++) {
  pthread_join(pline->slave[ci],NULL);
 }
 destroy_th_barrier(&(pline->gate1));
 destroy_th_barrier(&(pline->gate2));
 pthread_attr_destroy(&(pline->attr));
 for (ci=0; ci<pline->N; ci++) {
  free(pline->sd[ci]);
 }
 free(pline->sd);
 free(pline->slave);
 pline->data=NULL;
}


/* slave thread >2 GPU function */
static void *
pipeline_slave_slave_code(void *data)
{
 slave_slave_tdata *td=(slave_slave_tdata*)data;
 gb_slave_gdata *gd0=(gb_slave_gdata*)(td->pline->data);
 int tid=td->tid;
 gb_slave_gdata *gd=&gd0[tid];
 cuFloatComplex *devimg,*devpsf;
 cuFloatComplex *hostimg,*hostpsf;
 unsigned long int *devoff,*hostoff;

 /* for QUV imaging */
 cuFloatComplex *devimgQ,*devimgU,*devimgV;
 cuFloatComplex *hostimgQ,*hostimgU,*hostimgV;
 /* offsets at NxNy, 2NxNy, 3NxNy */

 cudaError_t err;
 int *doutpos;
 /* determine max support to allocate buffer
    not this is full support for +ve and -ve halves */
 int Maxs;
 int nThreads,nBloks;
 int *outpos;
 /* buffer length for host buffers > device buffers */
 int BL=0;
 int bfilled=0;

 while(1) {
  sync_barrier(&(td->pline->gate1)); /* stop at gate 1*/
  if(td->pline->terminate) break; /* if flag is set, break loop */
  sync_barrier(&(td->pline->gate2)); /* stop at gate 2 */
  if (gd->status==PT_DO_WORK_GRID && !gd->vis.flag) {
/************************* work *********************/
   cudaSetDevice(gd->card);
   /* find support for this w */
   //printf("nBlocks=%d nThreads=%d\n",nBloks,nThreads);
   cuda_search<<<nBloks, nThreads>>>(doutpos, gd->wparr, gd->vis.w, gd->lambda, gd->Nw+2);
   cudaDeviceSynchronize();
   err=cudaGetLastError();
   //checkCudaError(err,__FILE__,__LINE__);
   err=cudaMemcpy(outpos,doutpos,sizeof(int),cudaMemcpyDeviceToHost);
   checkCudaError(err,__FILE__,__LINE__);
   if (*outpos>=gd->Nw) { *outpos=gd->Nw-1; }
   int Mx=gd->wpsupportX[*outpos];
   //printf("card %d kernel %d buff %d Nx=%d Ny=%d Np=%d Nw=%d, plane %d support %d\n",gd->card,tid,Maxs,gd->Nx,gd->Ny,gd->Np,gd->Nw,*outpos-1,Mx);

   if (Mx>0) { /* only when support is finite */
   /* depending on the actual support, adjust dimensions 
      to cover 2D array of Mx x Mx */
    dim3 threadsPerBlock(8, 8);
    dim3 numBlocks((Mx+threadsPerBlock.x-1)/threadsPerBlock.x,
               (Mx+threadsPerBlock.y-1)/threadsPerBlock.y);
//printf("Mx=%d threads=%d,%d blocks=%d,%d\n",Mx,threadsPerBlock.x,threadsPerBlock.y,numBlocks.x,numBlocks.y);
    /* reset offset values to Nx x Ny  to eliminate points off the grid */
    cuda_resetmemory<<<numBlocks,threadsPerBlock>>>(devoff,Mx,gd->Nx,gd->Ny);
    cudaDeviceSynchronize();
    int MxM=2*Mx*Mx;
    /* max W is at gd->wparr[gd->Nw] */
    if (gd->imgmode>1) {
     cuda_griddata_iquv<<<numBlocks,threadsPerBlock>>>(devimg,devpsf,devoff,devimgQ,devimgU,devimgV,Mx,gd->vis.u,gd->vis.v,gd->vis.w,gd->vis.wt,make_cuFloatComplex(gd->vis.xx.x,gd->vis.xx.y),make_cuFloatComplex(gd->vis.yy.x,gd->vis.yy.y), make_cuFloatComplex(gd->vis.xy.x,gd->vis.xy.y),make_cuFloatComplex(gd->vis.yx.x,gd->vis.yx.y),gd->maxW,gd->expW,gd->lambda,gd->uvscale,gd->deltaU,gd->Nx,gd->Ny,gd->Nw,gd->Np,gd->Nz);
    } else {
     cuda_griddata<<<numBlocks,threadsPerBlock>>>(devimg,devpsf,devoff,Mx,gd->vis.u,gd->vis.v,gd->vis.w,gd->vis.wt,make_cuFloatComplex(gd->vis.xx.x,gd->vis.xx.y),make_cuFloatComplex(gd->vis.yy.x,gd->vis.yy.y), make_cuFloatComplex(gd->vis.xy.x,gd->vis.xy.y),make_cuFloatComplex(gd->vis.yx.x,gd->vis.yx.y), gd->maxW,gd->expW,gd->lambda,gd->uvscale,gd->deltaU,gd->Nx,gd->Ny,gd->Nw,gd->Np,gd->Nz);
    }
    cudaDeviceSynchronize();
    err=cudaGetLastError();
    checkCudaError(err,__FILE__,__LINE__);
    /* check if there is enough memory in buffers for copying to host */
    if (bfilled+MxM>=BL) { 
     /* write to output */
     if (gd->imgmode>1) {
      write_buffers_to_output_iquv(bfilled,gd->Nx,gd->Ny,gd->writelock_img,gd->writelock_psf,gd->uvgrid,gd->psfgrid,hostoff,hostimg,hostpsf,hostimgQ,hostimgU,hostimgV,rand()%MAX_GPU);
     } else {
      write_buffers_to_output(bfilled,gd->Nx,gd->Ny,gd->writelock_img,gd->writelock_psf,gd->uvgrid,gd->psfgrid,hostoff,hostimg,hostpsf,rand()%MAX_GPU);
     }
     bfilled=0;
    } 
     err=cudaMemcpy(&hostimg[bfilled],devimg,sizeof(cuFloatComplex)*MxM,cudaMemcpyDeviceToHost);
     checkCudaError(err,__FILE__,__LINE__);
     err=cudaMemcpy(&hostpsf[bfilled],devpsf,sizeof(cuFloatComplex)*MxM,cudaMemcpyDeviceToHost);
     checkCudaError(err,__FILE__,__LINE__);
     err=cudaMemcpy(&hostoff[bfilled],devoff,sizeof(unsigned long int)*MxM,cudaMemcpyDeviceToHost);
     checkCudaError(err,__FILE__,__LINE__);
     if (gd->imgmode>1) {
      err=cudaMemcpy(&hostimgQ[bfilled],devimgQ,sizeof(cuFloatComplex)*MxM,cudaMemcpyDeviceToHost);
      checkCudaError(err,__FILE__,__LINE__);
      err=cudaMemcpy(&hostimgU[bfilled],devimgU,sizeof(cuFloatComplex)*MxM,cudaMemcpyDeviceToHost);
      checkCudaError(err,__FILE__,__LINE__);
      err=cudaMemcpy(&hostimgV[bfilled],devimgV,sizeof(cuFloatComplex)*MxM,cudaMemcpyDeviceToHost);
      checkCudaError(err,__FILE__,__LINE__);
     }
     cudaDeviceSynchronize();
     /* advance buffer */
     bfilled+=MxM; 

   }
/************************* work *********************/
  } else if (gd->status==PT_DO_AGPU ) {
   cudaSetDevice(gd->card);
   Maxs=gd->maxsupport;
   int MxMs=2*Maxs*Maxs;
   nThreads=16;
   nBloks=(gd->Nw+2+nThreads-1)/nThreads;

   /* for Maxs ~ 512, Maxs*Maxs*2*sizeof(cuFloatComplex) ~= 2MB
      so make buffer size ~ 1 MB by dividing by 2 */
   if (MxMs<DATA_BUF_LEN) {
    BL=DATA_BUF_LEN;
   } else {
    BL=MxMs;
   }
   err=cudaMalloc((void**)&devimg,sizeof(cuFloatComplex)*MxMs);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaHostAlloc((void**)&hostimg,sizeof(cuFloatComplex)*BL,cudaHostAllocDefault);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaMalloc((void**)&devpsf,sizeof(cuFloatComplex)*MxMs);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaHostAlloc((void**)&hostpsf,sizeof(cuFloatComplex)*BL,cudaHostAllocDefault);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaMalloc((void**)&devoff,sizeof(unsigned long int)*MxMs);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaHostAlloc((void**)&hostoff,sizeof(unsigned long int)*BL,cudaHostAllocDefault);
   checkCudaError(err,__FILE__,__LINE__);
 
   err=cudaMalloc((void**)&doutpos,sizeof(int));
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaHostAlloc((void**)&outpos,sizeof(int),cudaHostAllocDefault);
   checkCudaError(err,__FILE__,__LINE__);

   if (gd->imgmode>1) {
    /* IQUV imaging */
    err=cudaMalloc((void**)&devimgQ,sizeof(cuFloatComplex)*MxMs);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaHostAlloc((void**)&hostimgQ,sizeof(cuFloatComplex)*BL,cudaHostAllocDefault);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaMalloc((void**)&devimgU,sizeof(cuFloatComplex)*MxMs);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaHostAlloc((void**)&hostimgU,sizeof(cuFloatComplex)*BL,cudaHostAllocDefault);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaMalloc((void**)&devimgV,sizeof(cuFloatComplex)*MxMs);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaHostAlloc((void**)&hostimgV,sizeof(cuFloatComplex)*BL,cudaHostAllocDefault);
    checkCudaError(err,__FILE__,__LINE__);
   }

   cudaDeviceSynchronize();
  } else if (gd->status==PT_DO_DGPU ) {
   /* write last buffer to output */
   if (gd->imgmode>1) {
    write_buffers_to_output_iquv(bfilled,gd->Nx,gd->Ny,gd->writelock_img,gd->writelock_psf,gd->uvgrid,gd->psfgrid,hostoff,hostimg,hostpsf,hostimgQ,hostimgU,hostimgV,rand()%MAX_GPU);
   } else {
    write_buffers_to_output(bfilled,gd->Nx,gd->Ny,gd->writelock_img,gd->writelock_psf,gd->uvgrid,gd->psfgrid,hostoff,hostimg,hostpsf,rand()%MAX_GPU);
   }
   bfilled=0;
   err=cudaFree(doutpos);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFree(devimg);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFree(devpsf);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFree(devoff);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFreeHost(hostimg);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFreeHost(hostpsf);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFreeHost(hostoff);
   checkCudaError(err,__FILE__,__LINE__);
   err=cudaFreeHost(outpos);
   checkCudaError(err,__FILE__,__LINE__);
   if (gd->imgmode>1) {
    err=cudaFree(devimgQ);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaFreeHost(hostimgQ);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaFree(devimgU);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaFreeHost(hostimgU);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaFree(devimgV);
    checkCudaError(err,__FILE__,__LINE__);
    err=cudaFreeHost(hostimgV);
    checkCudaError(err,__FILE__,__LINE__);
   }
   cudaDeviceSynchronize();
  }
 }
 return NULL;
}



/* initialize the pipeline
  and start the slaves rolling 
N: total slaves
*/
void
init_slave_pipeline(th_slave_pipeline *pline, int N,
     void *data)
{
 if ((pline->sd=(slave_slave_tdata**)malloc(sizeof(slave_slave_tdata*)*N))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }
 if ((pline->slave=(pthread_t*)malloc(sizeof(pthread_t)*N))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }

 int ci;
 pthread_attr_init(&(pline->attr));
 pthread_attr_setdetachstate(&(pline->attr),PTHREAD_CREATE_JOINABLE);

 init_th_barrier(&(pline->gate1),N+1); /* N+1 threads, including master */
 init_th_barrier(&(pline->gate2),N+1); /* N+1 threads, including master */
 pline->terminate=0;
 pline->data=data; /* data is an array of data for t1,t2,t3,... */

 for(ci=0; ci<N; ci++) {
  if ((pline->sd[ci]=(slave_slave_tdata*)malloc(sizeof(slave_slave_tdata)))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
  }
  pline->sd[ci]->pline=pline;
  pline->sd[ci]->tid=ci;
  pthread_create(&(pline->slave[ci]),&(pline->attr),pipeline_slave_slave_code,(void*)pline->sd[ci]);
 }
}

/* destroy the pipeline */
/* need to kill the slaves first */
void
destroy_slave_pipeline(th_slave_pipeline *pline,int N)
{

 pline->terminate=1;
 sync_barrier(&(pline->gate1));
 int ci;
 for(ci=0; ci<N; ci++) {
  pthread_join(pline->slave[ci],NULL);
  free(pline->sd[ci]);
 }
 destroy_th_barrier(&(pline->gate1));
 destroy_th_barrier(&(pline->gate2));
 pthread_attr_destroy(&(pline->attr));
 pline->data=NULL;
 free(pline->slave);
 free(pline->sd);
}

}
