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

#include "helper.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include "largefft.h"
#ifdef HAVE_CUDA
#include "cugridder.h"
#endif





/* note: all complex float arrays are passed as float 
  uvgrid
  psfgrid
  wkernel  
*/
/* darr: Nrows x 1  of data from MS 
   uvscale: scale uv coords to fit [-0.5,0.5]
   lambda: wavelength
   deltaU: uv pixel size
   uvgrid: array to store gridded data : Nx*Ny complex float
   psfgrid: array to store PSF (gridded 1) : Nx*Ny complex float
   Nx,Ny : grid size (also image size)
   expW: exponent in w spacing w^(expW) : 0.5 for sqrt() 1.0 for linear
   wparr: W plane values (with guard values) : Nw+2
   Nw: how many W planes
   wpsupport (X,Y): conv. kernel support for each W plane : Nw 
   maxsupport: max support in pixels
   wkernel: conv. kernels : each kernel NpxNp, Nw planes, complex float
   Np: conv. kernel width in uv pixels
   Nz: zero padded conv kernel length
   Nt: no of GPU threads  (per GPU)
   Ngpu: no of GPUs to use
   
   Note about weight calculation
   for each data point, and pixel i in the grid that has a contribution
   uvgrid +=conv_kernel*data
   psfgrid +=(conv_kernel)
*/

int
cuda_griddata(iodata *darr, int Nrows, float uvscale, float lambda, float deltaU, complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, float expW, float *wparr, int Nw, int *wpsupportX, int *wpsupportY, int maxsupport, complex float *wkernel, int Np, int Nz, int Nt, int Ngpu, int imgmode) {
 /********* thread data ******************/
  /* barrier */
  th_pipeline tp;
  gbgdata tpg;
  int N=Nt; /* threads per GPU so total threads is Ngpu x Nt */
  int NrowsPerThread=Nrows/(Ngpu*N); /* using both cards */
  /* if there are too few rows per thread, reduce no of threads */
  if (NrowsPerThread<2) {
    N=Nrows/16; /* FIXME: recheck this */
    NrowsPerThread=Nrows/(Ngpu*N);
  }
  int ci,i,j,cj;
/****************************************/
  pthread_mutex_t writelock_img,writelock_psf;
  pthread_mutex_init(&writelock_img, NULL);
  pthread_mutex_init(&writelock_psf, NULL);

/********** setup threads *******************************/
  if ((tpg.status=(int*)calloc((size_t)Ngpu,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((tpg.vis=(visdata*)calloc((size_t)Ngpu,sizeof(visdata)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  for (ci=0; ci<Ngpu; ci++) {
   tpg.vis[ci].N=N;
   tpg.status[ci]=PT_DO_NOTHING;
  }
  for (ci=0; ci<Ngpu; ci++) {
   if ((tpg.vis[ci].flag=(int*)calloc((size_t)N,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].u=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].v=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].w=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].wt=(float*)calloc((size_t)N,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].xx=(float2*)calloc((size_t)N,sizeof(float2)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].yy=(float2*)calloc((size_t)N,sizeof(float2)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].xy=(float2*)calloc((size_t)N,sizeof(float2)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((tpg.vis[ci].yx=(float2*)calloc((size_t)N,sizeof(float2)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
  }
  tpg.lambda=lambda;
  tpg.uvscale=uvscale;
  tpg.deltaU=deltaU;
  tpg.expW=expW;
  tpg.uvgrid=(float*)uvgrid;
  tpg.psfgrid=(float*)psfgrid;
  tpg.writelock_img=&writelock_img;
  tpg.writelock_psf=&writelock_psf;
  tpg.Nx=Nx;
  tpg.Ny=Ny;
  tpg.wparr=wparr;
  tpg.Nw=Nw;
  tpg.wpsupportX=wpsupportX;
  tpg.wpsupportY=wpsupportY;
  tpg.maxsupport=maxsupport;
  tpg.wkernel=(float*)wkernel;
  tpg.Np=Np;
  tpg.Nz=Nz;
  tpg.imgmode=imgmode;

  init_pipeline(&tp,Ngpu,&tpg);
  sync_barrier(&(tp.gate1)); /* sync at gate 1*/
  for (ci=0; ci<Ngpu; ci++) {
   tpg.status[ci]=PT_DO_AGPU;
  }
  sync_barrier(&(tp.gate2)); /* sync at gate 2*/

  sync_barrier(&(tp.gate1)); /* sync at gate 1*/
  for (ci=0; ci<Ngpu; ci++) {
   tpg.status[ci]=PT_DO_NOTHING;
  }
  sync_barrier(&(tp.gate2)); /* sync at gate 2*/
/********** done setup threads *******************************/
  ci=0;
  /* for printing progress */
  int lBar=0;
  int deltaBar=NrowsPerThread/10;
  while (ci<NrowsPerThread) {
   sync_barrier(&(tp.gate1)); /* sync at gate 1 */
   for (cj=0; cj<Ngpu; cj++) {
    tpg.status[cj]=PT_DO_WORK_GRID;
   }
   i=ci*Ngpu*N;
   /* update gridding related data in pipeline */
   for (j=0; j<N; j++) {
    for (cj=0; cj<Ngpu; cj++) {
    tpg.vis[cj].flag[j]=darr[i].flag;
    tpg.vis[cj].u[j]=darr[i].u;
    tpg.vis[cj].v[j]=darr[i].v;
    tpg.vis[cj].w[j]=darr[i].w;
    tpg.vis[cj].wt[j]=darr[i].wt;
    tpg.vis[cj].xx[j].x=darr[i].xx.x;
    tpg.vis[cj].xx[j].y=darr[i].xx.y;
    tpg.vis[cj].yy[j].x=darr[i].yy.x;
    tpg.vis[cj].yy[j].y=darr[i].yy.y;
    tpg.vis[cj].xy[j].x=darr[i].xy.x;
    tpg.vis[cj].xy[j].y=darr[i].xy.y;
    tpg.vis[cj].yx[j].x=darr[i].yx.x;
    tpg.vis[cj].yx[j].y=darr[i].yx.y;
    i++;
    }
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2 */

   ci++;
   /* redirect output to stderr, to remove junk from stdout */
   /* also only print every 1/10 of the values to reduce clutter */
   if (ci>lBar) {
    //fprintf(stderr,"\r%d/%d",ci,NrowsPerThread);
    printf(". ");
    lBar+=deltaBar;
    fflush(stdout);
   }

   sync_barrier(&(tp.gate1)); /* sync at gate 1 */
   for (cj=0; cj<Ngpu; cj++) {
    tpg.status[cj]=PT_DO_NOTHING;
   }
   sync_barrier(&(tp.gate2)); /* sync at gate 2 */
  
 }
 printf("\n");
  /******** free threads ***************/
  sync_barrier(&(tp.gate1)); /* sync at gate 1*/
  for (cj=0; cj<Ngpu; cj++) {
   tpg.status[cj]=PT_DO_DGPU;
  }
  sync_barrier(&(tp.gate2)); /* sync at gate 2*/


  destroy_pipeline(&tp);
  /******** done free threads ***************/

  for (ci=0; ci<Ngpu; ci++) {
    free(tpg.vis[ci].flag);
    free(tpg.vis[ci].u);
    free(tpg.vis[ci].v);
    free(tpg.vis[ci].w);
    free(tpg.vis[ci].wt);
    free(tpg.vis[ci].xx);
    free(tpg.vis[ci].yy);
    free(tpg.vis[ci].yx);
    free(tpg.vis[ci].xy);
  }
  pthread_mutex_destroy(&writelock_img);
  pthread_mutex_destroy(&writelock_psf);

  free(tpg.status);
  free(tpg.vis);

 return 0;
}
