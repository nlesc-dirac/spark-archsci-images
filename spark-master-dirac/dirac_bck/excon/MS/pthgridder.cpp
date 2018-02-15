
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

#include "pthgridder.h"
#include <sys/mman.h>

using namespace std;

typedef struct thread_data_grid_t_ {
  iodata *darr;
  int startrow,Nrows;
  float lambda;
  float uvscale;
  float deltaU;
  complex float *uvgrid;
  float *wgrid;
  complex float *psfgrid;
  int Nx,Ny;
  float *wparr;
  int Nw;
  int *wpsupportX;
  int *wpsupportY;
  complex float *wkernel;
  float *uvxgrid;
  int Np;
  pthread_mutex_t *writelock_img;
  pthread_mutex_t *writelock_psf;
} thread_data_grid_t;


/* write all buffers to output */
/* tid: 0,1,2,3 : change ordering depending on this */
static void
write_buffers_to_output_cpu(int bufflen,pthread_mutex_t *writelock_img,pthread_mutex_t *writelock_psf,complex float *uvgrid,complex float *psfgrid,complex float *buffimg,complex float *buffpsf,unsigned long int *buffoff,int tid) {
   if (tid%2==0) { /* even */
                  pthread_mutex_lock(writelock_img);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   uvgrid[buffoff[bidx]]+=buffimg[bidx];
                  }
                  pthread_mutex_unlock(writelock_img);
                  pthread_mutex_lock(writelock_psf);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   psfgrid[buffoff[bidx]]+=buffpsf[bidx];
                  }
                  pthread_mutex_unlock(writelock_psf);
   } else if (tid%2==1) { /* odd */
                  pthread_mutex_lock(writelock_psf);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   psfgrid[buffoff[bidx]]+=buffpsf[bidx];
                  }
                  pthread_mutex_unlock(writelock_psf);
                  pthread_mutex_lock(writelock_img);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   uvgrid[buffoff[bidx]]+=buffimg[bidx];
                  }
                  pthread_mutex_unlock(writelock_img);
   } 
}


static void
write_buffers_to_output_cpu_iquv(int bufflen,int Nx, int Ny, pthread_mutex_t *writelock_img,pthread_mutex_t *writelock_psf,complex float *uvgrid,complex float *psfgrid,complex float *buffimg,complex float *buffpsf,unsigned long int *buffoff,complex float *buffimgQ, complex float *buffimgU, complex float *buffimgV, int tid) {
   unsigned long int NN=Nx*Ny;
   if (tid%2==0) {
                  pthread_mutex_lock(writelock_img);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   uvgrid[buffoff[bidx]]+=buffimg[bidx];
                   uvgrid[buffoff[bidx]+NN]+=buffimgQ[bidx];
                   uvgrid[buffoff[bidx]+2*NN]+=buffimgU[bidx];
                   uvgrid[buffoff[bidx]+3*NN]+=buffimgV[bidx];
                  }
                  pthread_mutex_unlock(writelock_img);
                  pthread_mutex_lock(writelock_psf);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   psfgrid[buffoff[bidx]]+=buffpsf[bidx];
                  }
                  pthread_mutex_unlock(writelock_psf);
   } else if (tid%2==1) {
                  pthread_mutex_lock(writelock_psf);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   psfgrid[buffoff[bidx]]+=buffpsf[bidx];
                  }
                  pthread_mutex_unlock(writelock_psf);
                  pthread_mutex_lock(writelock_img);
                  for (int bidx=0; bidx<bufflen; bidx++) {
                   uvgrid[buffoff[bidx]]+=buffimg[bidx];
                   uvgrid[buffoff[bidx]+NN]+=buffimgQ[bidx];
                   uvgrid[buffoff[bidx]+2*NN]+=buffimgU[bidx];
                   uvgrid[buffoff[bidx]+3*NN]+=buffimgV[bidx];
                  }
                  pthread_mutex_unlock(writelock_img);
   } 
}


/* worker thread function for prediction */
static void *
gridder_threadfn(void *data) {
   thread_data_grid_t *t=(thread_data_grid_t*)data;
   complex float *buffimg;
   complex float *buffpsf;
   unsigned long int *buffoff;
   int BL=DATA_BUF_LEN;
   if ((buffoff=(unsigned long int*)calloc((size_t)BL,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffimg=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffpsf=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }

   int bfilled=0;

/*************************************************/
   for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
      if (!t->darr[i].flag) {
                float tempu = t->darr[i].u;
                float tempv = t->darr[i].v;
                float tempw = t->darr[i].w;
                fcomp tempxx=t->darr[i].xx;
                fcomp tempyy=t->darr[i].yy;
/* [-minu,-minv]->(0,0) [maxu,maxv]-> (Nx,Ny) */
/* move (0,0) in uv plane to [Nx/2,Ny/2] pixels (both odd,even sizes) */
/* full uv plane is gridded. 
     /|\    --> y      <--v  /|\  u
      |  x                    | 
     axis directions as above 

     size of d: Nx*Ny complex float values
     size of img: Nx*Ny complex float values 
*/

         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*t->uvscale;
          float vf=tempv*t->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)t->Nx); /*  width */
          float vi=(vf*(float)t->Ny); /*  width */
        /* also find the w plane */
          float z=(fabs(tempw));
          int k=bin_search(t->wparr,z,0,t->Nw+2);
          k--;
          if (k>t->Nw-2) {
            k=t->Nw-2;
          }
          if (k<0) { k=0; }
//cout<<"for (|w|)="<<z<<" wplane "<<k<<" ["<<t->wparr[k+1]<<","<<t->wparr[k+2]<<"]"<<endl;
          float poffX=-1.0f;
          float poffY=-1.0f;
          int Mx=MAX(t->wpsupportX[k],t->wpsupportX[k+1]);
          int My=MAX(t->wpsupportY[k],t->wpsupportY[k+1]);
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*t->Nx+poffX);
          int y=(int)round(-vi+0.5f*t->Ny+poffY);
          /* now find the valid pixel range [x-M/2,x+M/2] [y-M/2,y+M/2]
             that we need to calculate convolution kernel */
          int xlow=(x-Mx/2>=0?x-Mx/2:0);
          int xhigh=(x+Mx/2<t->Nx?x+Mx/2:t->Nx-1);
          int ylow=(y-My/2>=0?y-My/2:0);
          int yhigh=(y+My/2<t->Ny?y+My/2:t->Ny-1);
//cout<<"for ["<<x<<","<<y<<"] -> ("<<xlow<<","<<ylow<<") ("<<xhigh<<","<<yhigh<<")"<<endl;
          unsigned long int coffset=0;
          int x1,y1;
          complex float ppxy;
          complex float ppxy1;

          complex float sI=0.5f*(tempxx.x+_Complex_I*tempxx.y+tempyy.x+_Complex_I*tempyy.y);
          float invsigma=t->darr[i].wt;
          for (int cy=ylow; cy<=yhigh; cy++) {
               /* move (x,y) to ifftshifted value */
               y1=ifftshift_index(t->Ny,cy);
               for (int cx=xlow; cx<=xhigh; cx++) {
                x1=ifftshift_index(t->Nx,cx);
                /* find right offset (+1 not added because start from 0) */
                coffset=x1*(t->Ny)+y1;
                float del_u=(uf-((float)(cx-t->Nx/2))*t->deltaU);
                float del_v=(-vf-((float)(cy-t->Ny/2))*t->deltaU);
                ppxy=conv_eval(t->Np,t->uvxgrid,t->wkernel,t->wparr,z,t->Nw,k,tempw,del_u,del_v);
//printf("u,v=%f,%f uf,vf=%f,%f up,vp=%d,%d ui,vi=%f,%f x0,y0=%d,%d, x1,y1=%d,%d x,y=%f,%f pswf=%f,%f\n",tempu*t->lambda,tempv*t->lambda,uf,vf,x,y,ui,vi,cx,cy,x1,y1,del_u,del_v,crealf(ppxy),cimagf(ppxy));
                /* multiply pswf with inverse sigma */
                ppxy=ppxy*invsigma;
                ppxy1=ppxy*sI;
                if (bfilled<BL) {
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffpsf[bfilled++]=ppxy;
                } else {
                  /* empty buffer */
                  write_buffers_to_output_cpu(BL,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,rand()%MAX_GPU);
                  bfilled=0;
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffpsf[bfilled++]=ppxy;
                }

                } 
          } /* iteration end over convlution support +ve (u,v) */
 
          /* also add (-uf,-vf) : mirror image */
          uf=-uf; vf=-vf; tempw=-tempw;
          x=-(int)round(-(-ui+0.5f*t->Nx+poffX));
          y=-(int)round(-(vi+0.5f*t->Ny+poffY));
          xlow=(x-Mx/2>=0?x-Mx/2:0);
          xhigh=(x+Mx/2<t->Nx?x+Mx/2:t->Nx-1);
          ylow=(y-My/2>=0?y-My/2:0);
          yhigh=(y+My/2<t->Ny?y+My/2:t->Ny-1);


          for (int cy=ylow; cy<=yhigh; cy++) {
                y1=ifftshift_index(t->Ny,cy);
                /* we do not do anymore checks because of symmetry earlier checks should have tested the validity of this point */
                for (int cx=xlow; cx<=xhigh; cx++) {
                 x1=ifftshift_index(t->Nx,cx);
                /* find right offset */
                coffset=x1*t->Ny+y1;
                float del_u=(uf-((float)(cx-t->Nx/2))*t->deltaU);
                float del_v=(-vf-((float)(cy-t->Ny/2))*t->deltaU);
                ppxy=conv_eval(t->Np,t->uvxgrid,t->wkernel,t->wparr,z,t->Nw,k,tempw,del_u,del_v);
                /* multiply pswf with inverse sigma */
                ppxy=ppxy*invsigma;
                ppxy1=ppxy*conj(sI);
                if (bfilled<BL) {
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffpsf[bfilled++]=ppxy;
                } else {
                  /* empty buffer */
                  write_buffers_to_output_cpu(BL,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,rand()%MAX_GPU);
                  bfilled=0;
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffpsf[bfilled++]=ppxy;
                }
               }
          } /* iteration end over convlution support +ve (-u,-v) */
        } 
      } /* end data rows loop */

 /* write last buffer */
 if (bfilled>0) {
   write_buffers_to_output_cpu(bfilled,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,rand()%MAX_GPU);
 }

 free(buffimg);
 free(buffpsf);
 free(buffoff);
 return NULL;
}

/* worker thread function for prediction */
static void *
gridder_threadfn_iquv(void *data) {
   thread_data_grid_t *t=(thread_data_grid_t*)data;
   complex float *buffimg,*buffimgQ,*buffimgU,*buffimgV;
   complex float *buffpsf;
   unsigned long int *buffoff;
   int BL=DATA_BUF_LEN;
   if ((buffoff=(unsigned long int*)calloc((size_t)BL,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffimg=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffpsf=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffimgQ=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffimgU=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((buffimgV=(complex float*)calloc((size_t)BL,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   int bfilled=0;

/*************************************************/
   for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
      if (!t->darr[i].flag) {
                float tempu = t->darr[i].u;
                float tempv = t->darr[i].v;
                float tempw = t->darr[i].w;
                fcomp tempxx=t->darr[i].xx;
                fcomp tempyy=t->darr[i].yy;
                fcomp tempxy=t->darr[i].xy;
                fcomp tempyx=t->darr[i].yx;
/* [-minu,-minv]->(0,0) [maxu,maxv]-> (Nx,Ny) */
/* move (0,0) in uv plane to [Nx/2,Ny/2] pixels (both odd,even sizes) */
/* full uv plane is gridded. 
     /|\    --> y      <--v  /|\  u
      |  x                    | 
     axis directions as above 

     size of d: Nx*Ny complex float values
     size of img: Nx*Ny complex float values 
*/

         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*t->uvscale;
          float vf=tempv*t->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)t->Nx); /*  width */
          float vi=(vf*(float)t->Ny); /*  width */
        /* also find the w plane */
          float z=(fabs(tempw));
          int k=bin_search(t->wparr,z,0,t->Nw+2);
          k--;
          if (k>t->Nw-2) {
            k=t->Nw-2;
          }
          if (k<0) { k=0; }
//cout<<"for (|w|)="<<z<<" wplane "<<k<<" ["<<t->wparr[k+1]<<","<<t->wparr[k+2]<<"]"<<endl;
          float poffX=-1.0f;
          float poffY=-1.0f;
          int Mx=MAX(t->wpsupportX[k],t->wpsupportX[k+1]);
          int My=MAX(t->wpsupportY[k],t->wpsupportY[k+1]);
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*t->Nx+poffX);
          int y=(int)round(-vi+0.5f*t->Ny+poffY);
          /* now find the valid pixel range [x-M/2,x+M/2] [y-M/2,y+M/2]
             that we need to calculate convolution kernel */
          int xlow=(x-Mx/2>=0?x-Mx/2:0);
          int xhigh=(x+Mx/2<t->Nx?x+Mx/2:t->Nx-1);
          int ylow=(y-My/2>=0?y-My/2:0);
          int yhigh=(y+My/2<t->Ny?y+My/2:t->Ny-1);
//cout<<"for ["<<x<<","<<y<<"] -> ("<<xlow<<","<<ylow<<") ("<<xhigh<<","<<yhigh<<")"<<endl;
          unsigned long int coffset=0;
          int x1,y1;
          complex float ppxy;
          complex float ppxy1,ppxy2,ppxy3,ppxy4;

          /* Q =(XX-YY)/2 U=(XY+YX)/2 V=imag(YX-XY)/2 */
          complex float sI=0.5f*(tempxx.x+_Complex_I*tempxx.y+tempyy.x+_Complex_I*tempyy.y);
          complex float sQ=0.5f*(tempxx.x+_Complex_I*tempxx.y-tempyy.x-_Complex_I*tempyy.y);
          //complex float sU=0.5f*(tempxy.x+_Complex_I*tempxy.y-tempyx.y+_Complex_I*tempyx.x);
          //complex float sV=-0.5f*(-tempxy.y+_Complex_I*tempxy.x+tempyx.x+_Complex_I*tempyx.y);
          complex float sU=0.5f*(tempxy.x+_Complex_I*tempxy.y+tempyx.x+_Complex_I*tempyx.y);
          complex float sV=0.5f*(tempxy.y-_Complex_I*tempxy.x-tempyx.y+_Complex_I*tempyx.x);
          float invsigma=t->darr[i].wt;
          for (int cy=ylow; cy<=yhigh; cy++) {
               /* move (x,y) to ifftshifted value */
               y1=ifftshift_index(t->Ny,cy);
               for (int cx=xlow; cx<=xhigh; cx++) {
                x1=ifftshift_index(t->Nx,cx);
                /* find right offset (+1 not added because start from 0) */
                coffset=x1*(t->Ny)+y1;
                float del_u=(uf-((float)(cx-t->Nx/2))*t->deltaU);
                float del_v=(-vf-((float)(cy-t->Ny/2))*t->deltaU);
                ppxy=conv_eval(t->Np,t->uvxgrid,t->wkernel,t->wparr,z,t->Nw,k,tempw,del_u,del_v);
//printf("u,v=%f,%f uf,vf=%f,%f up,vp=%d,%d ui,vi=%f,%f x0,y0=%d,%d, x1,y1=%d,%d x,y=%f,%f pswf=%f,%f\n",tempu*t->lambda,tempv*t->lambda,uf,vf,x,y,ui,vi,cx,cy,x1,y1,del_u,del_v,crealf(ppxy),cimagf(ppxy));
                /* multiply pswf with inverse sigma */
                ppxy=ppxy*invsigma;
                ppxy1=ppxy*sI;
                ppxy2=ppxy*sQ;
                ppxy3=ppxy*sU;
                ppxy4=ppxy*sV;
                if (bfilled<BL) {
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffimgQ[bfilled]=ppxy2;
                  buffimgU[bfilled]=ppxy3;
                  buffimgV[bfilled]=ppxy4;
                  buffpsf[bfilled++]=ppxy;
                } else {
                  /* empty buffer */
                  write_buffers_to_output_cpu_iquv(BL,t->Nx,t->Ny,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,buffimgQ,buffimgU,buffimgV,rand()%MAX_GPU);
                  bfilled=0;
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffimgQ[bfilled]=ppxy2;
                  buffimgU[bfilled]=ppxy3;
                  buffimgV[bfilled]=ppxy4;
                  buffpsf[bfilled++]=ppxy;
                }

                } 
          } /* iteration end over convlution support +ve (u,v) */
 
          /* also add (-uf,-vf) : mirror image */
          uf=-uf; vf=-vf; tempw=-tempw;
          x=-(int)round(-(-ui+0.5f*t->Nx+poffX));
          y=-(int)round(-(vi+0.5f*t->Ny+poffY));
          xlow=(x-Mx/2>=0?x-Mx/2:0);
          xhigh=(x+Mx/2<t->Nx?x+Mx/2:t->Nx-1);
          ylow=(y-My/2>=0?y-My/2:0);
          yhigh=(y+My/2<t->Ny?y+My/2:t->Ny-1);


          for (int cy=ylow; cy<=yhigh; cy++) {
                y1=ifftshift_index(t->Ny,cy);
                /* we do not do anymore checks because of symmetry earlier checks should have tested the validity of this point */
                for (int cx=xlow; cx<=xhigh; cx++) {
                 x1=ifftshift_index(t->Nx,cx);
                /* find right offset */
                coffset=x1*t->Ny+y1;
                float del_u=(uf-((float)(cx-t->Nx/2))*t->deltaU);
                float del_v=(-vf-((float)(cy-t->Ny/2))*t->deltaU);
                ppxy=conv_eval(t->Np,t->uvxgrid,t->wkernel,t->wparr,z,t->Nw,k,tempw,del_u,del_v);
                /* multiply pswf with inverse sigma */
                ppxy=ppxy*invsigma;
                ppxy1=ppxy*conj(sI);
                ppxy2=ppxy*conj(sQ);
                ppxy3=ppxy*conj(sU);
                ppxy4=ppxy*conj(sV);
                if (bfilled<BL) {
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffimgQ[bfilled]=ppxy2;
                  buffimgU[bfilled]=ppxy3;
                  buffimgV[bfilled]=ppxy4;
                  buffpsf[bfilled++]=ppxy;
                } else {
                  /* empty buffer */
                  write_buffers_to_output_cpu_iquv(BL,t->Nx,t->Ny,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,buffimgQ,buffimgU,buffimgV,rand()%MAX_GPU);
                  bfilled=0;
                  buffoff[bfilled]=coffset;
                  buffimg[bfilled]=ppxy1;
                  buffimgQ[bfilled]=ppxy2;
                  buffimgU[bfilled]=ppxy3;
                  buffimgV[bfilled]=ppxy4;
                  buffpsf[bfilled++]=ppxy;
                }
               }
          } /* iteration end over convlution support +ve (-u,-v) */
        } 
      } /* end data rows loop */

 /* write last buffer */
 if (bfilled>0) {
   write_buffers_to_output_cpu_iquv(BL,t->Nx,t->Ny,t->writelock_img,t->writelock_psf,t->uvgrid,t->psfgrid,buffimg,buffpsf,buffoff,buffimgQ,buffimgU,buffimgV,rand()%MAX_GPU);
 }

 free(buffimg);
 free(buffpsf);
 free(buffoff);
 return NULL;
}


int
griddata(iodata *darr, int Nrows, float uvscale, float lambda, float deltaU, complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, float *wparr, int Nw, int *wpsupportX, int *wpsupportY, complex float *wkernel, float *uvxgrid, int Np, int Nt, int imgmode) {

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_grid_t *threaddata;
   pthread_mutex_t writelock_img,writelock_psf;

   /* divide data rows over threads */
   int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1,ci,Nthb;


   /* setup threads */
   pthread_mutex_init(&writelock_img, NULL);
   pthread_mutex_init(&writelock_psf, NULL);
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_grid_t*)malloc((size_t)Nt*sizeof(thread_data_grid_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }

   /* iterate over threads, allocating baselines per thread */
  ci=0;
  for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }
    threaddata[nth].darr=darr;
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].uvscale=uvscale;
    threaddata[nth].lambda=lambda;
    threaddata[nth].deltaU=deltaU;
    threaddata[nth].uvgrid=uvgrid;
    threaddata[nth].psfgrid=psfgrid;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].wparr=wparr;
    threaddata[nth].Nw=Nw;
    threaddata[nth].wpsupportX=wpsupportX;
    threaddata[nth].wpsupportY=wpsupportY;
    threaddata[nth].wkernel=wkernel;
    threaddata[nth].uvxgrid=uvxgrid;
    threaddata[nth].Np=Np;
    threaddata[nth].writelock_img=&writelock_img;
    threaddata[nth].writelock_psf=&writelock_psf;
    
    //printf("thread %d work on %d rows, starting from %d\n",nth, Nthb, ci);
    if (imgmode==IMG_IQUV0||imgmode==IMG_IQUV) {
     /* also grid Q,U,V */
     pthread_create(&th_array[nth],&attr,gridder_threadfn_iquv,(void*)(&threaddata[nth]));
    } else {
     pthread_create(&th_array[nth],&attr,gridder_threadfn,(void*)(&threaddata[nth]));
    }
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

 pthread_attr_destroy(&attr);
 pthread_mutex_destroy(&writelock_img);
 pthread_mutex_destroy(&writelock_psf);

 free(th_array);
 free(threaddata);

 return 0;
}

/* worker thread function */
static void *
weight_threadfn(void *data) {
  thread_data_weight_t *t=(thread_data_weight_t*)data;
  for(long int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
     t->wgrid[ci]=cabsf(t->psfgrid[ci]);
  } 
  return NULL;
}

int
weightdata(complex float *uvgrid, complex float *psfgrid, float *wgrid, int Nx, int Ny, float robust, int Nt) {

   return 0; /* FIXME: why return here? */
   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_weight_t *threaddata;

   /* divide data rows over threads */
   long int Nrows=Nx*Ny;
   long int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1;
   long int ci,Nthb;


   /* setup threads */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_weight_t*)malloc((size_t)Nt*sizeof(thread_data_weight_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }

   /* iterate over threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }

    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].wgrid=wgrid;
    threaddata[nth].psfgrid=psfgrid;
    
    //printf("thread %d work on %ld rows, starting from %ld\n",nth, Nthb, ci);
    pthread_create(&th_array[nth],&attr,weight_threadfn,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

 pthread_attr_destroy(&attr);

 free(th_array);
 free(threaddata);

 return 0;
}
