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

typedef struct thread_data_griduv_t_ {
  iodata *darr;
  int startrow,Nrows;
  float lambda;
  float uvscale;
  float *wgrid;
  int Nx,Ny;
  float sumwt;
  float f2,d2;
  pthread_mutex_t *writelock_w;
} thread_data_griduv_t;

static void *
weighteruv_threadfn(void *data) {
   thread_data_griduv_t *t=(thread_data_griduv_t*)data;
   unsigned long int *buffoff;
   int BL=DATA_BUF_LEN;
   if ((buffoff=(unsigned long int*)calloc((size_t)BL,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   int bfilled=0;

   
/*************************************************/
   for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
      t->darr[i].wt=0.0f;
      if (!t->darr[i].flag) {
          float tempu = t->darr[i].u;
          float tempv = t->darr[i].v;
         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*t->uvscale;
          float vf=tempv*t->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)t->Nx); /*  width */
          float vi=(vf*(float)t->Ny); /*  width */
          float poffX=-1.0f;
          float poffY=-1.0f;
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*t->Nx+poffX);
          int y=(int)round(-vi+0.5f*t->Ny+poffY);
          unsigned long int coffset=0;
          /* no FFTshift */
          if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
            coffset=ifftshift_index(t->Nx,x)*(t->Ny)+ifftshift_index(t->Ny,y);
            if (bfilled<BL) {
               buffoff[bfilled++]=coffset;
            } else {
              /* empty buffer */
              pthread_mutex_lock(t->writelock_w);
              for (int bidx=0; bidx<BL; bidx++) {
                t->wgrid[buffoff[bidx]]+=1.0f;
              }
              pthread_mutex_unlock(t->writelock_w);

              bfilled=0;
              buffoff[bfilled++]=coffset;
            }
            t->sumwt+=1.0f; /* increment total weight */
          }
          /* also add (-uf,-vf) : mirror image */
          x=(int)round(-ui+0.5f*t->Nx+poffX);
          y=(int)round(vi+0.5f*t->Ny+poffY);

          if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
            coffset=ifftshift_index(t->Nx,x)*(t->Ny)+ifftshift_index(t->Ny,y);
            if (bfilled<BL) {
               buffoff[bfilled++]=coffset;
            } else {
              /* empty buffer */
              pthread_mutex_lock(t->writelock_w);
              for (int bidx=0; bidx<BL; bidx++) {
                t->wgrid[buffoff[bidx]]+=1.0f;
              }
              pthread_mutex_unlock(t->writelock_w);
              bfilled=0;
              buffoff[bfilled++]=coffset;
            }
            t->sumwt+=1.0f; /* increment total weight */
          }
        }
     }
     /* write last buffer */
     if (bfilled>0) {
       /* empty buffer */
       pthread_mutex_lock(t->writelock_w);
       for (int bidx=0; bidx<bfilled; bidx++) {
            t->wgrid[buffoff[bidx]]+=1.0f;
       }
       pthread_mutex_unlock(t->writelock_w);
     }
 free(buffoff);
 return NULL;
}


static void *
weighteruvupdate_threadfn(void *data) {
   thread_data_griduv_t *t=(thread_data_griduv_t*)data;
   
/*************************************************/
   for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
      if (!t->darr[i].flag) {
          float tempu = t->darr[i].u;
          float tempv = t->darr[i].v;
         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*t->uvscale;
          float vf=tempv*t->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)t->Nx); /*  width */
          float vi=(vf*(float)t->Ny); /*  width */
          float poffX=-1.0f;
          float poffY=-1.0f;
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*t->Nx+poffX);
          int y=(int)round(-vi+0.5f*t->Ny+poffY);
          unsigned long int coffset=0;

          /* only the +ve u,v considered */
          if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
            coffset=ifftshift_index(t->Nx,x)*(t->Ny)+ifftshift_index(t->Ny,y);
            if (t->wgrid[coffset]>0.0f) {
             t->darr[i].wt=1.0f/(t->wgrid[coffset]*t->f2+t->d2);
            } 
          } 
        }
     }
 return NULL;
}


/* natural weights */
static void *
weighteruvupdate_natural_threadfn(void *data) {
   thread_data_griduv_t *t=(thread_data_griduv_t*)data;
   
/*************************************************/
   for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
      if (!t->darr[i].flag) {
          t->darr[i].wt=1.0f;
        }
     }
 return NULL;
}




typedef struct thread_data_weightsqr_t_ {
  long int startrow,Nrows;
  float *wgrid;
  float sumwsquare;
} thread_data_weightsqr_t;


static void *
zero_threadfn(void *data) {
  thread_data_weightsqr_t *t=(thread_data_weightsqr_t*)data;
  for(long int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
    t->wgrid[ci]=0.0f;
  }
  return NULL;
}

static void *
square_threadfn(void *data) {
  thread_data_weightsqr_t *t=(thread_data_weightsqr_t*)data;
  for(long int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
    t->sumwsquare+=t->wgrid[ci]*t->wgrid[ci];
  }
  return NULL;
}

int
weightuvdata(iodata *darr, int Nrows, float uvscale, float lambda, float robustval, float *wgrid, int Nx, int Ny, int Nt, int weightmode, float *sumweights, float *sumweights2) {

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_griduv_t *threaddata;
   thread_data_weightsqr_t *threaddata_sq;
   pthread_mutex_t writelock_w;

   /* divide data rows over threads */
   int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1,ci,Nthb;

   /* divide data rows over threads */
   long int NrowsGrid=Nx*Ny;
   long int NrowsNtGrid=(NrowsGrid+Nt-1)/Nt;
   long int ciGrid,NthbGrid;


   /* setup threads */
   pthread_mutex_init(&writelock_w, NULL);
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_griduv_t*)malloc((size_t)Nt*sizeof(thread_data_griduv_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((threaddata_sq=(thread_data_weightsqr_t*)malloc((size_t)Nt*sizeof(thread_data_weightsqr_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   /******************** zero weight array ***************************/
   ciGrid=0;
   for (nth=0;  nth<Nt && ciGrid<Nrows; nth++) {
    if (ciGrid+NrowsNtGrid<NrowsGrid) {
     NthbGrid=NrowsNtGrid;
    } else {
     NthbGrid=NrowsGrid-ciGrid;
    }

    threaddata_sq[nth].startrow=ciGrid;
    threaddata_sq[nth].Nrows=NthbGrid;
    threaddata_sq[nth].wgrid=wgrid;
    threaddata_sq[nth].sumwsquare=0.0f;

    pthread_create(&th_array[nth],&attr,zero_threadfn,(void*)(&threaddata_sq[nth]));
    /* next baseline set */
    ciGrid=ciGrid+NrowsNtGrid;
  }

  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

   /******************** calculate weight array ***************************/
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
    threaddata[nth].wgrid=wgrid;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].sumwt=0.0f;
    threaddata[nth].writelock_w=&writelock_w;
    
    pthread_create(&th_array[nth],&attr,weighteruv_threadfn,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* calculate sum (weight)^2 */

  /* now wait for threads to finish */
  /* also calculate the sum of weights */
  float sumwt=0.0f;
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
   sumwt+=threaddata[nth1].sumwt;
  }

  printf("Sum of weights=%f\n",sumwt);

   /******************** sum (weight)^2 array ***************************/
   ciGrid=0;
   for (nth=0;  nth<Nt && ciGrid<Nrows; nth++) {
    if (ciGrid+NrowsNtGrid<NrowsGrid) {
     NthbGrid=NrowsNtGrid;
    } else {
     NthbGrid=NrowsGrid-ciGrid;
    }

    threaddata_sq[nth].startrow=ciGrid;
    threaddata_sq[nth].Nrows=NthbGrid;
    threaddata_sq[nth].wgrid=wgrid;
    threaddata_sq[nth].sumwsquare=0.0f;

    pthread_create(&th_array[nth],&attr,square_threadfn,(void*)(&threaddata_sq[nth]));
    /* next baseline set */
    ciGrid=ciGrid+NrowsNtGrid;
  }

  /* now wait for threads to finish */
  /* also calculate sum of weight^2 */
  float sumwt2=0.0f;
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
   sumwt2+=threaddata_sq[nth1].sumwsquare;
  }

  printf("Sum of weights^2=%f\n",sumwt2);
  *sumweights=sumwt;
  *sumweights2=sumwt2;

   /******************** update uv weights ***************************/
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
    threaddata[nth].wgrid=wgrid;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    if (weightmode==WMODE_UNIFORM) {
     threaddata[nth].f2=1.0f;
     threaddata[nth].d2=0.0f;
    } else if (weightmode==WMODE_ROBUST) {
     threaddata[nth].f2=robustval/(sumwt2/sumwt);
     threaddata[nth].d2=1.0f;
    } else {
     threaddata[nth].f2=1.0f;
     threaddata[nth].d2=0.0f;
    }
    threaddata[nth].writelock_w=&writelock_w;
    
    if (weightmode==WMODE_NATURAL) {
     pthread_create(&th_array[nth],&attr,weighteruvupdate_natural_threadfn,(void*)(&threaddata[nth]));
    } else {
     pthread_create(&th_array[nth],&attr,weighteruvupdate_threadfn,(void*)(&threaddata[nth]));
    }
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* now wait for threads to finish */
  /* also calculate the sum of weights */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }


 pthread_attr_destroy(&attr);
 pthread_mutex_destroy(&writelock_w);

 free(th_array);
 free(threaddata);
 free(threaddata_sq);

 return 0;
}
