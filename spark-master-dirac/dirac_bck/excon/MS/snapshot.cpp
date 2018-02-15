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

typedef struct thread_data_wplane_{
  iodata *darr;
  int startrow,Nrows;
  float xx,xy,xz,yy,yz;
} thread_data_wplane_t;

typedef struct thread_data_wsub_{
  iodata *darr;
  int startrow,Nrows;
  float a,b,wmax,wmin;
} thread_data_wsub_t;


static void*
projectw_threadfn(void *data) {
 thread_data_wplane_t *t=(thread_data_wplane_t*)data;
 for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
  if (!t->darr[i].flag) {
    t->xx+=t->darr[i].u*t->darr[i].u;
    t->xy+=t->darr[i].u*t->darr[i].v;
    t->xz+=t->darr[i].u*t->darr[i].w;
    t->yy+=t->darr[i].v*t->darr[i].v;
    t->yz+=t->darr[i].v*t->darr[i].w;
  }
 }
 return NULL;
}

static void*
subtractw_threadfn(void *data) {
 thread_data_wsub_t *t=(thread_data_wsub_t*)data;
 for(int i = t->startrow; i < t->startrow+t->Nrows; i++) {
  if (!t->darr[i].flag) {
    t->darr[i].w-=t->a*t->darr[i].u+t->b*t->darr[i].v;
    if (t->darr[i].w>t->wmax) {
     t->wmax=t->darr[i].w;
    } else if (t->darr[i].w<t->wmin) {
     t->wmin=t->darr[i].w;
    }
  }
 }
 return NULL;
}

/* w-snapshot, project w coords to best fitting plane 
 w= a*u+b*v+delta_w 
 also find the new w-limits
*/
int
project_wplane(iodata *darr, int Nrows, int Nt, double *wpa,double *wpb, float *wmax, float *wmin) {

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_wplane_t *threaddata;
   thread_data_wsub_t *threaddata1;

   /* divide data rows over threads */
   int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1,ci,Nthb;

   /* setup threads */
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_wplane_t*)malloc((size_t)Nt*sizeof(thread_data_wplane_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
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
    threaddata[nth].xx=0.0f;
    threaddata[nth].xy=0.0f;
    threaddata[nth].xz=0.0f;
    threaddata[nth].yy=0.0f;
    threaddata[nth].yz=0.0f;
    
    pthread_create(&th_array[nth],&attr,projectw_threadfn,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  double xx=0.0,xy=0.0,xz=0.0,yy=0.0,yz=0.0;

  /* now wait for threads to finish */
  /* also calculate the sum of weights */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
   xx+=(double)threaddata[nth1].xx;
   xy+=(double)threaddata[nth1].xy;
   xz+=(double)threaddata[nth1].xz;
   yy+=(double)threaddata[nth1].yy;
   yz+=(double)threaddata[nth1].yz;
  }
  free(threaddata);
  /* invert matrix eqn |xx xy| |a| = |xz|
                      |xy yy| |b|   |yz|
  */
  *wpa=(xz*yy-yz*xy)/(xx*yy-xy*xy);
  *wpb=(xz-xx*(*wpa))/xy;

  /* now subtract a*u+b*v from baseline w, also find min,max */
  if ((threaddata1=(thread_data_wsub_t*)malloc((size_t)Nt*sizeof(thread_data_wsub_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  ci=0;
  for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }
    threaddata1[nth].darr=darr;
    threaddata1[nth].startrow=ci;
    threaddata1[nth].Nrows=Nthb;
    threaddata1[nth].a=(float)*wpa;
    threaddata1[nth].b=(float)*wpb;
    threaddata1[nth].wmax=-1e9f;
    threaddata1[nth].wmin=1e9f;

    pthread_create(&th_array[nth],&attr,subtractw_threadfn,(void*)(&threaddata1[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  float awmin,awmax;
  awmin=1e9f;
  awmax=-1e9f;
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
   if (awmin>threaddata1[nth1].wmin) {
     awmin=threaddata1[nth1].wmin;
   } else if (awmax<threaddata1[nth1].wmax) {
     awmax=threaddata1[nth1].wmax;
   }
 }
 *wmin=awmin;
 *wmax=awmax;
 free(threaddata1);
 pthread_attr_destroy(&attr);
 free(th_array);

 return 0;
}
