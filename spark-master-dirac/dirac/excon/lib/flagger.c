/* 
 *
   Copyright (C) 2016 Sarod Yatawatta <sarod@users.sf.net>  
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <glib.h>

//#define DEBUG

typedef struct thread_data_flagger_t_ {
  int startrow,Nrows;
  iodata *data;
  float sums[8]; //sum to find mean (real,imag) of XX,XY,YX,YY
  float vars[4]; //variances |XX|^2,|XY|^2,|YX|^2,|YY|^2
  int ncount; //how many points are unflagged
} thread_data_flagger_t;



static void *
estimate_std(void *data) {
   thread_data_flagger_t *t=(thread_data_flagger_t *)data;
   int nrow;
   t->sums[0]=t->sums[1]=t->sums[2]=t->sums[3]=t->sums[4]=t->sums[5]=t->sums[6]=t->sums[7]=0.0f;
   t->vars[0]=t->vars[1]=t->vars[2]=t->vars[3]=0.0f;
   t->ncount=0;
   for (nrow=t->startrow; nrow<t->Nrows+t->startrow; nrow++) {
     if (!t->data[nrow].flag) {
     t->sums[0]+=t->data[nrow].xx.x;
     t->sums[1]+=t->data[nrow].xx.y;
     t->sums[2]+=t->data[nrow].xy.x;
     t->sums[3]+=t->data[nrow].xy.y;
     t->sums[4]+=t->data[nrow].yx.x;
     t->sums[5]+=t->data[nrow].yx.y;
     t->sums[6]+=t->data[nrow].yy.x;
     t->sums[7]+=t->data[nrow].yy.y;
     t->vars[0]+=t->data[nrow].xx.x*t->data[nrow].xx.x+t->data[nrow].xx.y*t->data[nrow].xx.y;
     t->vars[1]+=t->data[nrow].xy.x*t->data[nrow].xy.x+t->data[nrow].xy.y*t->data[nrow].xy.y;
     t->vars[2]+=t->data[nrow].yx.x*t->data[nrow].yx.x+t->data[nrow].yx.y*t->data[nrow].yx.y;
     t->vars[3]+=t->data[nrow].yy.x*t->data[nrow].yy.x+t->data[nrow].yy.y*t->data[nrow].yy.y;
     t->ncount++;
     }
   }
   return NULL;
}


static void *
flagger_th(void *data) {
   thread_data_flagger_t *t=(thread_data_flagger_t *)data;
   int nrow;
   for (nrow=t->startrow; nrow<t->Nrows+t->startrow; nrow++) {
     if (!t->data[nrow].flag) {
     if (t->vars[0]<=t->data[nrow].xx.x*t->data[nrow].xx.x+t->data[nrow].xx.y*t->data[nrow].xx.y
     || t->vars[1]<=t->data[nrow].xy.x*t->data[nrow].xy.x+t->data[nrow].xy.y*t->data[nrow].xy.y
     || t->vars[2]<=t->data[nrow].yx.x*t->data[nrow].yx.x+t->data[nrow].yx.y*t->data[nrow].yx.y
     || t->vars[3]<=t->data[nrow].yy.x*t->data[nrow].yy.x+t->data[nrow].yy.y*t->data[nrow].yy.y) {
      t->data[nrow].flag=1;
     }
     }
   }
   return NULL;
}



/* flag data: first calculate standard deviation of XX,XY,YX,YY
              and flag data points that are XXs,XYs,YXs,YYs sigmas away
   iodata: Nx1 data points */
int
flag_data_sigma(iodata *data,int N, float XXs, float XYs, float YXs, float YYs, int Nthreads) {


  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_flagger_t *threaddata;


  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

  if ((th_array=(pthread_t*)malloc((size_t)Nthreads*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_flagger_t*)malloc((size_t)Nthreads*sizeof(thread_data_flagger_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* split data to threads */
  int NrowsNt=(N+Nthreads-1)/Nthreads;

  int ci,nth,nth1,Nthb;


  ci=0;
  for (nth=0;nth<Nthreads && ci<N; nth++) {
    if (ci+NrowsNt<N) {
      Nthb=NrowsNt;
    } else {
      Nthb=N-ci;
    }
    
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].data=data;

    pthread_create(&th_array[nth],&attr,estimate_std,(void*)(&threaddata[nth]));
   
    ci=ci+NrowsNt;
  }

  float XXm_re,XXm_im,XYm_re,XYm_im,YXm_re,YXm_im,YYm_re,YYm_im;
  float XXv,XYv,YXv,YYv;
  float invData;
  XXm_re=XXm_im=XYm_re=XYm_im=YXm_re=YXm_im=YYm_re=YYm_im=0.0f;
  XXv=XYv=YXv=YYv=0.0f;
  int Ncount=0;
  for (nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
   /* calculate mean,variance */
   XXm_re+=threaddata[nth1].sums[0];
   XXm_im+=threaddata[nth1].sums[1];
   XYm_re+=threaddata[nth1].sums[2];
   XYm_im+=threaddata[nth1].sums[3];
   YXm_re+=threaddata[nth1].sums[4];
   YXm_im+=threaddata[nth1].sums[5];
   YYm_re+=threaddata[nth1].sums[6];
   YYm_im+=threaddata[nth1].sums[7];

   XXv+=threaddata[nth1].vars[0];
   XYv+=threaddata[nth1].vars[1];
   YXv+=threaddata[nth1].vars[2];
   YYv+=threaddata[nth1].vars[3];
   Ncount+=threaddata[nth1].ncount;
  }
  invData=(Ncount>0?1.0f/(float)Ncount:1.0f);
  XXm_re*=invData;
  XXm_im*=invData;
  XYm_re*=invData;
  XYm_im*=invData;
  YXm_re*=invData;
  YXm_im*=invData;
  YYm_re*=invData;
  YYm_im*=invData;
  XXv*=invData;
  XYv*=invData;
  YXv*=invData;
  YYv*=invData;
#ifdef DEBUG
  printf("mean %f,%f,%f,%f %f,%f,%f,%f\n",XXm_re,XYm_re,YXm_re,YYm_re,XXm_im,XYm_im,YXm_im,YYm_im);
  printf("var %f,%f,%f,%f\n",XXv,XYv,YXv,YYv);
#endif
  /* find standard deviations */
  XXv=XXv-(XXm_re*XXm_re+XXm_im*XXm_im);
  XYv=XYv-(XYm_re*XYm_re+XYm_im*XYm_im);
  YXv=YXv-(YXm_re*YXm_re+YXm_im*YXm_im);
  YYv=YYv-(YYm_re*YYm_re+YYm_im*YYm_im);

#ifdef DEBUG
  printf("std^2 %f,%f,%f,%f\n",XXv,XYv,YXv,YYv);
#endif

  /* update the clipping value^2 = (factor)^2 * std^2 */
  XXv *=XXs*XXs;
  XYv *=XYs*XYs;
  YXv *=YXs*YXs;
  YYv *=YYs*YYs;
  for (nth1=0; nth1<nth; nth1++) {
     threaddata[nth1].vars[0]=XXv;    
     threaddata[nth1].vars[1]=XYv;    
     threaddata[nth1].vars[2]=YXv;    
     threaddata[nth1].vars[3]=YYv;    
     /* flag data */
     pthread_create(&th_array[nth1],&attr,flagger_th,(void*)(&threaddata[nth1]));
  }

  for (nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  free(threaddata);
  pthread_attr_destroy(&attr);
  free(th_array);

  return 0;
}
