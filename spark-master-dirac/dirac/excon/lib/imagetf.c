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


#include <wcs.h>
#include <math.h>
#include <pthread.h>

#include "largefft.h"


typedef struct thread_data_lmwcs_t_ {
  int startrow,Nrows;
  float *din;
  float *dout;
  struct wcsprm* mwcs;
  int Nx,Ny;
  double a,b;
} thread_data_lmwcs_t;


/* transform image pixels
  l' = l +a(sqrt(1-l^2-m^2)-1)
  m' = m +b(sqrt(1-l^2-m^2)-1)
  din is in l',m' grid, with pixel values (column major order)
  dout is l,m grid, with no pixel value
  update dout pixel values using din pixel values and right coords
  both din,dout size Nx x Ny

  delta_l: pixel size (rad)
  ra0,dec0: phase center (rad)
  Nt: no. threads
*/

static void *
imtrans_threadfn(void *data) {
  thread_data_lmwcs_t *t=(thread_data_lmwcs_t*)data;
  int ci,cj;

 /* coordinate conversion */
 double *pixelc, *imgc, *worldc, *phic, *thetac;
 int *statc,status;
 double myl,mym,ll,mm;

 int ncoord=1;
 if ((pixelc=(double*)calloc((size_t)ncoord*2,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((imgc=(double*)calloc((size_t)ncoord*2,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((worldc=(double*)calloc((size_t)ncoord*2,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((phic=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((thetac=(double*)calloc((size_t)ncoord,sizeof(double)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }
 if ((statc=(int*)calloc((size_t)ncoord,sizeof(int)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
 }

  unsigned long int NxNy=t->Nx*t->Ny;
  unsigned long int poffin;
  unsigned long int poffout;
  /* iterate over pixels of the output image, 
    and read right value from right image */
  for (cj=0; cj<t->Ny; cj++) {
   for (ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
       pixelc[0]=(double)(ci+1); /* need to do this here because we update it later */
       pixelc[1]=(double)(cj+1);
       //printf("========== pixel %lf,%lf\n",pixelc[0],pixelc[1]);
       status=wcsp2s(t->mwcs, ncoord, 2, pixelc, imgc, phic, thetac, worldc, statc); 

     //printf("%d x,y %lf %lf -> l,m %lf %lf\n",status,pixelc[0],pixelc[1],imgc[0],imgc[1]);
     if (!statc[0] && !status) {
       myl=(imgc[0])*CPI_180; mym=(imgc[1])*CPI_180;
       /* now calculate l',m' in old image */
       /* no need to check if -ve because wcs status will be >0 later */
       double tmps=(sqrt(1.0-myl*myl-mym*mym)-1.0);
       ll=myl+t->a*tmps;
       mm=mym+t->b*tmps;
      imgc[0]=ll*C180_PI;
      imgc[1]=mm*C180_PI;
      status=celx2s(&(t->mwcs->cel), 1,1,1,1, &imgc[0], &imgc[1], phic, thetac, &worldc[0], &worldc[1], statc);

      //printf("%d ll,mm %lf %lf -> ra,dec %lf %lf\n",status, imgc[0],imgc[1],worldc[0],worldc[1]);
      status=wcss2p(t->mwcs, ncoord, 2, worldc, phic, thetac, imgc, pixelc, statc); 
      //printf("%d ra,dec %lf %lf -> x,y %lf %lf\n",status, worldc[0],worldc[1],pixelc[0],pixelc[1]);
    
      if (!status) {
       /* update this pixel with value of orig image */
       poffin=t->Nx*((int)round(pixelc[1])-1)+(int)round(pixelc[0])-1;
       poffout=t->Nx*(cj)+ci;
       if (poffin<NxNy && poffout<NxNy) {
        t->dout[poffout]=t->din[poffin]; 
       } else {
      //printf("Error %ld->%ld pix0=%d,%d, lm=%lf,%lf, lm'=%lf,%lf, pix'=%lf,%lf\n",poffin,poffout,ci+1,cj+1,myl,mym,ll,mm,pixelc[0],pixelc[1]);
        t->dout[poffout]=0.0f;
       }
      }
     }
    }
  }
  free(pixelc);
  free(imgc);
  free(worldc);
  free(phic);
  free(thetac);
  free(statc);

  return NULL;
}



int
do_image_lm_transform(float *din,float *dout,int Nx,int Ny,double a,double b,double delta_l,double ra0, double dec0, int Nt) {
  int ci;
  struct wcsprm *mwcs;
  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_lmwcs_t *threaddata;

  /* divide image rows over threads */
  int NrowsNt=(Nx+Nt-1)/Nt;
  int nth,nth1,Nthb;


  if ((mwcs=(struct wcsprm*)calloc((size_t)Nt,sizeof(struct wcsprm)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  for (ci=0; ci<Nt; ci++) {
    generate_def_wcs(&mwcs[ci],delta_l*C180_PI,(double)(Nx/2+1),ra0*C180_PI,dec0*C180_PI,0);
  }

  /* setup threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

  if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_lmwcs_t*)malloc((size_t)Nt*sizeof(thread_data_lmwcs_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* iterate over threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Nx; nth++) {
    if (ci+NrowsNt<Nx) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nx-ci;
    }
    threaddata[nth].din=din;
    threaddata[nth].dout=dout;
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].mwcs=&mwcs[nth];
    threaddata[nth].a=a;
    threaddata[nth].b=b;
    //printf("thread %d work on %d rows, starting from %d\n",nth, Nthb, ci);
    pthread_create(&th_array[nth],&attr,imtrans_threadfn,(void*)(&threaddata[nth]));
    /* next row set */
    ci=ci+NrowsNt;
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  free(threaddata);
  free(th_array);

  for (ci=0; ci<Nt; ci++) {
    wcsfree(&mwcs[ci]);
  }
  free(mwcs);
  return 0;
}
