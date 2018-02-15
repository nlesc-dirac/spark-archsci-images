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

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "largefft.h"
#include "gridder.h"

#ifdef HAVE_CUDA
#include <cuda.h>
#include <cuComplex.h>
#endif /* HAVE_CUDA */

typedef struct thread_data_wcs_t_ {
  int startrow,Nrows;
  double *lmgrid;
  double *denom;
  struct wcsprm* mwcs;
  int Np;
  int snapshot;
  double wpa,wpb;
} thread_data_wcs_t;

typedef struct thread_data_lmpswf_t_ {
  int startrow,Nrows;
  int Np,Npad,Nz0;
  complex float *din;
  float *pswfxy;
  float w;
  double *lmgrid;
  double *denom;
} thread_data_lmpswf_t;


static void *
lmpswf_threadfn(void *data) {
  thread_data_lmpswf_t *t=(thread_data_lmpswf_t*)data;

  float phaseterm;
  float cosp,sinp;
  int ci,nrow;
  for (nrow=t->startrow; nrow<t->startrow+t->Nrows; nrow++) {
    for (ci=0; ci<t->Np; ci++) { 
      phaseterm=(float)(t->lmgrid[(nrow-t->Npad)*t->Np+ci]*t->w);
      sincosf(phaseterm,&sinp,&cosp);
      float invdenom=t->pswfxy[nrow*t->Nz0+t->Npad+ci]/(float)t->denom[(nrow-t->Npad)*t->Np+ci];
      t->din[nrow*t->Nz0+t->Npad+ci]=cosp*invdenom+_Complex_I*sinp*invdenom;
     }
  } 

  return NULL;
}

static void *
wcsgrid_threadfn(void *data) {
 thread_data_wcs_t *t=(thread_data_wcs_t*)data;

 /* coordinate conversion */
 double *pixelc, *imgc, *worldc, *phic, *thetac;
 int *statc,status;

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

// wcsprt(t->mwcs);
 double myl,mym;
 int nrow,ci;
 if (t->snapshot) {
 for (nrow=t->startrow; nrow<t->Nrows+t->startrow; nrow++) {
    pixelc[0]=(double)(nrow+1);
    for(ci=0; ci<t->Np; ci++) {
      pixelc[1]=(double)(ci+1);
      
     status = wcsp2s(t->mwcs, ncoord, t->mwcs->naxis, pixelc, imgc, phic, thetac, worldc, statc);
     if (!statc[0] && !status) {
       myl=(imgc[0])*CPI_180; mym=(imgc[1])*CPI_180;
       /* if snapshot mode, l'<=l+a(sqrt(1-l^2-m^2)-1) m'<=m+b(sqrt(1-l^2-m^2)-1) */
       double tmpa=1.0-myl*myl-mym*mym;
       /* no need to check tmpa>0 because wcs status would have indicated this*/
       tmpa=sqrt(tmpa)-1.0;
       double myl1,mym1;
       myl1=myl+t->wpa*tmpa;
       mym1=mym+t->wpb*tmpa;
       /* sqrt(1-l'^2-m'^2)-1 */
       double ltr=1.0-myl1*myl1-mym1*mym1;
       if (ltr>0.0) {
       t->lmgrid[nrow*t->Np+ci]=(sqrt(ltr)-1.0);
       /* sqrt(1-l^2-m^2) -a l - b m */
       t->denom[nrow*t->Np+ci]=(tmpa+1.0-t->wpa*myl-t->wpb*mym);
       } else {
        t->lmgrid[nrow*t->Np+ci]=-1.0;
        t->denom[nrow*t->Np+ci]=1.0;
       }
      } else { 
     // else lmgrid[nrow*Np+ci]=0.0; which is automatically done using calloc()
       t->denom[nrow*t->Np+ci]=1.0;
      }
     }
    }
 } else {
 for (nrow=t->startrow; nrow<t->Nrows+t->startrow; nrow++) {
    pixelc[0]=(double)(nrow+1);
    for(ci=0; ci<t->Np; ci++) {
      pixelc[1]=(double)(ci+1);
      
     status = wcsp2s(t->mwcs, ncoord, t->mwcs->naxis, pixelc, imgc, phic, thetac, worldc, statc);
     if (!statc[0] && !status) {
       myl=(imgc[0])*CPI_180; mym=(imgc[1])*CPI_180;
       double tmpa=1.0-myl*myl-mym*mym;
       /* no need to check tmpa>0 because wcs status would have indicated this*/
       tmpa=sqrt(tmpa)-1.0;
       /* sqrt(1-l^2-m^2)-1 */
       t->lmgrid[nrow*t->Np+ci]=tmpa;
        /* sqrt(1-l^2-m^2) */
       t->denom[nrow*t->Np+ci]=tmpa+1.0;
     } else {
     // else lmgrid[nrow*Np+ci]=0.0; which is automatically done using calloc()
      t->denom[nrow*t->Np+ci]=1.0;
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

/* 
 1)generate l,m grid size NpxNp in image FOV
 2)generate PSWF grid pswf[] in [-0.5,0.5]x[-0.5,0.5] for the pixel grid
 3)multiply pswf[] with 
   e^{-2 pi j(w (sqrt(1-r^2) -1))}
   r^2=l^2+m^2 : l,m in l,m grid of image (or support of PSWF)
   w : taken from each element of wparr[]
 4)zero pad multiplied value, FFTshift and then do complex to complex FFT
 5)FFTshift back the transformed value
 6)copy back the original length NpxNp as output array
 
  wparr: Nwx1 array of W values for w projection
  
  wkernel: (NpxNp)*Nw x 1 array of complex float w kernels 

  Nz: zero padded length Nz>Np
  Np: actual no of points which samples range [0,1]
  M: convolution support in uv pixels, for the bandwidth of PSWF
  delta_l,delta_m: pixel width for conv. kernel evaluation (radians), delta_l*Np=FOV
  tol: cutoff to determine support length of each kernel
  supportX,Y: Nwx1 support length for each w plane
  maxsupport: returns the max support value

  ra0,dec0: coordinates of reference pixel (rad)
  Nt: no. of threads
*/ 
int
generate_2d_wconv_kernels(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel, float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot, int Nt) {
  /* zero pixel
    <---------------> Np/2+1 values
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero  for even Np : each half Np/2 pixels
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero for odd Np : left Np/2 right Np/2+1 pixels
  */
  int ci,nw;
  /* arrays to calculate PSWF values in [0,1] 
     [0,1] is sampled by Np values, the values beyond that upto Nz
     are additions
  */
  float *pswf_x,*pswf_y,*pswf_yl,*pswfxy;
  float pswf_c=M_PI*0.5f*(float)M;
  int Np2,nrow;
  /* get right half length for even/odd Np */
  if(Np%2) { /* odd */
   Np2=Np/2;
  } else { /* even : note right side is shorter as 0 pixel is taken from left */
   Np2=Np/2-1;
  }

  int Nl=Np/2+1;
//printf("Np=%d Np2=%d Nl=%d\n",Np,Np2,Nl);
  if ((pswf_x=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_y=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_yl=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }


  float scaleval;
  float pswf_delta=1.0f/((float)(Nl-1)); /* 1,2,..Nl : Nl-1 intervals */
  pswf_x[0]=0.0f;//+0.5*pswf_delta; FIXME why 0.5pixel?
  for (ci=1; ci<Nl; ci++) {
      pswf_x[ci]=pswf_x[ci-1]+pswf_delta;
  }
  pswf_nn(pswf_c,Nl,pswf_x,pswf_y);

  /* copy right hand side to left hand side*/
  for (ci=0; ci<Nl; ci++) {
    pswf_yl[ci]=pswf_y[Np/2-ci];
  }

  //for(ci=0; ci<Nl; ci++) {
  //printf("%d %f %f\n",ci,pswf_y[ci],pswf_yl[ci]);
  //}
  //exit(1);
  
  int Npad=(Nz-Np)/2; /* same for even/odd */
  int Nz0=Np+2*Npad;
  /* now allocate zero padded memory for 2D PSWF, row major order */
  if ((pswfxy=(float*)calloc((size_t)Nz0*Nz0,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* pswfxy is zero, valid range for copying
    |000000000000..............000000000000|
      (Nz-Np)/2  ^             ^  (Nz-Np)/2
      =Npad           =Np          =Npad
                 0...........Np-1 
    zero pixel:
                      ^Npad+Np/2 : even Np ,left (Npad+Np/2) right (Npad+Np/2)
                      ^Npad+Np/2 : odd Np, left (Npad+Np/2) right (Npad+Np/2+1)
  */
  for (nrow=Npad; nrow<Npad+Np; nrow++) {
    /* copy to right hand side Np2 pixels,omit 0 pixel */
    memcpy(&pswfxy[nrow*Nz0+Npad+Np/2+1],&pswf_y[1],(size_t)Np2*sizeof(float));
    /* copy to left hand side Np/2+1 pixels */
    memcpy(&pswfxy[nrow*Nz0+Npad],pswf_yl,(size_t)(Np/2+1)*sizeof(float));
    /* now scale this row */
    ci=nrow-Npad;
    if (ci<Nl) {
     scaleval=pswf_yl[ci];
    } else {
     ci=ci-Nl;
     scaleval=pswf_y[ci+1];
    }
    my_fscal(Np,scaleval,&pswfxy[nrow*Nz0+Npad]);
    //printf("row %d ci %d scale %f\n",nrow,ci,scaleval);
  }
  free(pswf_x);
  free(pswf_y);
  free(pswf_yl);
/*  for (ci=0; ci<Nz0*Nz0; ci++) {
    printf("%f\n",pswfxy[ci]);
  }  */
/*  for (nrow=Npad; nrow<Npad+Np; nrow++) {
   for (ci=Npad; ci<Np+Npad; ci++) {
    printf("%f %f\n",(pswfxy[nrow*Nz+ci]),0.0f);
   }
  } */ 


  struct wcsprm *mwcs;
  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_wcs_t *threaddata;
  /* divide data rows over threads */
  int NrowsNt=(Np+Nt-1)/Nt;
  int nth,nth1,Nthb;

  if ((mwcs=(struct wcsprm*)calloc((size_t)Nt,sizeof(struct wcsprm)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  for (ci=0; ci<Nt; ci++) {
    generate_def_wcs(&mwcs[ci],(double)delta_l*C180_PI,(double)Np*0.5,ra0*C180_PI,dec0*C180_PI,0);
  }
  /* setup threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

  if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_wcs_t*)malloc((size_t)Nt*sizeof(thread_data_wcs_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* generate l,m grid NpxNp 
    (sqrt(1-r^2) -1)
    r^2=l^2+m^2  
    Np*delta_l = image FOV = Nx*delta_l0, delta_l0: image pixel width 
    so delta_l=Nx*delta_l0/Np

    if snapshot mode,
    l<= l+a*(sqrt(1-r^2)-1)
    m<= m+b*(sqrt(1-r^2)-1)
 */
  double *lmgrid,*denom;
  if ((lmgrid=(double*)calloc((size_t)Np*Np,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((denom=(double*)calloc((size_t)Np*Np,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* iterate over threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Np; nth++) {
    if (ci+NrowsNt<Np) {
     Nthb=NrowsNt;
    } else {
     Nthb=Np-ci;
    }
    threaddata[nth].lmgrid=lmgrid;
    threaddata[nth].denom=denom;
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].Np=Np;
    threaddata[nth].mwcs=&mwcs[nth];
    threaddata[nth].snapshot=snapshot;
    threaddata[nth].wpa=wpa;
    threaddata[nth].wpb=wpb;
    //printf("thread %d work on %d rows, starting from %d\n",nth, Nthb, ci);
    pthread_create(&th_array[nth],&attr,wcsgrid_threadfn,(void*)(&threaddata[nth]));
    /* next row set */
    ci=ci+NrowsNt;
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  free(threaddata);

  for (ci=0; ci<Nt; ci++) {
   wcsfree(&mwcs[ci]);
  }
  free(mwcs);

  /* scale by -2pi */
  my_dscal(Np*Np,-2.0*M_PI,lmgrid);


/*  for (ci=0; ci<Np*Np; ci++) {
    printf("%lf\n",lmgrid[ci]);
  }  
*/
  complex float *din;
  if ((din=(complex float*)calloc((size_t)Nz0*Nz0,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  fftwf_plan plan0;
  do_create_fftw_plan_2d(&plan0, din, din, Nz0, Nz0, Nt);

  thread_data_lmpswf_t *threaddataN;
  if ((threaddataN=(thread_data_lmpswf_t*)malloc((size_t)Nt*sizeof(thread_data_lmpswf_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }

/* DEBUG WPROJ */
  FILE *bfilep=fopen("debug_w_in","w+");
  fprintf(bfilep,"%d %d\n",Np,Nw);

  /* form the product pswf[] x exp{j2 pi w (...) } */
  for (nw=0; nw<Nw;nw++) {
/*  for (ci=0; ci<Nz0*Nz0; ci++) {
    printf("%f %f\n",crealf(din[ci]),cimagf(din[ci]));
  }  */

  ci=0;
  for (nth=0;  nth<Nt && ci<Np; nth++) {
    if (ci+NrowsNt<Np) {
     Nthb=NrowsNt;
    } else {
     Nthb=Np-ci;
    }
    threaddataN[nth].w=wparr[nw];
    threaddataN[nth].lmgrid=lmgrid;
    threaddataN[nth].denom=denom;
    threaddataN[nth].pswfxy=pswfxy;
    threaddataN[nth].din=din;
    threaddataN[nth].startrow=ci+Npad;
    threaddataN[nth].Nrows=Nthb;
    threaddataN[nth].Np=Np;
    threaddataN[nth].Npad=Npad;
    threaddataN[nth].Nz0=Nz0;
    pthread_create(&th_array[nth],&attr,lmpswf_threadfn,(void*)(&threaddataN[nth]));
    /* next row set */
    ci=ci+NrowsNt;
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }
/* DEBUG WPROJ */
  for (nrow=Npad; nrow<Npad+Np; nrow++) {
   for (ci=Npad; ci<Npad+Np; ci++) {
    fprintf(bfilep,"%f %f\n",crealf(din[nrow*Nz0+ci]),cimagf(din[nrow*Nz0+ci]));
   }
  }

    /* do fftshift */
   do_fftshift_inplace(din,Nz0,Npad+Np/2);

   /* do fft */
   fftwf_execute(plan0);
   do_ifftshift_inplace(din,Nz0,Npad+Np/2);
   
  /* ignore padded values */
   /* copy back to the kernel array */
   for (nrow=Npad; nrow<Npad+Np; nrow++) {
    memcpy(&wkernel[nw*Np*Np+(nrow-Npad)*Np],&din[nrow*Nz0+Npad],(size_t)Np*sizeof(complex float));
   }
   /* reset din to zero */
   memset(din,0,(size_t)Nz0*Nz0*sizeof(complex float));

  }

/* DEBUG WPROJ */
  fclose(bfilep);

  free(threaddataN);
  pthread_attr_destroy(&attr);
  free(th_array);

  do_destroy_fftw_plan(plan0);
  free(pswfxy);
  free(din);
  free(lmgrid);
  free(denom);

  /* scale such that zero pixel of 1st w plane is 1 */
  float abs_re;
  abs_re=cabsf(wkernel[(Np/2)*Np+Np/2]); 
  my_csscal(Np*Np*Nw,1.0f/abs_re,wkernel);

/* DEBUG WPROJ */
  FILE *cfilep=fopen("debug_w","w+");
  fprintf(cfilep,"%d %d\n",Np,Nw);
  for (ci=0; ci<Np*Np*Nw; ci++) {
    fprintf(cfilep,"%f %f\n",crealf(wkernel[ci]),cimagf(wkernel[ci]));
  }
  fclose(cfilep);


  float tol0;
  *maxsupport=0;
  /* scan the zero pixel row to find support in X direction */
  /* scan the zero pixel column to find support in Y direction */
  /* support determined by edge/center > tol */
  for (nw=0; nw<Nw; nw++) {
    supportX[nw]=0;
    tol0=tol*cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2]);
    for (ci=Np/2; ci<Np; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+ci]);
      if (tol0>abs_re) {
        break;
      } else {
        supportX[nw]++;
      }
    }

    supportY[nw]=0;
    for (ci=0; ci<Np/2; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2+ci*Np]);
      if (tol0>abs_re) {
        break;
      } else {
        supportY[nw]++;
      }
    }


    /* scale back support to original pixels (without padding) */
    abs_re=(float)(supportX[nw]*Np)/(float)Nz;
    supportX[nw]=(int)(abs_re*2.0f)+1;
    if (supportX[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->X  plane %d is %d, exceeds %d\n",nw,supportX[nw],Np);
     supportX[nw]=Np;
    }
    abs_re=(float)(supportY[nw]*Np)/(float)Nz;
    supportY[nw]=(int)(abs_re*2.0f)+1;
    if (supportY[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->Y plane %d is %d, exceeds %d\n",nw,supportY[nw],Np);
     supportY[nw]=Np;
    }
    /* also find max support */
    if (supportX[nw]>*maxsupport) {
      *maxsupport=supportX[nw];
    }
  }
  printf("Maximum convolution kernel support is %d pixels but allowed is %d.\n",*maxsupport*Nz/Np,Np);

  printf("[");
  for(nw=0; nw<Nw; nw++) { 
   printf("(%d %d %d %f) ",nw,supportX[nw],supportY[nw],wparr[nw]);
  } 
  printf("]\n");

  /* DEBUG WPROJ */
  /* the w-planes where support length has exceeded, oversample 
    and recalculate the FFT */



  return 0;
}



/* oversample when generating w kernels, but still select Np pixels
   r: oversample ratio >1.0
   both zero padding and lm-pswf evaluation is oversampled by this ratio (Np x r)
   but after FFT, Np pixels again copied back */
int
generate_2d_wconv_kernels_over(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel, float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot,  double r, int Nt) {
  /* zero pixel
    <---------------> Np/2+1 values
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero  for even Np : each half Np/2 pixels
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero for odd Np : left Np/2 right Np/2+1 pixels
  */
  int ci,nw;
  /* arrays to calculate PSWF values in [0,1] 
     [0,1] is sampled by Np values, the values beyond that upto Nz
     are additions
  */
  float *pswf_x,*pswf_y,*pswf_yl,*pswfxy;
  float pswf_c=M_PI*0.5f*(float)M;
  int Np2,nrow;
  /* oversample Np */
  int Npad1=(int)((r-1.0)*(double)Np*0.5);
  int Npd=Np+2*Npad1;
  /* oversample Nz */
  int Nzd=(int)(r*(double)Nz);
  /* undersample delta_l */
  float delta_ld=delta_l/(float)r;

  /* get right half length for even/odd Np */
  if(Npd%2) { /* odd */
   Np2=Npd/2;
  } else { /* even : note right side is shorter as 0 pixel is taken from left */
   Np2=Npd/2-1;
  }

  int Nl=Npd/2+1;
  if ((pswf_x=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_y=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_yl=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }


  float scaleval;
  float pswf_delta=1.0f/((float)(Nl-1)); /* 1,2,..Nl : Nl-1 intervals */
  pswf_x[0]=0.0f;//+0.5*pswf_delta; FIXME why 0.5pixel?
  for (ci=1; ci<Nl; ci++) {
      pswf_x[ci]=pswf_x[ci-1]+pswf_delta;
  }
  pswf_nn(pswf_c,Nl,pswf_x,pswf_y);

  /* copy right hand side to left hand side*/
  for (ci=0; ci<Nl; ci++) {
    pswf_yl[ci]=pswf_y[Npd/2-ci];
  }
 
  int Npad=(Nzd-Npd)/2; /* same for even/odd */
  int Nz0=Npd+2*Npad;
  /* now allocate zero padded memory for 2D PSWF, row major order */
  if ((pswfxy=(float*)calloc((size_t)Nz0*Nz0,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* pswfxy is zero, valid range for copying
    |000000000000..............000000000000|
      (Nz-Np)/2  ^             ^  (Nz-Np)/2
      =Npad           =Np          =Npad
                 0...........Np-1 
    zero pixel:
                      ^Npad+Np/2 : even Np ,left (Npad+Np/2) right (Npad+Np/2)
                      ^Npad+Np/2 : odd Np, left (Npad+Np/2) right (Npad+Np/2+1)
    after oversampling
          Np=(int)((r-1)*Np*0.5) * 2+Np
          so Npad<=Npad+(int)((r-1)*Np*0.5) for copying back
  */
  for (nrow=Npad; nrow<Npad+Npd; nrow++) {
    /* copy to right hand side Np2 pixels,omit 0 pixel */
    memcpy(&pswfxy[nrow*Nz0+Npad+Npd/2+1],&pswf_y[1],(size_t)Np2*sizeof(float));
    /* copy to left hand side Np/2+1 pixels */
    memcpy(&pswfxy[nrow*Nz0+Npad],pswf_yl,(size_t)(Npd/2+1)*sizeof(float));
    /* now scale this row */
    ci=nrow-Npad;
    if (ci<Nl) {
     scaleval=pswf_yl[ci];
    } else {
     ci=ci-Nl;
     scaleval=pswf_y[ci+1];
    }
    my_fscal(Npd,scaleval,&pswfxy[nrow*Nz0+Npad]);
  }
  free(pswf_x);
  free(pswf_y);
  free(pswf_yl);

  struct wcsprm *mwcs;
  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_wcs_t *threaddata;
  /* divide data rows over threads */
  int NrowsNt=(Npd+Nt-1)/Nt;
  int nth,nth1,Nthb;

  if ((mwcs=(struct wcsprm*)calloc((size_t)Nt,sizeof(struct wcsprm)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  for (ci=0; ci<Nt; ci++) {
    generate_def_wcs(&mwcs[ci],(double)delta_ld*C180_PI,(double)Npd*0.5,ra0*C180_PI,dec0*C180_PI,0);
  }
  /* setup threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

  if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_wcs_t*)malloc((size_t)Nt*sizeof(thread_data_wcs_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* generate l,m grid NpxNp 
    (sqrt(1-r^2) -1)
    r^2=l^2+m^2  
    Np*delta_l = image FOV = Nx*delta_l0, delta_l0: image pixel width 
    so delta_l=Nx*delta_l0/Np

    if snapshot mode,
    l<= l+a*(sqrt(1-r^2)-1)
    m<= m+b*(sqrt(1-r^2)-1)
 */
  double *lmgrid,*denom;
  if ((lmgrid=(double*)calloc((size_t)Npd*Npd,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((denom=(double*)calloc((size_t)Npd*Npd,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* iterate over threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Npd; nth++) {
    if (ci+NrowsNt<Npd) {
     Nthb=NrowsNt;
    } else {
     Nthb=Npd-ci;
    }
    threaddata[nth].lmgrid=lmgrid;
    threaddata[nth].denom=denom;
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].Np=Npd;
    threaddata[nth].mwcs=&mwcs[nth];
    threaddata[nth].snapshot=snapshot;
    threaddata[nth].wpa=wpa;
    threaddata[nth].wpb=wpb;
    pthread_create(&th_array[nth],&attr,wcsgrid_threadfn,(void*)(&threaddata[nth]));
    /* next row set */
    ci=ci+NrowsNt;
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  free(threaddata);

  for (ci=0; ci<Nt; ci++) {
   wcsfree(&mwcs[ci]);
  }
  free(mwcs);

  /* scale by -2pi */
  my_dscal(Npd*Npd,-2.0*M_PI,lmgrid);


  complex float *din;
  if ((din=(complex float*)calloc((size_t)Nz0*Nz0,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  fftwf_plan plan0;
  do_create_fftw_plan_2d(&plan0, din, din, Nz0, Nz0, Nt);

  thread_data_lmpswf_t *threaddataN;
  if ((threaddataN=(thread_data_lmpswf_t*)malloc((size_t)Nt*sizeof(thread_data_lmpswf_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }

/* DEBUG WPROJ */
//  FILE *bfilep=fopen("debug_w_in_r","w+");
//  fprintf(bfilep,"%d %d\n",Npd,Nw);

  /* form the product pswf[] x exp{j2 pi w (...) } */
  /* initialize threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Npd; nth++) {
    if (ci+NrowsNt<Npd) {
     Nthb=NrowsNt;
    } else {
     Nthb=Npd-ci;
    }
    threaddataN[nth].lmgrid=lmgrid;
    threaddataN[nth].denom=denom;
    threaddataN[nth].pswfxy=pswfxy;
    threaddataN[nth].din=din;
    threaddataN[nth].startrow=ci+Npad;
    threaddataN[nth].Nrows=Nthb;
    threaddataN[nth].Np=Npd;
    threaddataN[nth].Npad=Npad;
    threaddataN[nth].Nz0=Nz0;
    /* next row set */
    ci=ci+NrowsNt;
  }

  for (nw=0; nw<Nw;nw++) {
   for(nth1=0; nth1<nth; nth1++) {
    threaddataN[nth1].w=wparr[nw];
    pthread_create(&th_array[nth1],&attr,lmpswf_threadfn,(void*)(&threaddataN[nth1]));
   }
   /* now wait for threads to finish */
   for(nth1=0; nth1<nth; nth1++) {
    pthread_join(th_array[nth1],NULL);
   }
/* DEBUG WPROJ */
//  for (nrow=Npad; nrow<Npad+Npd; nrow++) {
//   for (ci=Npad; ci<Npad+Npd; ci++) {
//    fprintf(bfilep,"%f %f\n",crealf(din[nrow*Nz0+ci]),cimagf(din[nrow*Nz0+ci]));
//   }
//  }

    /* do fftshift */
   do_fftshift_inplace(din,Nz0,Npad+Npd/2);

   /* do fft */
   fftwf_execute(plan0);
   do_ifftshift_inplace(din,Nz0,Npad+Npd/2);
   
  /* ignore padded values */
   /* copy back to the kernel array */
   for (nrow=Npad+Npad1; nrow<Npad+Npad1+Np; nrow++) {
    memcpy(&wkernel[nw*Np*Np+(nrow-Npad-Npad1)*Np],&din[nrow*Nz0+Npad+Npad1],(size_t)Np*sizeof(complex float));
   }
   /* reset din to zero */
   memset(din,0,(size_t)Nz0*Nz0*sizeof(complex float));

  }

  /* also calculate kernel for w=0, regardless if it is given as input,
     this is needed to scale all kernels */
  for(nth1=0; nth1<nth; nth1++) {
    threaddataN[nth1].w=0.0f;
    pthread_create(&th_array[nth1],&attr,lmpswf_threadfn,(void*)(&threaddataN[nth1]));
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
    pthread_join(th_array[nth1],NULL);
  }
  /* do fftshift */
  do_fftshift_inplace(din,Nz0,Npad+Npd/2);
  /* do fft */
  fftwf_execute(plan0);
  do_ifftshift_inplace(din,Nz0,Npad+Npd/2);
 
  /* find peak value 
   scale such that zero pixel of 1st w plane is 1 */
  float abs_re;
  abs_re=cabsf(din[(Nz0/2)*Nz0+Nz0/2]); 
//printf("peak value=%f\n",abs_re);
  my_csscal(Np*Np*Nw,1.0f/abs_re,wkernel);

  
/* DEBUG WPROJ */
//  fclose(bfilep);

  free(threaddataN);
  pthread_attr_destroy(&attr);
  free(th_array);

  do_destroy_fftw_plan(plan0);
  free(pswfxy);
  free(din);
  free(lmgrid);
  free(denom);

  float tol0;
  *maxsupport=0;
  /* scan the zero pixel row to find support in X direction */
  /* scan the zero pixel column to find support in Y direction */
  /* support determined by edge/center > tol */
  for (nw=0; nw<Nw; nw++) {
    supportX[nw]=0;
    tol0=tol*cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2]);
    for (ci=Np/2; ci<Np; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+ci]);
      if (tol0>abs_re) {
        break;
      } else {
        supportX[nw]++;
      }
    }

    supportY[nw]=0;
    for (ci=0; ci<Np/2; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2+ci*Np]);
      if (tol0>abs_re) {
        break;
      } else {
        supportY[nw]++;
      }
    }


    /* scale back support to original pixels (without padding) */
    abs_re=(float)(supportX[nw]*Np)/(float)Nz;
    supportX[nw]=(int)(abs_re*2.0f)+1;
    if (supportX[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->X  plane %d is %d, exceeds %d\n",nw,supportX[nw],Np);
     supportX[nw]=Np;
    }
    abs_re=(float)(supportY[nw]*Np)/(float)Nz;
    supportY[nw]=(int)(abs_re*2.0f)+1;
    if (supportY[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->Y plane %d is %d, exceeds %d\n",nw,supportY[nw],Np);
     supportY[nw]=Np;
    }
    /* also find max support */
    if (supportX[nw]>*maxsupport) {
      *maxsupport=supportX[nw];
    }
  }
  printf("Maximum convolution kernel support is %d pixels but allowed is %d.\n",*maxsupport*Nz/Np,Np);

  printf("[");
  for(nw=0; nw<Nw; nw++) { 
   printf("(%d %d %d %f) ",nw,supportX[nw],supportY[nw],wparr[nw]);
  } 
  printf("]\n");


  return 0;
}


#ifdef HAVE_CUDA
/* oversample when generating w kernels, but still select Np pixels
   r: oversample ratio >1.0
   both zero padding and lm-pswf evaluation is oversampled by this ratio (Np x r)
   but after FFT, Np pixels again copied back
  GPU version */
int
generate_2d_wconv_kernels_over_gpu(int M, int Np, float delta_l, float *wparr, int Nw, int Nz, complex float *wkernel, float tol, int *supportX, int *supportY, int *maxsupport, double ra0, double dec0, double wpa, double wpb, int snapshot,  double r, int Nt) {
  /* zero pixel
    <---------------> Np/2+1 values
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero  for even Np : each half Np/2 pixels
    |0,1,....,Np/2-1,Np/2,Np/2+1,.....,Np-1|
                     ^zero for odd Np : left Np/2 right Np/2+1 pixels
  */
  int ci,nw;
  /* arrays to calculate PSWF values in [0,1] 
     [0,1] is sampled by Np values, the values beyond that upto Nz
     are additions
  */
  float *pswf_x,*pswf_y,*pswf_yl,*pswfxy;
  float pswf_c=M_PI*0.5f*(float)M;
  int Np2,nrow;
  /* oversample Np */
  int Npad1=(int)((r-1.0)*(double)Np*0.5);
  int Npd=Np+2*Npad1;
  /* oversample Nz */
  int Nzd=(int)(r*(double)Nz);
  /* undersample delta_l */
  float delta_ld=delta_l/(float)r;

  /* get right half length for even/odd Np */
  if(Npd%2) { /* odd */
   Np2=Npd/2;
  } else { /* even : note right side is shorter as 0 pixel is taken from left */
   Np2=Npd/2-1;
  }

  int Nl=Npd/2+1;
  if ((pswf_x=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_y=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_yl=(float*)calloc((size_t)Nl,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }


  float scaleval;
  float pswf_delta=1.0f/((float)(Nl-1)); /* 1,2,..Nl : Nl-1 intervals */
  pswf_x[0]=0.0f;//+0.5*pswf_delta; FIXME why 0.5pixel?
  for (ci=1; ci<Nl; ci++) {
      pswf_x[ci]=pswf_x[ci-1]+pswf_delta;
  }
  pswf_nn(pswf_c,Nl,pswf_x,pswf_y);

  /* copy right hand side to left hand side*/
  for (ci=0; ci<Nl; ci++) {
    pswf_yl[ci]=pswf_y[Npd/2-ci];
  }
 
  int Npad=(Nzd-Npd)/2; /* same for even/odd */
  int Nz0=Npd+2*Npad;
  /* now allocate zero padded memory for 2D PSWF, row major order */
  if ((pswfxy=(float*)calloc((size_t)Nz0*Nz0,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* pswfxy is zero, valid range for copying
    |000000000000..............000000000000|
      (Nz-Np)/2  ^             ^  (Nz-Np)/2
      =Npad           =Np          =Npad
                 0...........Np-1 
    zero pixel:
                      ^Npad+Np/2 : even Np ,left (Npad+Np/2) right (Npad+Np/2)
                      ^Npad+Np/2 : odd Np, left (Npad+Np/2) right (Npad+Np/2+1)
    after oversampling
          Np=(int)((r-1)*Np*0.5) * 2+Np
          so Npad<=Npad+(int)((r-1)*Np*0.5) for copying back
  */
  for (nrow=Npad; nrow<Npad+Npd; nrow++) {
    /* copy to right hand side Np2 pixels,omit 0 pixel */
    memcpy(&pswfxy[nrow*Nz0+Npad+Npd/2+1],&pswf_y[1],(size_t)Np2*sizeof(float));
    /* copy to left hand side Np/2+1 pixels */
    memcpy(&pswfxy[nrow*Nz0+Npad],pswf_yl,(size_t)(Npd/2+1)*sizeof(float));
    /* now scale this row */
    ci=nrow-Npad;
    if (ci<Nl) {
     scaleval=pswf_yl[ci];
    } else {
     ci=ci-Nl;
     scaleval=pswf_y[ci+1];
    }
    my_fscal(Npd,scaleval,&pswfxy[nrow*Nz0+Npad]);
  }
  free(pswf_x);
  free(pswf_y);
  free(pswf_yl);

  struct wcsprm *mwcs;
  pthread_attr_t attr;
  pthread_t *th_array;
  thread_data_wcs_t *threaddata;
  /* divide data rows over threads */
  int NrowsNt=(Npd+Nt-1)/Nt;
  int nth,nth1,Nthb;

  if ((mwcs=(struct wcsprm*)calloc((size_t)Nt,sizeof(struct wcsprm)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  for (ci=0; ci<Nt; ci++) {
    generate_def_wcs(&mwcs[ci],(double)delta_ld*C180_PI,(double)Npd*0.5,ra0*C180_PI,dec0*C180_PI,0);
  }
  /* setup threads */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

  if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  if ((threaddata=(thread_data_wcs_t*)malloc((size_t)Nt*sizeof(thread_data_wcs_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* generate l,m grid NpxNp 
    (sqrt(1-r^2) -1)
    r^2=l^2+m^2  
    Np*delta_l = image FOV = Nx*delta_l0, delta_l0: image pixel width 
    so delta_l=Nx*delta_l0/Np

    if snapshot mode,
    l<= l+a*(sqrt(1-r^2)-1)
    m<= m+b*(sqrt(1-r^2)-1)
 */
  double *lmgrid,*denom;
  if ((lmgrid=(double*)calloc((size_t)Npd*Npd,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((denom=(double*)calloc((size_t)Npd*Npd,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* iterate over threads */
  ci=0;
  for (nth=0;  nth<Nt && ci<Npd; nth++) {
    if (ci+NrowsNt<Npd) {
     Nthb=NrowsNt;
    } else {
     Nthb=Npd-ci;
    }
    threaddata[nth].lmgrid=lmgrid;
    threaddata[nth].denom=denom;
    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].Np=Npd;
    threaddata[nth].mwcs=&mwcs[nth];
    threaddata[nth].snapshot=snapshot;
    threaddata[nth].wpa=wpa;
    threaddata[nth].wpb=wpb;
    pthread_create(&th_array[nth],&attr,wcsgrid_threadfn,(void*)(&threaddata[nth]));
    /* next row set */
    ci=ci+NrowsNt;
  }
  /* now wait for threads to finish */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }

  free(threaddata);

  for (ci=0; ci<Nt; ci++) {
   wcsfree(&mwcs[ci]);
  }
  free(mwcs);

  pthread_attr_destroy(&attr);
  free(th_array);


/************************************************************************/
  /* copy arrays to GPU  pswfxy (Nz0xNz0)
   lmgrid,denom (NpdxNpd) */
  
  float abs_re;
  evaluate_wplane_fft(Nz0,Np,Npd,Npad,Npad1,pswfxy,lmgrid,denom,(float*)wkernel, wparr, Nw, &abs_re);
//printf("peak value=%f\n",abs_re);
  my_csscal(Np*Np*Nw,1.0f/abs_re,wkernel);


  free(pswfxy);
  free(lmgrid);
  free(denom);

  float tol0;
  *maxsupport=0;
  /* scan the zero pixel row to find support in X direction */
  /* scan the zero pixel column to find support in Y direction */
  /* support determined by edge/center > tol */
  for (nw=0; nw<Nw; nw++) {
    supportX[nw]=0;
    tol0=tol*cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2]);
    for (ci=Np/2; ci<Np; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+ci]);
      if (tol0>abs_re) {
        break;
      } else {
        supportX[nw]++;
      }
    }

    supportY[nw]=0;
    for (ci=0; ci<Np/2; ci++) {
      abs_re=cabsf(wkernel[nw*Np*Np+(Np/2)*Np+Np/2+ci*Np]);
      if (tol0>abs_re) {
        break;
      } else {
        supportY[nw]++;
      }
    }


    /* scale back support to original pixels (without padding) */
    abs_re=(float)(supportX[nw]*Np)/(float)Nz;
    supportX[nw]=(int)(abs_re*2.0f)+1;
    if (supportX[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->X  plane %d is %d, exceeds %d\n",nw,supportX[nw],Np);
     supportX[nw]=Np;
    }
    abs_re=(float)(supportY[nw]*Np)/(float)Nz;
    supportY[nw]=(int)(abs_re*2.0f)+1;
    if (supportY[nw]>Np) {
     fprintf(stderr,"Warning: support of W ->Y plane %d is %d, exceeds %d\n",nw,supportY[nw],Np);
     supportY[nw]=Np;
    }
    /* also find max support */
    if (supportX[nw]>*maxsupport) {
      *maxsupport=supportX[nw];
    }
  }
  printf("Maximum convolution kernel support is %d pixels but allowed is %d.\n",*maxsupport*Nz/Np,Np);

  printf("[");
  for(nw=0; nw<Nw; nw++) { 
   printf("(%d %d %d %f) ",nw,supportX[nw],supportY[nw],wparr[nw]);
  } 
  printf("]\n");


  return 0;
}
#endif /* HAVE_CUDA */
