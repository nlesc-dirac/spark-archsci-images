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



#include <complex.h>
#include <fftw3.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "largefft.h"

/* create,destroy plans */
void
init_fftw(int Nt) {
  /* threads */
  fftwf_init_threads();
  fftwf_plan_with_nthreads(Nt);
}

void
stop_fftw(void) {
  fftwf_cleanup_threads();
}

/* create,destroy plans */
int
do_create_fftw_plan(fftwf_plan *p, complex float *d, complex float *im, int n0, int n1) {
  *p=fftwf_plan_dft_2d(n0, n1, d, im,
          FFTW_BACKWARD,FFTW_ESTIMATE|FFTW_UNALIGNED);
  return 0;
}

int
do_destroy_fftw_plan(fftwf_plan p) {
  fftwf_destroy_plan(p);
  return 0;
}


/* do FFTshift of in to out, both size Nx1 arrays 
   N even 
     [N/2,N-1] -> [0,N/2-1] (N/2 values) and [0,N/2-1] -> [N/2,N-1]  (N/2 values)
   N odd
     [N/2+1,N-1] -> [0,N/2-1] (N/2 values) and [0,N/2] -> [N/2,N-1]  (N/2+1 values)
*/
static void
do_fftshift_1d(float *in, float *out, int N) {
  int n2=N/2;
  if (N%2) { /* odd */
    memcpy((void*)&out[0],(void*)&in[n2+1],(size_t)(n2)*sizeof(float));
    memcpy((void*)&out[n2],(void*)&in[0],(size_t)(n2+1)*sizeof(float));
  } else { /* even */
    memcpy((void*)&out[0],(void*)&in[n2],(size_t)(n2)*sizeof(float));
    memcpy((void*)&out[n2],(void*)&in[0],(size_t)(n2)*sizeof(float));
  }
}

/* IFFTshift to undo the FFTshift */
/* N even
    [0,N/2-1] -> [N/2,N-1]  (N/2 values) and [N/2,N-1] -> [0,N/2-1] (N/2 values)
   N odd
    [0,N/2-1] -> [N/2+1,N-1]  (N/2 values) and [N/2,N-1] -> [0,N/2]  (N/2+1 values)
*/
int
ifftshift_index(int N, int x) {
 if (N%2) { /* odd */
  return((x+N/2)%N);
 } else { /* even */
  return((x+N/2+1)%N);
 }
 return 0;
}


/* do FFT shift of input complex array, (row major order)
   also scale data points by scalefac
   and write back in column major order,
   real and imaginary parts seperately
   also performa apodization correction  
 
   din: input data array, size n0xn1 complex float
   dout: output data array, size 2*n0xn1 float
   din: row major, dout: col major 
   scalefac: scale for correction of FFT scaling
   M:  PSWF bandwidth in uv pixels, used for calculation of apodization correction
   Nx0: unpadded image size
   Npad: padding, so n0=n1=Nx0+2*Npad

   tol: cutoff tolerance for apodization correction

   output, the first Nx0*Nx0 values starting at dout[0] is the true image (real)
   and Nx0*Nx0 values starting at dout[2*n0*n1] is  the imaginary image
 */
int
do_data_fftshift(complex float *din, float *dout, float scalefac, int n0, int n1,int M, int Nx0, int Npad, float tol) {
 /* we to the shift in several steps */
  /* true image => |A B|    din => |D C|
                   |C D|           |B A|
    1) din: read rows (size n1 x 1) and perform shift, write back to din
       as the same row
        din => |C D|
               |A B|
    2) read columns of din (size n0 x 1) step size n1, perform shift
       write back to dout, as a row, each row of size n0 x 1 
        dout => |A C|
                |B D|
    3) now when dout is read in column major order, true image recovered */


  float *img;
  float *xtmp,*ytmp;
  int nmax,ci;
  if (n0>n1) { nmax=n0; } else { nmax=n1; }
  /* allocate memory */
  if ((xtmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ytmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* PSWF calculation */
  /* PSWF 1....Np|0.....0| zero padding Np*Nzc
     so total = 1....Np|0.....0|0.....0|Np-1....1 = 2*Np+2*Nzc
     do FFT,
     truncate again 1...Np|0....... 0|Np-1....1
     do IFFT
     extract 1....Np
  */
  float pswf_c=M_PI*0.5f*(float)M;
  int N0=nmax/2; /* 0 pixel */
  int Np=nmax/2+1;
  float pswf_delta=1.0f/((float)(nmax/2));
  float *pswf_x,*pswf_y;
  if ((pswf_x=(float*)calloc((size_t)Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((pswf_y=(float*)calloc((size_t)Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  pswf_x[0]=0.0f;
  for (ci=1; ci<Np; ci++) {
      pswf_x[ci]=pswf_x[ci-1]+pswf_delta;
  }
  //float frat=(float)(M_PI*2.0)/(float)(M)*(float)nmax/(float)(nmax-2*Npad);
  float tmpf;
  pswf_nn(pswf_c,Np,pswf_x,pswf_y);
  for (ci=1; ci<Np; ci++) { /* assert 0-th value is 1.0 */
    if (fabs(pswf_y[ci])>tol) {
     //pswf_y[ci]=1.0f/(pswf_y[ci]);
     tmpf=pswf_x[ci]; //*frat;
     ///* sinc() function is because of the pixel width */
     pswf_y[ci]=tmpf/(sinf(tmpf)*pswf_y[ci]);
    } else {
     pswf_y[ci]=1.0f;
    }
  }
  float ipswf; /* scale factor for apodization correction */
/*
  for (ci=0; ci<Np; ci++) {
   printf("%f %f\n",pswf_x[ci],pswf_y[ci]);
  }
*/
  free(pswf_x);

  for (ci=0; ci<n0; ci++) {
   /* read rows of size n1 x 1, stride is 2 because its complex */
   img=(float*)&din[ci*n1];
   my_fcopy(n1,img,2,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n1);
   /* scale data : with FFT scale and inv PSWF scale */
   if (ci<N0) {
    ipswf=pswf_y[ci+1];
   } else {
    ipswf=pswf_y[n0-ci-1];
   }
   my_fscal(n1,scalefac*ipswf,xtmp);
   my_fcopy(n1,xtmp,1,img,2);

   /* now the imaginary part */
   my_fcopy(n1,&img[1],2,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n1);
   my_fscal(n1,scalefac*ipswf,xtmp);
   my_fcopy(n1,xtmp,1,&img[1],2);
  }

  /* note: we remove padded values */  
  unsigned long int writeoff=0; /* remember next offset to write */
  /* only start from columns Npad */
  for (ci=Npad; ci<Npad+Nx0; ci++) {
   /* read columns of size n0 x 1 */
   img=(float*)&din[ci];
   my_fcopy(n0,img,2*n1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n0);
   /* scale data : with FFT scale and inv PSWF scale */
   if (ci<N0) {
    ipswf=pswf_y[N0-ci];
   } else {
    ipswf=pswf_y[ci-N0];
   }
   my_fscal(n0,ipswf,xtmp);
   /* only write values [Npad,Nx0+Npad-1] */
   my_fcopy(Nx0,&xtmp[Npad],1,&dout[writeoff],1);

   /* now the imaginary part */
   my_fcopy(n0,&img[1],2*n1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n0);
   my_fscal(n0,ipswf,xtmp);
   my_fcopy(Nx0,&xtmp[Npad],1,&dout[n1*n0+writeoff],1);

   writeoff+=Nx0;
  } 

  free(xtmp);
  free(ytmp);
  free(pswf_y);
  return 0;
}


int
do_weight_fftshift(float *din, float *dout, int n0, int n1, int Nx0, int Npad) {
 /* we to the shift in several steps */
  /* true image => |A B|    din => |D C|
                   |C D|           |B A|
    1) din: read rows (size n1 x 1) and perform shift, write back to din
       as the same row
        din => |C D|
               |A B|
    2) read columns of din (size n0 x 1) step size n1, perform shift
       write back to dout, as a row, each row of size n0 x 1 
        dout => |A C|
                |B D|
    3) now when dout is read in column major order, true image recovered */


  float *img;
  float *xtmp,*ytmp;
  int nmax,ci;
  if (n0>n1) { nmax=n0; } else { nmax=n1; }
  /* allocate memory */
  if ((xtmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ytmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  for (ci=0; ci<n0; ci++) {
   /* read rows of size n1 x 1 */
   img=(float*)&din[ci*n1];
   my_fcopy(n1,img,1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n1);
   my_fcopy(n1,xtmp,1,img,1);
  }

  /* note: we remove padded values */  
  unsigned long int writeoff=0; /* remember next offset to write */
  /* only start from columns Npad */
  for (ci=Npad; ci<Npad+Nx0; ci++) {
   /* read columns of size n0 x 1 */
   img=(float*)&din[ci];
   my_fcopy(n0,img,n1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n0);
   /* only write values [Npad,Nx0+Npad-1] */
   my_fcopy(Nx0,&xtmp[Npad],1,&dout[writeoff],1);

   writeoff+=Nx0;
  } 

  free(xtmp);
  free(ytmp);
  return 0;
}


int
do_grid_fftshift(complex float *din, float *dout, float scalefac, float *scratch, int n0, int n1, int Nx0, int Npad) {
 /* we to the shift in several steps */
  /* true image => |A B|    din => |D C|
                   |C D|           |B A|
    1) din: read rows (size n1 x 1) and perform shift, write back to din
       as the same row
        din => |C D|
               |A B|
    2) read columns of din (size n0 x 1) step size n1, perform shift
       write back to dout, as a row, each row of size n0 x 1 
        dout => |A C|
                |B D|
    3) now when dout is read in column major order, true image recovered */


  float *img,*outimg;
  float *xtmp,*ytmp;
  int nmax,ci;
  if (n0>n1) { nmax=n0; } else { nmax=n1; }
  /* allocate memory */
  if ((xtmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ytmp=(float*)calloc((size_t)nmax,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  /* real part */
  for (ci=0; ci<n0; ci++) {
   /* read rows of size n1 x 1, stride is 2 because its complex */
   img=(float*)&din[ci*n1];
   outimg=(float*)&scratch[ci*n1];
   my_fcopy(n1,img,2,ytmp,1);
   my_fscal(n1,scalefac,ytmp);
   do_fftshift_1d(ytmp,xtmp,n1);
   my_fcopy(n1,xtmp,1,outimg,1);
  }

  /* note: we remove padded values */  
  unsigned long int writeoff=0; /* remember next offset to write */
  /* only start from columns Npad */
  for (ci=Npad; ci<Npad+Nx0; ci++) {
   /* read columns of size n0 x 1 */
   img=(float*)&scratch[ci];
   my_fcopy(n0,img,n1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n0);
   /* only write values [Npad,Nx0+Npad-1] */
   my_fcopy(Nx0,&xtmp[Npad],1,&dout[writeoff],1);
   writeoff+=Nx0;
  } 
  /* imaginary part */
  for (ci=0; ci<n0; ci++) {
   /* read rows of size n1 x 1, stride is 2 because its complex */
   img=(float*)&din[ci*n1];
   outimg=(float*)&scratch[ci*n1];
   my_fcopy(n1,&img[1],2,ytmp,1);
   my_fscal(n1,scalefac,ytmp);
   do_fftshift_1d(ytmp,xtmp,n1);
   my_fcopy(n1,xtmp,1,outimg,1);
  }

  /* only start from columns Npad */
  for (ci=Npad; ci<Npad+Nx0; ci++) {
   /* read columns of size n0 x 1 */
   img=(float*)&scratch[ci];
   my_fcopy(n0,img,n1,ytmp,1);
   do_fftshift_1d(ytmp,xtmp,n0);
   /* only write values [Npad,Nx0+Npad-1] */
   my_fcopy(Nx0,&xtmp[Npad],1,&dout[writeoff],1);
   writeoff+=Nx0;
  } 
  free(xtmp);
  free(ytmp);
  return 0;
}



/* create,destroy plans */
int
do_create_fftw_plan_2d(fftwf_plan *p, complex float *din, complex float *dout, int n0, int n1, int Nt) {
  /* threads */
  fftwf_plan_with_nthreads(Nt);
  *p=fftwf_plan_dft_2d(n0, n1, din, dout,
          FFTW_BACKWARD,FFTW_ESTIMATE);
  return 0;
}

/* in-place fftshift */
/* in: NxN array, row major order
   N0: location of zero pixel
   [N0,N-1] -> [0,N-1-N0] and [0,N0-1] -> [N-N0,N-1]
*/
int
do_fftshift_inplace(complex float *in, int N, int N0) {
  complex float *tmp;
  int ci;
  if ((tmp=(complex float*)calloc((size_t)N,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* shift rows first */
  for (ci=0; ci<N; ci++) {
   /* 0  ot N-N0 */
   memcpy(&tmp[N-N0],&in[ci*N+0],(size_t)(N0)*sizeof(complex float));
   /* N0 to 0 */
   memcpy(&tmp[0],&in[ci*N+N0],(size_t)(N-N0)*sizeof(complex float));
   memcpy(&in[ci*N+0],tmp,(size_t)N*sizeof(complex float));
  } 
  /* shift columns now */
  for (ci=0; ci<N; ci++) {
   /* 0 to N-N0 */
   my_ccopy(N0,&in[ci],N,&tmp[N-N0],1);
   /* N0 to 0 */
   my_ccopy(N-N0,&in[N0*N+ci],N,&tmp[0],1);
   my_ccopy(N,tmp,1,&in[ci],N);
  } 
  free(tmp);
  return 0;
}

/* in-place ifftshift */
/* in: NxN array, row major order
   N0: location of zero pixel
   [0,N-1-N0] -> [N0,N-1]  and [N-N0,N-1] -> [0,N0-1]  
*/
int
do_ifftshift_inplace(complex float *in, int N, int N0) {
  complex float *tmp;
  int ci;
  if ((tmp=(complex float*)calloc((size_t)N,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* shift rows first */
  for (ci=0; ci<N; ci++) {
   /* 0 to N0 */
   memcpy(&tmp[N0],&in[ci*N+0],(size_t)(N-N0)*sizeof(complex float));
   /* N-N0 to 0 */
   memcpy(&tmp[0],&in[ci*N+N-N0],(size_t)(N0)*sizeof(complex float));
   memcpy(&in[ci*N+0],tmp,(size_t)N*sizeof(complex float));
  } 
  /* shift columns now */
  for (ci=0; ci<N; ci++) {
   /* 0 to N0 */
   my_ccopy(N-N0,&in[ci],N,&tmp[N0],1);
   /* N-N0 to 0 */
   my_ccopy(N0,&in[(N-N0)*N+ci],N,&tmp[0],1);
   my_ccopy(N,tmp,1,&in[ci],N);
  } 
  free(tmp);
  return 0;
}
