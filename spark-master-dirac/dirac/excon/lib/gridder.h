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


#ifndef GRIDDER_H
#define GRIDDER_H
#ifdef __cplusplus
        extern "C" {
#endif

#ifdef HAVE_CUDA
#include "cuda.h"
#include <cuComplex.h>
#endif
#include <stdio.h>
#include <pthread.h>
#include <glib.h>

/* max possible GPUs */
#ifndef MAX_GPU
#define MAX_GPU 4 
#endif


/* pipeline state values */
#ifndef PT_DO_NOTHING
#define PT_DO_NOTHING 0
#endif
#ifndef PT_DO_AGPU
#define PT_DO_AGPU 1 /* allocate GPU memory, attach GPU */
#endif
#ifndef PT_DO_DGPU
#define PT_DO_DGPU 2 /* free GPU memory, detach GPU */
#endif
#ifndef PT_DO_WORK_GRID /* grid work */
#define PT_DO_WORK_GRID 3
#endif

/* data size (how many 4 bytes) */
/* buffer length for slave threads : use ~1 MB*/
#ifndef DATA_BUF_LEN
#define DATA_BUF_LEN 262144
#endif

/* different weighting modes */
#ifndef WMODE_UNIFORM
#define WMODE_UNIFORM 0
#endif
#ifndef WMODE_ROBUST
#define WMODE_ROBUST 1
#endif
#ifndef WMODE_NATURAL
#define WMODE_NATURAL 2
#endif
#ifndef WMODE_PIPEMENON
#define WMODE_PIPEMENON 3
#endif
#ifndef WMODE_PMTAPER
#define WMODE_PMTAPER 4 /* uniform weight+taper gridded data */
#endif


/* which global weight function to use in adaptive weighting */
#ifndef CONV_MODE_UNIFORM
#define CONV_MODE_UNIFORM 0
#endif
#ifndef CONV_MODE_NCP
#define CONV_MODE_NCP 1 /* NCP 50-250lambda imaging (use 30-300 limits) */
#endif
#ifndef CONV_MODE_NCP_HD
#define CONV_MODE_NCP_HD 2 /* NCP highres imaging 0-40klambda */
#endif
#ifndef CONV_MODE_50_250
#define CONV_MODE_50_250 3 /* 50-250 lambda flat, with edges planck taper */
#endif


/******************** barrier.c **********************/
typedef struct t_barrier_ {
  int tcount; /* current no. of threads inside barrier */
  int nthreads; /* the no. of threads the barrier works
                with. This is a constant */
  pthread_mutex_t enter_mutex;
  pthread_mutex_t exit_mutex;
  pthread_cond_t lastthread_cond;
  pthread_cond_t exit_cond;
} th_barrier;


/* initialize barrier */
/* N - no. of accomodated threads */
extern void
init_th_barrier(th_barrier *barrier, int N);

/* destroy barrier */
extern void
destroy_th_barrier(th_barrier *barrier);

/* the main operation of the barrier */
extern void
sync_barrier(th_barrier *barrier);
/******************** gridder.cu **********************/
#ifdef HAVE_CUDA
/* data used for gridding GPU master threads */
typedef struct visdata_ {
  int N; /* how many points: how many slave threads in GPU */ 
  float *u,*v,*w; /* Nx1 */
  float *wt; /* data weights N x 1 */
  float2 *xx,*yy,*xy,*yx; /* Nx1 */
  int *flag;
} visdata;

/* data used for gridding GPU slave threads */
typedef struct visdata_slave_ {
  float u,v,w,wt;
  int flag;
  float2 xx,yy,xy,yx;
} visdata_slave;


/* data struct shared by all threads */
typedef struct gb_gdata_ {
  int *status; /* 0: do nothing, 
              1: allocate GPU  memory, attach GPU
              2: free GPU memory, detach GPU 
              3,4..: do work on GPU 
              99: reset GPU memory (memest all memory) */

  visdata *vis;
  float lambda;
  float uvscale;
  float deltaU;
  float expW;
  /* output pointers */
  float *uvgrid; /* complex */
  float *psfgrid; /* complex */
  /* mutex for writing output */
  pthread_mutex_t *writelock_img; 
  pthread_mutex_t *writelock_psf; 

  int Nx,Ny;
  float *wparr;
  int Nw;
  int *wpsupportX;
  int *wpsupportY;
  int maxsupport;
  float *wkernel; /* complex */
  int Np;
  int Nz;
  int imgmode;
  /* GPU related info */
  void *carray[MAX_GPU]; /* keep 4 here because only 2 to 4 will be used */
} gbgdata;

/* structs for thread pool (reusable), using a barrier */
/* slave thread data struct */
typedef struct slave_tdata_ {
  struct pipeline_ *pline; /* forward declaration */
  int tid; /* 0,1 used as the GPU card */
} slave_tdata;

/* pipeline struct */
typedef struct pipeline_ {
  void *data; /* all data needed by N threads */
  int terminate; /* 1: terminate, default 0*/
  int N; /* no. of slaves */
  pthread_t *slave;
  slave_tdata **sd; /* note recursive types */
  th_barrier gate1;
  th_barrier gate2;
  pthread_attr_t attr;
} th_pipeline;


extern void
init_pipeline(th_pipeline *pline, int Ngpu, void *data);
extern void
destroy_pipeline(th_pipeline *pline);


/* data struct shared by all threads */
typedef struct gb_slave_gdata_ {
  int status; /* 0: do nothing, 
              1: allocate GPU  memory, attach GPU
              2: free GPU memory, detach GPU 
              3,4..: do work on GPU 
              99: reset GPU memory (memest all memory) */
  int card;
  visdata_slave vis;
  float lambda;
  float uvscale;
  float deltaU;
  /* output pointers */
  float *uvgrid; /* complex */
  float *psfgrid; /* complex */
 /* mutex for writing output */
  pthread_mutex_t *writelock_img;
  pthread_mutex_t *writelock_psf;

  int Nx,Ny;
  float maxW; /* max value of W planes */
  float expW; /* 0.5 for sqrt() spacing, 1 for linear spacing w^(expW) space */
  float *wparr; /* note: this is a GPU pointer */
  int Nw;
  int *wpsupportX;
  int *wpsupportY;
  int maxsupport;
  int Np;
  int Nz;
  int imgmode;
} gb_slave_gdata;

/* structs for thread pool (reusable), using a barrier */
/* slave thread data struct */
typedef struct slave_slave_tdata_ {
  struct slave_pipeline_ *pline; /* forward declaration */
  int tid; /* 0,1,... */
} slave_slave_tdata;

/* pipeline struct */
typedef struct slave_pipeline_ {
  void *data; /* all data needed by N threads */
  int terminate; /* 1: terminate, default 0*/
  pthread_t *slave;
  slave_slave_tdata **sd; /* note recursive types */
  th_barrier gate1;
  th_barrier gate2;
  pthread_attr_t attr;
} th_slave_pipeline;

extern void
init_slave_pipeline(th_slave_pipeline *pline, int N, void *data);
extern void
destroy_slave_pipeline(th_slave_pipeline *pline,int N);
#endif
/************************** wkcuda.cu ***********************************/

#ifdef HAVE_CUDA
extern void
evaluate_wplane_fft(int Nz0, int Np, int Npd, int Npad, int Npad1, float *pswfxy, double *lmgrid, double *denom, float *wkernel, float *wparr, int Nw, float *peakval);
#endif

/******************** pipe_menon.cu **********************/
typedef struct fcomp_ {
 float x,y;
} fcomp;

/* struct for one row of data array data::writedata */
typedef struct iodata_ {
  fcomp xx,xy,yx,yy;
  float u,v,w;
  float wt;
  int flag;
} iodata;

/* list of points in a bucket */
typedef struct uvlist_ {
 GList *pix; /* list of points: */
 /* following used after the list is built */
 int P; /* total no of points */
 unsigned long int *id; /* Px1 arrays: row id of each point */
 float *u,*v,*w; /* Px1 arrays: scaled u,v coord of each point */

} uvlist;

#ifdef HAVE_CUDA
/* data struct shared by all threads */
typedef struct gb_pmgdata_ {
  int *status; /* 0: do nothing, 
              1: allocate GPU  memory, attach GPU
              2: free GPU memory, detach GPU 
              3,4..: do work on GPU 
              99: reset GPU memory (memest all memory) */
  /* data array */
  iodata *darr;
  unsigned long int *startrow,*Nrows;


  float uvscale;
  int convmode;
  int Nx,Ny;
  int B;
  pthread_mutex_t *writelock_hash;
  GHashTable *ht;

  /* output pointers */
  float *wold,*wnew; /* weight values */

  float *wkernel; /* zero w plane */
  int Np;
  int M;
  /* GPU related info */
  void **carray;
} gbpmgdata;

/* structs for thread pool (reusable), using a barrier */
/* slave thread data struct */
typedef struct slave_pm_tdata_ {
  struct pipeline_pm_ *pline; /* forward declaration */
  int tid; /* 0,1,.. */
} slave_pmtdata;

/* pipeline struct */
typedef struct pipeline_pm_ {
  void *data; /* all data needed by two threads */
  int terminate; /* 1: terminate, default 0*/
  pthread_t *slave;
  slave_pmtdata *sd; /* note recursive types */
  th_barrier gate1;
  th_barrier gate2;
  pthread_attr_t attr;
} th_pmpipeline;


extern void
init_pm_pipeline(th_pmpipeline *pline, void *data, int N);

extern void
destroy_pm_pipeline(th_pmpipeline *pline, int N);
#endif /* HAVE_CUDA */


/******************** ptree.c **********************/
/* Tree for spatially partitioning 2D point set into blocks 
  of equal points 
  But we do not use actual tree stucture (only the leaves), 
  so a sorted list is fine
  But But to preserve locality, traversing the list in the order of 
  nodes in the tree is better.
*/

/* node */
typedef struct ptreenode__ {
   int id; /* unique node no */
   int level; /* even/odd : split in x or y direction */
   int N; /* how many points are included in this node */
   int *index; /* array of row numbers of points included */
   double x_min,x_max,y_min,y_max; /* bounding rectangle of points */

} ptreenode;


/* tree, with auxiliary data */
typedef struct ptree__ {
  GList *leaves; /* list of leaf nodes, ordered insert with N */
  int count; /* how many leaf nodes, increment with each insert */
  int *unflagged; /* original indices of unflagged data, Dx1 array */
  int D; 
} ptree;

extern int
init_ptree(ptree *p, iodata *data, int N);

extern int
destroy_ptree(ptree *p);

extern int
split_ptree(ptree *p, iodata *data, int N);

extern int
sort_data_with_ptree(ptree *p, iodata *data, int N);

extern int
sort_data(iodata *data, int N, int M, int Nthreads);



/******************** flagger.c **********************/
/* flag data: first calculate standard deviation of XX,XY,YX,YY
              and flag data points that are XXs,XYs,YXs,YYs sigmas away
   iodata: Nx1 data points */
extern int
flag_data_sigma(iodata *data,int N, float XXs, float XYs, float YXs, float YYs, int Nthreads);
#ifdef __cplusplus
     } /* extern "C" */
#endif
#endif /* GRIDDER_H */
