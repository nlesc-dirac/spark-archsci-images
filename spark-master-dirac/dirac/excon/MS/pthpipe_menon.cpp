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
#include <glib.h>

using namespace std;

typedef struct thread_data_gridhash_t_ {
  iodata *darr;
  unsigned long int startrow,Nrows;
  float lambda;
  float uvscale;
  float deltaU;
  int Nx,Ny; /* image pixel size */
  int B; /* bucket pixel size  BxB */
  pthread_mutex_t *writelock_hash;
  GHashTable *ht;
  float *wgrid; /* array storing gridded weights (binary file) */
  float *wold,*wnew; /* weight vectors used in the update, size equal to darr */

  /* for scale matched weights */
  float imscale;
  /* which global density function to use */
  int convmode;

  /* for evaluation of PSWF convolution */
  float *kernel;  /* NpxNp kernel */
  float *uvgrid; /* Np+2 valus of grid points (two guard values) */
  int Np;
  int M; /* conv. kernel support is [-M,M] pixels (buckets) */

} thread_data_gridhash_t;

/********************************/

/* data needed for each point in a list */
typedef struct gdpoint_ {
 unsigned long int id;
 float u,v,w; /* scaled uv */
} gdpoint;

/* hash table functions */
/* key destroy function */
static void
destroy_key(gpointer data) {
 free((gint64*)data);
}
/* value destroy function */
static void
destroy_value(gpointer data) {
  uvlist *uv=(uvlist*)data;
  if (uv->pix) {
   GList *li;
   for(li=uv->pix; li!=NULL; li=g_list_next(li)) {
     gdpoint *ii=(gdpoint*)li->data;
     g_free(ii);
   }
   g_list_free(uv->pix);
  }
  free(uv->id); 
  free(uv->u); 
  free(uv->v); 
  free(uv->w); 
  free(uv);
}
/********************************/

/* evaluate desired density function */
/* x = sqrt(u^2+v^2), x in [30,300] -> effective [50,250] */
static float
ncp_weight(float x) {
/*    fo(x) = 
              a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + 
              a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + 
              a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2) + 
              a7*exp(-((x-b7)/c7)^2)
    mean(fo(x)) is about 1
*/
 if (x<30.0f || x>300.0f) { return 0.0f; }
 if (x<82.51f) {  /* [30,82.51) */
  return 2.08f/(1.0f+expf(-(x-50.0f)*0.333333333f)); 
 }
 if (x>247.5f) {  /* (247.5,300] */
  float x2=(x-247.5f);
  return 1.115f*expf(-x2*x2*0.005f);
 }
 /* else [82.51,247.5] */
 float r[7];
 float a1 =-19.69f;
 float b1 =69.68f;
 float c1 =20.33f;
 float a2 =1.428f;
 float b2 =37.94f;
 float c2 =8.104f;
 float a3 =21.42f;
 float b3 =69.34f;
 float c3 =21.35f;
 float a4 =0.5998f;
 float b4 =255.7f;
 float c4 =35.21f;
 float a5 =1.231f;
 float b5 =194.4f;
 float c5 =57.97f;
 float a6 =0.64f;
 float b6 =284.1f;
 float c6 =18.88f;
 float a7 =1.115f;
 float b7 =110.5f;
 float c7 =36.72f;

 r[0]=(x-b1)/c1;
 r[1]=(x-b2)/c2;
 r[2]=(x-b3)/c3;
 r[3]=(x-b4)/c4;
 r[4]=(x-b5)/c5;
 r[5]=(x-b6)/c6;
 r[6]=(x-b7)/c7;
 r[0]*=-r[0];
 r[1]*=-r[1];
 r[2]*=-r[2];
 r[3]*=-r[3];
 r[4]*=-r[4];
 r[5]*=-r[5];
 r[6]*=-r[6];
 float sum=0.0f;
 sum+=a1*expf(r[0]);
 sum+=a2*expf(r[1]);
 sum+=a3*expf(r[2]);
 sum+=a4*expf(r[3]);
 sum+=a5*expf(r[4]);
 sum+=a6*expf(r[5]);
 sum+=a7*expf(r[6]);
 return sum;
}
/* x = sqrt(u^2+v^2), x in [40,285/800] */
static float
ncp_weight00(float x) {
/*    fo(x) = 
              a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + 
              a3*exp(-((x-b3)/c3)^2) + a4*exp(-((x-b4)/c4)^2) + 
              a5*exp(-((x-b5)/c5)^2) + a6*exp(-((x-b6)/c6)^2)
    mean(fo(x)) is about 1
*/
 if (x<40.0f) { return 0.0f; }
 if (x>800.0f) { 
  return 1.0f;
 }
 /* else [40,285] */
 float r[6];
 float a1 =-0.9415f;
 float b1 =117.1f;
 float c1 =15.08f;
 float a2 =5.231f;
 float b2 =49.57f;
 float c2 =13.79f;
 float a3 =2.209f;
 float b3 =67.29f;
 float c3 =14.86f;
 float a4 =10.43f;
 float b4 =72.19f;
 float c4 =200.8f;
 float a5 =104.9f;
 float b5 =98.72f;
 float c5 =65.8f;
 float a6 =-101.3f;
 float b6 =101.2f;
 float c6 =66.63f;

 r[0]=(x-b1)/c1;
 r[1]=(x-b2)/c2;
 r[2]=(x-b3)/c3;
 r[3]=(x-b4)/c4;
 r[4]=(x-b5)/c5;
 r[5]=(x-b6)/c6;
 r[0]*=-r[0];
 r[1]*=-r[1];
 r[2]*=-r[2];
 r[3]*=-r[3];
 r[4]*=-r[4];
 r[5]*=-r[5];
 float sum=0.0f;
 sum+=a1*expf(r[0]);
 sum+=a2*expf(r[1]);
 sum+=a3*expf(r[2]);
 sum+=a4*expf(r[3]);
 sum+=a5*expf(r[4]);
 sum+=a6*expf(r[5]);
 return (sum+1.0f); /* as x-> inf, goes to 1 */
}



/* evaluate desired density function */
/* x = sqrt(u^2+v^2), x in [50,40k]
 only zero baselines in[1400,2000]  */
static float
ncp_weight_highres(float x) {
 if (x<50.0f || x>40000.0f) { return 0.0f; }
 if (x<=2000.0f && x>=1400.0f) { return 0.0f; }
 /* else */
 return 1.0f;
}

/* x = sqrt(u^2+v^2), box in low,high 
  add outer and inner taper (planck-taper window)
  [t1...,t2 (inner) ,t3...,t4 (outer)]  
  */
static float
box_weight(float x,float t1, float t2, float t3, float t4) {
 if (x<=t1 || x>=t4) { return 0.0f; }
 if (x>t1 && x<t2) { 
   float z=(t2-t1)*(1.0f/(x-t1) + 1.0f/(x-t2));
   return (1.0f/(expf(z)+1.0f));
 }
 if (x>t3 && x<t4) { 
   float z=(t3-t4)*(1.0f/(x-t3) + 1.0f/(x-t4));
   return (1.0f/(expf(z)+1.0f));
 }
 /* else */
 return 1.0f;
}


/* build hash table, put each data point u,v and -u,-v into a bucket */
static void *
buildhash_threadfn(void *data) {
  thread_data_gridhash_t *t=(thread_data_gridhash_t*)data;

  uvlist *ll,*llnew;
  unsigned long int *key;
  gdpoint *idx;

  for(unsigned long int i = t->startrow; i < t->startrow+t->Nrows; i++) {
    if (!t->darr[i].flag) {
         /* if not flagged set initial weight to 1 */
          t->wold[i]=t->wnew[i]=1.0f;
    /* calculate right bucket index */
          float tempu = t->darr[i].u;
          float tempv = t->darr[i].v;
          float tempw = t->darr[i].w;
         /* do all computations using float, till the last moment */
         /* scale to [-0.5,0.5] */
          float uf=tempu*t->uvscale;
          float vf=tempv*t->uvscale;
          float wf=tempw*t->uvscale;
         /* scale to pixel values (never will be > 1/2 image size) */
          float ui=(uf*(float)t->Nx); /*  width */
          float vi=(vf*(float)t->Ny); /*  width */
          float poffX=-1.0f;
          float poffY=-1.0f;
        /* use relations  y+vi=(Ny/2-1) and x=ui+Nx/2-1 */
          int x=(int)round(ui+0.5f*t->Nx+poffX);
          int y=(int)round(-vi+0.5f*t->Ny+poffY);
         if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
           int bx=x;
           int by=y;
           /* bucket index is bx,by */
           unsigned long int bi=bx*t->B+by;
           if ((idx= (gdpoint*)malloc(sizeof(gdpoint)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
           }
           idx->id=i;
           idx->u=uf;
           idx->v=-vf;
           idx->w=wf;

           /* thread safe code */
           pthread_mutex_lock(t->writelock_hash);
           ll=(uvlist*)g_hash_table_lookup(t->ht,&bi);
           if (!ll) { /* new bucket */
              if ((key = (unsigned long int*)malloc(sizeof(unsigned long int)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
              }
              *key=bi;
              if ((llnew=(uvlist*)malloc(sizeof(uvlist)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
              }
              llnew->pix=NULL;
              llnew->P=0;
              llnew->id=NULL;
              llnew->u=NULL;
              llnew->v=NULL;
              llnew->w=NULL;
              llnew->pix=g_list_prepend(llnew->pix,idx);
              g_hash_table_insert(t->ht,(gpointer)key,(gpointer)llnew);
           } else {
             /* old bucket */
             ll->pix=g_list_prepend(ll->pix,idx);
           }
           pthread_mutex_unlock(t->writelock_hash);
         }
         /* also add the -ve point */
         x=-(int)round(-(-ui+0.5f*t->Nx+poffX));
         y=-(int)round(-(vi+0.5f*t->Ny+poffY));
         if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
           int bx=x;
           int by=y;
           /* bucket index is bx,by */
           unsigned long int bi=bx*t->B+by;
           if ((idx= (gdpoint*)malloc(sizeof(gdpoint)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
           }
           idx->id=i;
           idx->u=-uf;
           idx->v=vf;
           idx->w=-wf;

           /* thread safe code */
           pthread_mutex_lock(t->writelock_hash);
           ll=(uvlist*)g_hash_table_lookup(t->ht,&bi);
           if (!ll) { /* new bucket */
              if ((key = (unsigned long int*)malloc(sizeof(unsigned long int)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
              }
              *key=bi;
              if ((llnew=(uvlist*)malloc(sizeof(uvlist)))==0) {
                fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
                exit(1);
              }
              llnew->pix=NULL;
              llnew->P=0;
              llnew->id=NULL;
              llnew->u=NULL;
              llnew->v=NULL;
              llnew->w=NULL;
              llnew->pix=g_list_prepend(llnew->pix,idx);
              g_hash_table_insert(t->ht,(gpointer)key,(gpointer)llnew);
           } else {
             /* old bucket */
             ll->pix=g_list_prepend(ll->pix,idx);
           }
           pthread_mutex_unlock(t->writelock_hash);
         }
    }
  }
  return NULL;
}


/* same as above, only using a Gaussian as convolution kernel */
static float
convolve_neighbours_gauss(int P,unsigned long int *id,float *u,float *v,float *wt, float uf, float vf, float deltaU) {
    float *deluv,*delw;
    if ((deluv=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((delw=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    float invdeltaU=2.0f/(M_PI*deltaU);
    /* vectorizabe loop */
    for (int ci=0; ci<P; ci++) {
      float tmp1,tmp2;
      tmp1=(u[ci]-uf)*invdeltaU;
      tmp1*=tmp1;
      tmp2=(v[ci]+vf)*invdeltaU;
      tmp2*=tmp2;
      /* -(u-u')^2-(v-v')^2 */
      deluv[ci]=(-tmp1-tmp2);
    }
    for (int ci=0; ci<P; ci++) {
      deluv[ci]=expf(deluv[ci]);
      delw[ci]=wt[id[ci]];
    }
    float sump=my_fdot(P,deluv,delw);
    
    free(deluv);
    free(delw);
    return sump;
}

/* convolve P points (u,v,w)  positions, wt : weight
   uf,vf: zero location 
   zero location for w is dependent on uf,vf */
static float
convolve_neighbours_3d(int P,unsigned long int *id,float *u,float *v,float *w,float *wt, float uf, float vf,float deltaU) {
    float *deluv,*delw;
    if ((deluv=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((delw=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    /* coordinate for the cone u^2+v^2 */
    float uvcone=1.0f*sqrtf(uf*uf+vf*vf);
    float invdeltaU=2.0f/(M_PI*deltaU);
    /* vectorizabe loop */
    for (int ci=0; ci<P; ci++) {
      float tmp1,tmp2,tmp3;
      tmp1=(u[ci]-uf)*invdeltaU;
      tmp1*=tmp1;
      tmp2=(v[ci]+vf)*invdeltaU;
      tmp2*=tmp2;
      /* -(u-u')^2-(v-v')^2-(w-w')^2
          w'=alpha sqrt(u^2+v^2) */
      tmp3=(uvcone-fabsf(w[ci]));
      tmp3*=tmp3;
      deluv[ci]=-tmp1-tmp2-tmp3;
    }
    for (int ci=0; ci<P; ci++) {
      deluv[ci]=expf(deluv[ci]);
      delw[ci]=wt[id[ci]];
    }
    float sump=my_fdot(P,deluv,delw);
    
    free(deluv);
    free(delw);
    return sump;
}

/* convolve P points (u,v)  positions, wt : weight
   uf,vf: zero location + half pixel width (deltaU)
   kernel: NpxNp convolution kernel
   uvgrid: Np+2 grid points
   Np: kernel size
 */
static float 
convolve_neighbours_pswf(int P,unsigned long int *id,float *u,float *v,float *wt, float uf, float vf, float deltaU, float *kernel, float *uvgrid, int Np) {
    float *deluv,*delw;
    if ((deluv=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((delw=(float*)malloc((size_t)P*sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }

//printf("uf,vf=%f,%f\n",uf,-vf);
    for (int ci=0; ci<P; ci++) {
      float tmp1,tmp2;
      tmp1=fabsf(-u[ci]+uf); /* use +ve value since PSWF is symmetric */
      tmp2=fabsf(v[ci]+vf); /* note: vf is read from pixel list itself, so no sign change */
      deluv[ci]=conv_plane_eval(Np, uvgrid, kernel, tmp1, tmp2);

      delw[ci]=wt[id[ci]];
//printf("uf,vf=%f,%f ui,vi=%f,%f delta=%f,%f pswf=%f\n",uf,vf,u[ci],v[ci],tmp1,tmp2,deluv[ci]);
//printf("W=%f pswf=%f\n",delw[ci],deluv[ci]);
    }
    float sump=my_fdot(P,deluv,delw);

    free(deluv);
    free(delw);
    return sump;
}

static void *
updateweight_threadfn(void *data) {
  thread_data_gridhash_t *t=(thread_data_gridhash_t*)data;

  uvlist *ll;
  float invscale=1.0f/t->uvscale; 
  for(unsigned long int i = t->startrow; i < t->startrow+t->Nrows; i++) {
    if (!t->darr[i].flag) {
//printf("i=%ld u,v %f,%f wold=%f\n",i,t->darr[i].u,t->darr[i].v,t->wold[i]);
    /* calculate right bucket index */
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
         if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
           float sump=0.0f;
           int xlow=(x-t->M>=0?x-t->M:0);
           int xhigh=(x+t->M>t->Nx?t->Nx:x+t->M);
           int ylow=(y-t->M>=0?y-t->M:0);
           int yhigh=(y+t->M>t->Ny?t->Ny:y+t->M);
           for (int bx=xlow; bx<xhigh; bx++ ) 
           for (int by=ylow; by<yhigh; by++ ) {
//printf("U,V %f,%f (u,v) (%f,%f) (x,y) (%d,%d) -> (bx,by) (%d,%d) [xlow,xhigh] (%d,%d) [ylow,yhigh] (%d,%d)\n",tempu,tempv,uf,vf,x,y,bx,by,xlow,xhigh,ylow,yhigh);
           /* bucket index is bx,by, also neighbour buckets */
           unsigned long int bi=bx*t->B+by;
           pthread_mutex_lock(t->writelock_hash);
           ll=(uvlist*)g_hash_table_lookup(t->ht,&bi);
           pthread_mutex_unlock(t->writelock_hash);
            if (ll) { /* found neighbour pixel list */
             //if (ll->convmode==CONV_MODE_GAUSS2D) {
             // sump+=convolve_neighbours_gauss(ll->P,ll->id,ll->u,ll->v,t->wold,uf,vf,t->deltaU);
             //} else {
              //sump+=convolve_neighbours_3d(ll->P,ll->id,ll->u,ll->v,ll->w,t->wold,uf,vf,t->deltaU);
              sump+=convolve_neighbours_pswf(ll->P,ll->id,ll->u,ll->v,t->wold,uf,vf,t->deltaU,t->kernel,t->uvgrid,t->Np);
             //}
            }
           }
           /* W_k+1 <= (W_k x G_k) / (W_k \odot C_k) */
           /* G_k = 1 for uniform weights */
           /* since u,v are normalized to [-0.5,0.5], multiply by 2*800 */
           if (t->convmode==CONV_MODE_NCP) { 
            t->wnew[i]=t->wold[i]*ncp_weight(sqrtf(uf*uf+vf*vf)*invscale)/(sump+TOL);
           } else { /* default is uniform weight */
            t->wnew[i]=t->wold[i]/(sump+TOL);
           }
//printf("i=%ld oldw=%f neww=%f\n",i,t->wold[i],t->wnew[i]);
         }
    }
  }
  return NULL;
}



static void *
normalizeweights_threadfn(void *data) {
  thread_data_gridhash_t *t=(thread_data_gridhash_t*)data;

  uvlist *ll;
  float *wvec;
  /* iterate over rows of the buckets, 
  that way most threads will get some work to do */
  for (unsigned int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
    for (int cj=t->startrow; cj<t->B; cj++) { /* only upper half is needed */
      unsigned long int bi=ci*t->B+cj;
      pthread_mutex_lock(t->writelock_hash);
      ll=(uvlist*)g_hash_table_lookup(t->ht,&bi);
      pthread_mutex_unlock(t->writelock_hash);
      if (ll) {
       if ((wvec=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
       }
       for (int nrow=0; nrow<ll->P; nrow++) {
         wvec[nrow]=t->wnew[ll->id[nrow]];
       }
       /* calculate the sum (assume positive) */
       float wsum1=my_fasum(ll->P, wvec);
       /* normalize such that sum of weights is always 1 */
       if (wsum1>0.0f) {
        my_fscal(ll->P, 1.0f/wsum1, wvec);
       }
       /* copy back to weight array */
       for (int nw=0; nw<ll->P; nw++) {
         t->wnew[ll->id[nw]]=wvec[nw];
       }
       free(wvec);
      }
    }
  }


  return NULL;
}

/* write back the updated weight to the data array,
   also write weight to the weight grid for imaging */
static void *
writeweight_threadfn(void *data) {
  thread_data_gridhash_t *t=(thread_data_gridhash_t*)data;
  unsigned long int *buffoff;
  int BL=DATA_BUF_LEN;
  if ((buffoff=(unsigned long int*)calloc((size_t)BL,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  int bfilled=0;


  for(unsigned long int i = t->startrow; i < t->startrow+t->Nrows; i++) {
    if (!t->darr[i].flag) {
       /* update weight of data point */
       t->darr[i].wt=t->wold[i];
       /* also write weight to data grid */
    /* calculate right bucket index */
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
         if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
           coffset=ifftshift_index(t->Nx,x)*(t->Ny)+ifftshift_index(t->Ny,y);
           if (bfilled<BL) {
               buffoff[bfilled++]=coffset;
           } else {
            /* empty buffer */
            pthread_mutex_lock(t->writelock_hash);
            for (int bidx=0; bidx<BL; bidx++) {
                t->wgrid[buffoff[bidx]]+=1.0f;
            }
            pthread_mutex_unlock(t->writelock_hash);
            bfilled=0;
            buffoff[bfilled++]=coffset;
           }
         }
         /* also add -u,-v point */
         x=(int)round(-ui+0.5f*t->Nx+poffX);
         y=(int)round(vi+0.5f*t->Ny+poffY);
         if (x>=0 && x<t->Nx && y>=0 && y<t->Ny) {
           coffset=ifftshift_index(t->Nx,x)*(t->Ny)+ifftshift_index(t->Ny,y);
           if (bfilled<BL) {
               buffoff[bfilled++]=coffset;
           } else {
            /* empty buffer */
            pthread_mutex_lock(t->writelock_hash);
            for (int bidx=0; bidx<BL; bidx++) {
                t->wgrid[buffoff[bidx]]+=1.0f;
            }
            pthread_mutex_unlock(t->writelock_hash);
            bfilled=0;
            buffoff[bfilled++]=coffset;
           }

         }
    }
  }
  /* write last buffer */
  if (bfilled>0) {
    /* empty buffer */
    pthread_mutex_lock(t->writelock_hash);
    for (int bidx=0; bidx<bfilled; bidx++) {
        t->wgrid[buffoff[bidx]]+=1.0f;
    }
    pthread_mutex_unlock(t->writelock_hash);
  }

  free(buffoff);
  return NULL;
}

int
weightuvdata_pipe_menon(iodata *darr, int Nrows, float uvscale, float deltaU, float lambda, float *wgrid, int Nx, int Ny, int Nt, int convmode,  int witermax, float *sumweights, float *sumweights2, float imscale, complex float *wkernel, float *uvxgrid, int Np, int M) {
   /* wkernel: extract 0 plane (real part) NpxNp pixels,
      grid values are given by Np+2 array uvxgrid */

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_gridhash_t *threaddata;
   pthread_mutex_t writelock_hash;

   GHashTable *bucket;

   float *wold,*wnew,*werr;

   /* divide data rows over threads */
   int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1,ci,Nthb;
   int Mx=Nx; /* number of bucket of width M pixles */

   /* setup threads */
   pthread_mutex_init(&writelock_hash, NULL);
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_gridhash_t*)malloc((size_t)Nt*sizeof(thread_data_gridhash_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   /* allocate weight vectors */
   if ((wnew=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((wold=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((werr=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }

   /* setup hash table */
   bucket=g_hash_table_new_full(g_int64_hash, g_int64_equal,destroy_key,destroy_value);
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
    threaddata[nth].deltaU=deltaU;
    threaddata[nth].lambda=lambda;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].B=Mx;
    threaddata[nth].M=M;
    threaddata[nth].writelock_hash=&writelock_hash;
    threaddata[nth].ht=bucket;
    threaddata[nth].wgrid=wgrid;
    threaddata[nth].wnew=wnew;
    threaddata[nth].wold=wold;
    threaddata[nth].imscale=imscale;
    threaddata[nth].convmode=convmode;
    
    pthread_create(&th_array[nth],&attr,buildhash_threadfn,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* now wait for threads to finish */
  /* also calculate the sum of weights */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }
printf("done building hash table\n");
 /* convert lists in hash table to arrays, destroy lists */
 uvlist *ll;
 GHashTableIter iter;
 gpointer key, value;
 g_hash_table_iter_init(&iter, bucket);
 while (g_hash_table_iter_next (&iter, &key, &value)){
  ll=(uvlist*)value;  
  ll->P=g_list_length(ll->pix);
//  printf("bucket %ld has %d pixels\n",*idx,ll->P);
  /* allocate memory for this bucket into u,v,wt arrays */
  if ((ll->id=(unsigned long int*)calloc((size_t)ll->P,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  if ((ll->u=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ll->v=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ll->w=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  GList *li;
  int k=0;
  for(li=ll->pix; li!=NULL; li=g_list_next(li)) {
     gdpoint *ii=(gdpoint*)li->data;
     ll->id[k]=ii->id;
     ll->u[k]=ii->u;
     ll->v[k]=ii->v;
     ll->w[k]=ii->w;
     g_free(ii);
     k++;
  }
  g_list_free(ll->pix);
  ll->pix=NULL;
 }

/* wkernel: extract 0 plane (real part) NpxNp pixels,
      grid values are given by Np+2 array uvxgrid */
  float *wzeroplane;
  if ((wzeroplane=(float*)calloc((size_t)Np*Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* copy real part of 0 plane, from complex data, using stride 1 */
  my_fcopy(Np*Np, (float *)wkernel, 2, wzeroplane, 1);

printf("starting weight update\n");
 /******** weight iteration ********************/
 float woldnorm=my_fnrm2(Nrows,wold); /* keep norm fixed at this value */
 float wnewnorm=woldnorm;

 /* 1) iterate over rows and calculate new weight (convolution) */
 for (int niter=0; niter<witermax; niter++) {
   NrowsNt=(Nrows+Nt-1)/Nt;
   ci=0;
   for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }

    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;
    threaddata[nth].kernel=wzeroplane;
    threaddata[nth].uvgrid=uvxgrid;
    threaddata[nth].Np=Np;

    pthread_create(&th_array[nth],&attr,updateweight_threadfn,(void*)(&threaddata[nth]));
    ci=ci+NrowsNt;
  }
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }
  /* werr=wold-wnew */
  my_fcopy(Nrows,wold,1,werr,1);
  my_faxpy(Nrows,wnew,1,-1.0f,werr,1);
  float werrnorm=my_fnrm2(Nrows,werr);
  printf("PM iter %d error ||Wold-Wnew||/||Wold|| %f\n",niter,werrnorm/wnewnorm);

  /* wold <- wnew, only if new weight is finite*/
  if (isfinite(werrnorm)) {
   wnewnorm=my_fnrm2(Nrows,wnew);
   my_fcopy(Nrows,wnew,1,wold,1);
  } else {
   fprintf(stderr,"PM weight not finite, stopping\n");
   break;
  }
 }

 free(wzeroplane);
 /* preserve norm of weights to original value */
 my_fscal(Nrows, woldnorm/wnewnorm, wnew);
  
 /* calculate sum of weights (all positive) */
 float sumwt=my_fasum(Nrows,wold);
 /* calculate norm sum(sqr(.)) */
 float sumwt2=woldnorm*woldnorm;
 printf("Sum of weights=%f\n",sumwt);
 printf("Sum of weights^2=%f\n",sumwt2);
 *sumweights=sumwt;
 *sumweights2=sumwt2;

 /* 2) copy back updated weights to data array, also write weight to 
    binary file */
 memset((void*)wgrid,0,sizeof(float)*Nx*Ny);
 NrowsNt=(Nrows+Nt-1)/Nt;
 ci=0;
 for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }

    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;

    pthread_create(&th_array[nth],&attr,writeweight_threadfn,(void*)(&threaddata[nth]));
    ci=ci+NrowsNt;
 }
 for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
 }

 pthread_attr_destroy(&attr);
 pthread_mutex_destroy(&writelock_hash);

 free(th_array);
 free(threaddata);

 /* write a summary of weights in each bucket */
 /*FILE *cfilep=fopen("bucket_weights","w+");
 g_hash_table_iter_init(&iter, bucket);
 while (g_hash_table_iter_next(&iter, &key, &value)){
  ll=(uvlist*)value;
  float meanu=0.0f;
  float meanv=0.0f;
  float mu=0.0f,mv=0.0f,sumw=0.0f;
  for (int nw=0; nw<ll->P; nw++) {
    mu+=ll->u[nw]*wnew[ll->id[nw]];
    mv+=ll->v[nw]*wnew[ll->id[nw]];
    sumw+=wnew[ll->id[nw]];
    meanu+=ll->u[nw];
    meanv+=ll->v[nw];
  }
  meanu/=(float)ll->P;
  meanv/=(float)ll->P;
  mu/=sumw;
  mv/=sumw;
  float uf=round(meanu*(float)Nx+0.5f*Nx-1.0f)*deltaU-0.5f;
  float vf=round(meanv*(float)Nx+0.5f*Nx-1.0f)*deltaU-0.5f;
  fprintf(cfilep,"%f %f %f %f %f %f\n",uf,vf,meanu,meanv,mu,mv);
 }
 fclose(cfilep);
*/

 free(wnew);
 free(wold);
 free(werr);



 /* destroy hash table */
 g_hash_table_destroy(bucket);
 return 0;
}



#ifdef HAVE_CUDA
/* GPU version */
int
weightuvdata_pipe_menon_gpu(iodata *darr, int Nrows, float uvscale, float deltaU, float lambda, float *wgrid, int Nx, int Ny, int Nt, int convmode,  int witermax, float *sumweights, float *sumweights2, float imscale, complex float *wkernel, float *uvxgrid, int Np, int M) {
   /* wkernel: extract 0 plane (real part) NpxNp pixels,
      grid values are given by Np+2 array uvxgrid */

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_gridhash_t *threaddata;
   pthread_mutex_t writelock_hash;

   GHashTable *bucket;

   float *wold,*wnew,*werr;

   /* divide data rows over threads */
   int NrowsNt=(Nrows+Nt-1)/Nt;
   int nth,nth1,ci,Nthb;
   int Mx=Nx; /* number of bucket of width M pixles */

   /* setup threads */
   pthread_mutex_init(&writelock_hash, NULL);
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

   if ((th_array=(pthread_t*)malloc((size_t)Nt*sizeof(pthread_t)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
   }
   if ((threaddata=(thread_data_gridhash_t*)malloc((size_t)Nt*sizeof(thread_data_gridhash_t)))==0) {
     fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   /* allocate weight vectors */
   if ((wnew=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((wold=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }
   if ((werr=(float*)calloc((size_t)Nrows,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
   }

   /* setup hash table */
   bucket=g_hash_table_new_full(g_int64_hash, g_int64_equal,destroy_key,destroy_value);
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
    threaddata[nth].deltaU=deltaU;
    threaddata[nth].lambda=lambda;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].B=Mx;
    threaddata[nth].M=M;
    threaddata[nth].writelock_hash=&writelock_hash;
    threaddata[nth].ht=bucket;
    threaddata[nth].wgrid=wgrid;
    threaddata[nth].wnew=wnew;
    threaddata[nth].wold=wold;
    threaddata[nth].imscale=imscale;
    threaddata[nth].convmode=convmode;
    
    pthread_create(&th_array[nth],&attr,buildhash_threadfn,(void*)(&threaddata[nth]));
    /* next baseline set */
    ci=ci+NrowsNt;
  }

  /* now wait for threads to finish */
  /* also calculate the sum of weights */
  for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
  }
printf("done building hash table\n");
 /* convert lists in hash table to arrays, destroy lists */
 uvlist *ll;
 GHashTableIter iter;
 gpointer key, value;
 g_hash_table_iter_init(&iter, bucket);
 while (g_hash_table_iter_next (&iter, &key, &value)){
  ll=(uvlist*)value;  
  ll->P=g_list_length(ll->pix);
//  printf("bucket %ld has %d pixels\n",*idx,ll->P);
  /* allocate memory for this bucket into u,v,wt arrays */
  if ((ll->id=(unsigned long int*)calloc((size_t)ll->P,sizeof(unsigned long int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }

  if ((ll->u=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ll->v=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  if ((ll->w=(float*)calloc((size_t)ll->P,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  GList *li;
  int k=0;
  for(li=ll->pix; li!=NULL; li=g_list_next(li)) {
     gdpoint *ii=(gdpoint*)li->data;
     ll->id[k]=ii->id;
     ll->u[k]=ii->u;
     ll->v[k]=ii->v;
     ll->w[k]=ii->w;
     g_free(ii);
     k++;
  }
  g_list_free(ll->pix);
  ll->pix=NULL;
 }

/* wkernel: extract 0 plane (real part) NpxNp pixels,
      grid values are given by Np+2 array uvxgrid */
  float *wzeroplane;
  if ((wzeroplane=(float*)calloc((size_t)Np*Np,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* copy real part of 0 plane, from complex data, using stride 1 */
  my_fcopy(Np*Np, (float *)wkernel, 2, wzeroplane, 1);

printf("starting weight update\n");
 /******** weight iteration - use GPUs ********************/
 float woldnorm=my_fnrm2(Nrows,wold); /* keep norm fixed at this value */
 float wnewnorm=woldnorm;
 gbpmgdata tpg;
 int Ntg=4*Nt;
 if ((tpg.status=(int*)malloc(sizeof(int)*2*Ntg))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }
 if ((tpg.carray=(void**)malloc(sizeof(void*)*2*Ntg))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }
 if ((tpg.startrow=(unsigned long int*)malloc(sizeof(unsigned long int)*2*Ntg))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }
 if ((tpg.Nrows=(unsigned long int*)malloc(sizeof(unsigned long int)*2*Ntg))==0) {
    fprintf(stderr,"no free memory\n");
    exit(1);
 }

 tpg.wkernel=wzeroplane;
 tpg.Np=Np;
 tpg.wold=wold;
 tpg.wnew=wnew;
 tpg.uvscale=uvscale;
 tpg.darr=(iodata*)darr;
 tpg.Nx=Nx;
 tpg.Ny=Ny;
 tpg.B=Mx;
 tpg.M=M;
 tpg.writelock_hash=&writelock_hash;
 tpg.ht=bucket;
 tpg.convmode=convmode;

 th_pmpipeline tp;
 init_pm_pipeline(&tp,&tpg,Ntg);
 sync_barrier(&(tp.gate1));
 for (nth=0; nth<2*Ntg; nth++) {
  tpg.status[nth]=PT_DO_AGPU;
 }
 sync_barrier(&(tp.gate2));
 sync_barrier(&(tp.gate1));
 for (nth=0; nth<2*Ntg; nth++) {
  tpg.status[nth]=PT_DO_NOTHING;
 }
 NrowsNt=(Nrows+2*Ntg-1)/(2*Ntg);
 ci=0;
 for (nth=0;  nth<2*Ntg && ci<Nrows; nth++) {
  if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
  } else {
     Nthb=Nrows-ci;
  }
  tpg.startrow[nth]=ci;
  tpg.Nrows[nth]=Nthb;
  ci=ci+NrowsNt;
 }
 sync_barrier(&(tp.gate2));

 /* 1) iterate over rows and calculate new weight (convolution) */
 for (int niter=0; niter<witermax; niter++) {
  sync_barrier(&(tp.gate1));
  for (nth=0; nth<2*Ntg; nth++) {
   tpg.status[nth]=PT_DO_WORK_GRID;
  }
  sync_barrier(&(tp.gate2));
  sync_barrier(&(tp.gate1));
  for (nth=0; nth<2*Ntg; nth++) {
   tpg.status[nth]=PT_DO_NOTHING;
  }
  sync_barrier(&(tp.gate2));
  /* werr=wold-wnew */
  my_fcopy(Nrows,wold,1,werr,1);
  my_faxpy(Nrows,wnew,1,-1.0f,werr,1);
  float werrnorm=my_fnrm2(Nrows,werr);
  printf("PM iter %d error ||Wold-Wnew||/||Wold|| %f\n",niter,werrnorm/wnewnorm);

  /* wold <- wnew, only if new weight is finite*/
  if (isfinite(werrnorm)) {
   /* normalize new weight vector */
   wnewnorm=my_fnrm2(Nrows,wnew);
   my_fcopy(Nrows,wnew,1,wold,1);

  } else {
   fprintf(stderr,"PM weight not finite, stopping\n");
   break;
  }
 }

 sync_barrier(&(tp.gate1));
 for (nth=0; nth<2*Ntg; nth++) {
  tpg.status[nth]=PT_DO_DGPU;
 }
 sync_barrier(&(tp.gate2));
 destroy_pm_pipeline(&tp,Ntg);
 free(tpg.status);
 free(tpg.carray);
 free(tpg.startrow);
 free(tpg.Nrows);
printf("destroyed texture\n");
 free(wzeroplane);

 /* preserve norm of the weights to original value */
 my_fscal(Nrows, woldnorm/wnewnorm, wnew);

 /* calculate sum of weights (all positive) */
 float sumwt=my_fasum(Nrows,wold);
 /* calculate norm sum(sqr(.)) */
 float sumwt2=woldnorm*woldnorm;
 printf("Sum of weights=%f\n",sumwt);
 printf("Sum of weights^2=%f\n",sumwt2);
 *sumweights=sumwt;
 *sumweights2=sumwt2;

 /* 2) copy back updated weights to data array, also write weight to 
    binary file */
 memset((void*)wgrid,0,sizeof(float)*Nx*Ny);
 NrowsNt=(Nrows+Nt-1)/Nt;
 ci=0;
 for (nth=0;  nth<Nt && ci<Nrows; nth++) {
    if (ci+NrowsNt<Nrows) {
     Nthb=NrowsNt;
    } else {
     Nthb=Nrows-ci;
    }

    threaddata[nth].startrow=ci;
    threaddata[nth].Nrows=Nthb;

    pthread_create(&th_array[nth],&attr,writeweight_threadfn,(void*)(&threaddata[nth]));
    ci=ci+NrowsNt;
 }
 for(nth1=0; nth1<nth; nth1++) {
   pthread_join(th_array[nth1],NULL);
 }

 pthread_attr_destroy(&attr);
 pthread_mutex_destroy(&writelock_hash);

 free(th_array);
 free(threaddata);

 /* write a summary of weights in each bucket */
/* FILE *cfilep=fopen("bucket_weights_gpu","w+");
 g_hash_table_iter_init(&iter, bucket);
 while (g_hash_table_iter_next(&iter, &key, &value)){
  ll=(uvlist*)value;
  float meanu=0.0f;
  float meanv=0.0f;
  float mu=0.0f,mv=0.0f,sumw=0.0f;
  for (int nw=0; nw<ll->P; nw++) {
    mu+=ll->u[nw]*wnew[ll->id[nw]];
    mv+=ll->v[nw]*wnew[ll->id[nw]];
    sumw+=wnew[ll->id[nw]];
    meanu+=ll->u[nw];
    meanv+=ll->v[nw];
  }
  meanu/=(float)ll->P;
  meanv/=(float)ll->P;
  mu/=sumw;
  mv/=sumw;
  float uf=round(meanu*(float)Nx+0.5f*Nx-1.0f)*deltaU-0.5f;
  float vf=round(meanv*(float)Nx+0.5f*Nx-1.0f)*deltaU-0.5f;
  fprintf(cfilep,"%f %f %f %f %f %f\n",uf,vf,meanu,meanv,mu,mv);
 }
 fclose(cfilep);
*/

 free(wnew);
 free(wold);
 free(werr);



 /* destroy hash table */
 g_hash_table_destroy(bucket);
 return 0;
}
#endif /* HAVE_CUDA */

static void *
weight_threadfn(void *data) {
  thread_data_weight_t *t=(thread_data_weight_t*)data;
  for(int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
    float uf=(ci-t->pix0)*t->deltaU;
    for (int cj=0; cj<t->Ny; cj++) { 
     float vf=(cj-t->pix0)*t->deltaU;
     unsigned long int idx=ifftshift_index(t->Nx,ci)*t->Ny+ifftshift_index(t->Ny,cj);
     float wt=0.0f;
     if (t->taperfunc==CONV_MODE_NCP) { 
      wt=ncp_weight(sqrtf(uf*uf+vf*vf));
     } else if (t->taperfunc==CONV_MODE_NCP_HD) {
      wt=ncp_weight_highres(sqrtf(uf*uf+vf*vf));
     } else if (t->taperfunc==CONV_MODE_50_250) {
      wt=box_weight(sqrtf(uf*uf+vf*vf),30.0f,50.0f,250.0f,270.0f);
     } 
     //printf("pix %d,%d, uv=%f,%f, wt=%f\n",ci,cj,uf,vf,wt);
     //printf("before psf=(%f,%f) dat=(%f,%f)\n",crealf(t->psfgrid[idx]),cimagf(t->psfgrid[idx]),crealf(t->uvgrid[idx]),cimagf(t->uvgrid[idx]));
     t->psfgrid[idx]*=wt;
     t->uvgrid[idx]*=wt;
    }
  }
  return NULL;
}


static void *
weight_threadfn_one(void *data) {
  thread_data_weight_t *t=(thread_data_weight_t*)data;
  for(int ci=t->startrow; ci<t->startrow+t->Nrows; ci++) {
    float uf=(ci-t->pix0)*t->deltaU;
    for (int cj=0; cj<t->Ny; cj++) { 
     float vf=(cj-t->pix0)*t->deltaU;
     unsigned long int idx=ifftshift_index(t->Nx,ci)*t->Ny+ifftshift_index(t->Ny,cj);
     
     float wt=0.0f;
     if (t->taperfunc==CONV_MODE_NCP) { 
      wt=ncp_weight(sqrtf(uf*uf+vf*vf));
     } else if (t->taperfunc==CONV_MODE_NCP_HD) {
      wt=ncp_weight_highres(sqrtf(uf*uf+vf*vf));
     } else if (t->taperfunc==CONV_MODE_50_250) {
      wt=box_weight(sqrtf(uf*uf+vf*vf),30.0f,50.0f,250.0f,270.0f);
     }
     t->uvgrid[idx]*=wt;
    }
  }
  return NULL;
}




/* taper both uvgrid and psfgrid
 if psfgrid=0, only taper uvgrid */
int
tapergrid(complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, int pix0, float deltaU, int taperfunc, int Nt) {

   pthread_attr_t attr;
   pthread_t *th_array;
   thread_data_weight_t *threaddata;

   /* divide data rows over threads */
   long int Nrows=Nx;
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
    threaddata[nth].uvgrid=uvgrid;
    threaddata[nth].psfgrid=psfgrid;
    threaddata[nth].pix0=pix0;
    threaddata[nth].Nx=Nx;
    threaddata[nth].Ny=Ny;
    threaddata[nth].deltaU=deltaU;
    threaddata[nth].taperfunc=taperfunc;
   

    //printf("thread %d work on %ld rows, starting from %ld\n",nth, Nthb, ci);
    if (psfgrid) {
     pthread_create(&th_array[nth],&attr,weight_threadfn,(void*)(&threaddata[nth]));
    } else {
     pthread_create(&th_array[nth],&attr,weight_threadfn_one,(void*)(&threaddata[nth]));
    }
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
