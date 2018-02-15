/* 
 *
   Copyright (C) 2015 Sarod Yatawatta <sarod@users.sf.net>  
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

/* we split points in u more often than in v
   because data is stored in col major order
   the freq. of splitting in u is SPLIT_LEVEL-1: 1
  */
#ifndef SPLIT_LEVEL
#define SPLIT_LEVEL 2
#endif
/* comparison function to compare two nodes */
static gint 
compare_nodes(gconstpointer a, gconstpointer b)
{
 ptreenode *p1=(ptreenode *)a;
 ptreenode *p2=(ptreenode *)b;

 return ((p1->N <= p2->N)? 1 : -1);
}


/* init tree */
/* input : 
  tree : p
  data : N x 1 (x,y) data array 
  
  will add root node to tree, create index array [1..N]
  add this to leaf list
*/
int
init_ptree(ptree *p, iodata *data, int N) {

  ptreenode *root;
  if ((root=(ptreenode *)calloc(1,sizeof(ptreenode)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  root->id=0;
  root->level=0;
  root->N=N;
  if ((root->index=(int*)calloc(root->N,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  /* go through the data and add indices, find limits */
  root->x_min=root->y_min=1e6;
  root->x_max=root->y_max=-1e6;
  int ci;
  int total=0;
  for (ci=0; ci<N; ci++) {
   if (!data[ci].flag) { /* only include unflagged data */
    root->index[total++]=ci;
    if (root->x_min > data[ci].u) {
     root->x_min=data[ci].u;
    }
    if (root->x_max < data[ci].u) {
     root->x_max=data[ci].u;
    }
    if (root->y_min > data[ci].v) {
     root->y_min=data[ci].v;
    }
    if (root->y_max < data[ci].v) {
     root->y_max=data[ci].v;
    }
   }
  }
  /* resize memory to only include unflagged data indices */
  if ((root->index=(int*)realloc(root->index,total*sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  root->N=total;

  /* add the node to the lists */
  p->leaves=NULL;
  p->leaves=g_list_insert_sorted(p->leaves,root,compare_nodes);

  p->count=1;
  if ((p->unflagged=(int*)calloc(root->N,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
  memcpy(p->unflagged,root->index,sizeof(int)*root->N);
  p->D=root->N;
#ifdef DEBUG
  printf("unflagged data %d out of %d total\n",p->D,N);
#endif

  return 0;
}



/* deallocate all memory and delete tree */
int 
destroy_ptree(ptree *p) {

  /* go through each node */
  GList *li;
  for(li=p->leaves; li!=NULL; li=g_list_next(li)) {
     ptreenode *nn= li->data;
     /* free index array */
     free(nn->index);
     free(nn);
  }
  g_list_free(p->leaves);
  free(p->unflagged);

  return 0;
}


/* find the node with the largest index array, split it into two */
/* data: Nx1 data array */
int
split_ptree(ptree *p, iodata *data, int N) {

 /* get node with highest index array */
 GList *li=g_list_first(p->leaves);
 ptreenode *parent=(ptreenode *)li->data;
 /* remove this from list */
 p->leaves=g_list_remove_link(p->leaves,li);
 g_list_free(li);

 /* now split parent to 2 */
 ptreenode *child1,*child2;
 if ((child1=(ptreenode *)calloc(1,sizeof(ptreenode)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 if ((child2=(ptreenode *)calloc(1,sizeof(ptreenode)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 child1->id=p->count++;
 child1->level=parent->level+1;
 child1->N=parent->N;
 if ((child1->index=(int*)calloc(parent->N,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 child2->id=p->count++;
 child2->level=parent->level+1;
 child2->N=N;
 if ((child2->index=(int*)calloc(parent->N,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
 }
 /* go through the data and add indices, find limits */
 child1->x_min=child1->y_min=1e6;
 child2->x_min=child2->y_min=1e6;
 child1->x_max=child1->y_max=-1e6;
 child2->x_max=child2->y_max=-1e6;
 double boundary;
 if (parent->level%SPLIT_LEVEL) {
   boundary=(parent->x_min+parent->x_max)*0.5;
#ifdef DEBUG
   printf("splitting on x boundary %lf\n",boundary);
#endif
 } else {
   boundary=(parent->y_min+parent->y_max)*0.5;
#ifdef DEBUG
   printf("splitting on y boundary %lf\n",boundary);
#endif
 }
 int ci,n1,n2;
 n1=n2=0;
 for (ci=0; ci<parent->N; ci++) {
   int ii=parent->index[ci];
   /* split on x or y depending on level */
   if (parent->level%SPLIT_LEVEL) {
    /* split on x coordinate */
    if (data[ii].u<=boundary) {
     child1->index[n1++]=ii;
     child1->x_min=(child1->x_min>data[ii].u?data[ii].u:child1->x_min);
     child1->y_min=(child1->y_min>data[ii].v?data[ii].v:child1->y_min);
     child1->x_max=(child1->x_max<data[ii].u?data[ii].u:child1->x_max);
     child1->y_max=(child1->y_max<data[ii].v?data[ii].v:child1->y_max);
    } else {
     child2->index[n2++]=ii;
     child2->x_min=(child2->x_min>data[ii].u?data[ii].u:child2->x_min);
     child2->y_min=(child2->y_min>data[ii].v?data[ii].v:child2->y_min);
     child2->x_max=(child2->x_max<data[ii].u?data[ii].u:child2->x_max);
     child2->y_max=(child2->y_max<data[ii].v?data[ii].v:child2->y_max);
    }
   } else {
    /* split on y coordinate */
    if (data[ii].v<=boundary) {
     child1->index[n1++]=ii;
     child1->x_min=(child1->x_min>data[ii].u?data[ii].u:child1->x_min);
     child1->y_min=(child1->y_min>data[ii].v?data[ii].v:child1->y_min);
     child1->x_max=(child1->x_max<data[ii].u?data[ii].u:child1->x_max);
     child1->y_max=(child1->y_max<data[ii].v?data[ii].v:child1->y_max);
    } else {
     child2->index[n2++]=ii;
     child2->x_min=(child2->x_min>data[ii].u?data[ii].u:child2->x_min);
     child2->y_min=(child2->y_min>data[ii].v?data[ii].v:child2->y_min);
     child2->x_max=(child2->x_max<data[ii].u?data[ii].u:child2->x_max);
     child2->y_max=(child2->y_max<data[ii].v?data[ii].v:child2->y_max);
    }
   }
 }
#ifdef DEBUG
 printf("split node %d to %d and %d\n",parent->N,n1,n2);
 printf("boundary [(%lf,%lf),(%lf,%lf)]->[(%lf,%lf),(%lf,%lf)]-[(%lf,%lf),(%lf,%lf)]\n",parent->x_min,parent->x_max,parent->y_min,parent->y_max,child1->x_min,child1->x_max,child1->y_min,child1->y_max,child2->x_min,child2->x_max,child2->y_min,child2->y_max);
#endif
 /* reclaim unused memory */
 if (n1>0) {
  if ((child1->index=(int*)realloc(child1->index,n1*sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
 } else {
  free(child1->index);
  child1->index=NULL;
 }
 child1->N=n1;
 if (n2>0) {
  if ((child2->index=(int*)realloc(child2->index,n2*sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
  }
 } else {
  free(child2->index);
  child2->index=NULL;
 }
 child2->N=n2;

 /* insert children to tree */
 if (n1>0) {
  p->leaves=g_list_insert_sorted(p->leaves,child1,compare_nodes);
 } else {
  free(child1);
 }
 if (n2>0) {
  p->leaves=g_list_insert_sorted(p->leaves,child2,compare_nodes);
 } else {
  free(child2);
 }

 /* free parent storage */
 free(parent->index);
 parent->index=NULL;
 free(parent);

 return 0;
}


/* hash table routines */
/* key destroy function */
static void
destroy_key(gpointer data) {
 free((int*)data);
}
/* value destroy function */
static void
destroy_value(gpointer data) {
 free((int*)data);
}


/* in-place sort the data to the order given by the 
   ptree, using a hash table to cache */
/* data: Nx1 data array */
int
sort_data_with_ptree(ptree *p, iodata *data, int N) {

  GHashTable *ht=g_hash_table_new_full(g_int_hash, g_int_equal,destroy_key,destroy_value);

  int gid=0;
  GList *li;
  iodata tmp;
  int *key,*value,*k,ci;
  /* go through the nodes */
  for(li=p->leaves; li!=NULL; li=g_list_next(li)) {
     ptreenode *nn= li->data;
     for (ci=0; ci<nn->N; ci++) {
       /* now switch index[ci] value with gid  index[ci] <-> gid */
       
       /* first look up if this row (gid) is already switched
          before, if so, then 
          tmp <- val(gid)
          val(gid) <- val(hash(index[ci]) == val(k)
          val(hash(index[ci])) <- tmp
          remove hash: index[ci]->k
          add hash: gid -> k */
       k=(int*)g_hash_table_lookup(ht,&nn->index[ci]);
       if (k) {
        memcpy(&tmp,&data[p->unflagged[gid]],sizeof(iodata));
        memcpy(&data[p->unflagged[gid]],&data[*k],sizeof(iodata));
        memcpy(&data[*k],&tmp,sizeof(iodata));
        if ((key=(int*)calloc(1,sizeof(int)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
        }
        *key=gid;
        if ((value=(int*)calloc(1,sizeof(int)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
        }
        *value=*k;

        g_hash_table_remove(ht,&nn->index[ci]);
        g_hash_table_insert(ht,(gpointer)key,(gpointer)value);
        
       } else {
       /* tmp <- val(gid)
          val(gid) <- val(index[ci])
          val(index[ci]) <- tmp 
          hash: gid->index[ci] == hash(gid) <- index[ci] */
        memcpy(&tmp,&data[p->unflagged[gid]],sizeof(iodata));
        memcpy(&data[p->unflagged[gid]],&data[nn->index[ci]],sizeof(iodata));
        memcpy(&data[nn->index[ci]],&tmp,sizeof(iodata));
        if ((key=(int*)calloc(1,sizeof(int)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
        }
        *key=gid;
        if ((value=(int*)calloc(1,sizeof(int)))==0) {
         fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
         exit(1);
        }
        *value=nn->index[ci];
        g_hash_table_insert(ht,(gpointer)key,(gpointer)value);
       }
       gid++;
       
     }
#ifdef DEBUG
   printf("Hash table size %d\n",g_hash_table_size(ht));
#endif
  }

 g_hash_table_destroy(ht);
 return 0;
}

/* data: Nx1 array of datapoints
   M : no of buckets to use
*/
int
sort_data(iodata *data, int N, int M, int Nthreads) {
  ptree p;
  init_ptree(&p, data, N);
  int ci;
  for (ci=1; ci<M; ci++) {
   split_ptree(&p, data, N);
  }

  sort_data_with_ptree(&p, data, N);
  destroy_ptree(&p);
  return 0;
}
