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


/* this for mmap large files */
#define _LARGEFILE64_SOURCE

#include <stdio.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>

#include "largefft.h"



int
open_binary_file(const char *uvgridname, int *fid0, complex float **d, const char *imname, int *fid1, complex float **im, const char *wtgridname, int *fid2, complex float **wgt, int Nx, int Ny,int nx, int ny, int imgmode, int snapshot) {


  long int len_file;
  struct stat statbuf;
  
  /* flags to tune mmap */
  int mmap_flags=MAP_SHARED|MAP_NONBLOCK|MAP_POPULATE;
  //int mmap_flags=MAP_SHARED|MAP_HUGETLB|MAP_POPULATE;
  
  /* check if files already exist, and exit after printing a warning */
  if (access(uvgridname,F_OK)!=-1) {
    fprintf(stderr,"%s: %d: file %s already exists. Remove this file or see -A,-B,-C options.\n",__FILE__,__LINE__,uvgridname);
    exit(1);
  }
  if (access(imname,F_OK)!=-1) {
    fprintf(stderr,"%s: %d: file %s already exists. Remove this file or see -A,-B,-C options.\n",__FILE__,__LINE__,imname);
    exit(1);
  }
  if (access(wtgridname,F_OK)!=-1) {
    fprintf(stderr,"%s: %d: file %s already exists. Remove this file or see -A,-B,-C options.\n",__FILE__,__LINE__,wtgridname);
    exit(1);
  }
  /* open/create file */
  if ((*fid0 = open(uvgridname,O_RDWR | O_CREAT, S_IRWXU | S_IRGRP | S_IROTH| O_LARGEFILE)) < 0){
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  /* file size in bytes */
  if (imgmode==IMG_I0||imgmode==IMG_I) {
   len_file=Nx*Ny*sizeof(complex float);
  } else {
   len_file=4*Nx*Ny*sizeof(complex float);
  }
  /* get correct size */
  if (ftruncate(*fid0,len_file)!=0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  /* make sure size is ok */
  /* find the file size */
  if (fstat (*fid0,&statbuf) < 0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  //printf("asked %ld, got %ld\n",len_file,(long int)statbuf.st_size);
  /* map memory */
  *d= (complex float*)mmap(NULL, len_file, PROT_READ|PROT_WRITE, mmap_flags, *fid0, 0);
  if ( !(*d)) {
     fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
     exit(1);
  }
  /* set to zero */
  memset((void*)*d,0,len_file);

  if ((*fid1 = open(imname,O_RDWR | O_CREAT, S_IRWXU | S_IRGRP | S_IROTH)) < 0){
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  len_file=Nx*Ny*sizeof(complex float);
  /* get correct size */
  if (ftruncate(*fid1,len_file)!=0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  /* make sure size is ok */
  /* find the file size */
  if (fstat (*fid1,&statbuf) < 0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  //printf("asked %ld, got %ld\n",len_file,(long int)statbuf.st_size);
  /* map memory */
  *im= (complex float*)mmap(NULL, len_file, PROT_READ|PROT_WRITE, mmap_flags, *fid1, 0);
  if (!(*im)) {
     fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
     exit(1);
  }
  /* set to zero */
  memset((void*)*im,0,len_file);
  /* add more if need to store images */
  if (snapshot) {
   len_file+=nx*ny*sizeof(float);
  }
  if ((*fid2 = open(wtgridname,O_RDWR | O_CREAT, S_IRWXU | S_IRGRP | S_IROTH)) < 0){
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  /* get correct size */
  if (ftruncate(*fid2,len_file)!=0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  /* make sure size is ok */
  /* find the file size */
  if (fstat (*fid2,&statbuf) < 0) {
    fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
    exit(1);
  }
  //printf("asked %ld, got %ld\n",len_file,(long int)statbuf.st_size);
  /* map memory */
  *wgt= (complex float*)mmap(NULL, len_file, PROT_READ|PROT_WRITE, mmap_flags, *fid2, 0);
  if (!(*wgt)) {
     fprintf(stderr,"%s: %d: no file open\n",__FILE__,__LINE__);
     exit(1);
  }
  /* set to zero */
  memset((void*)*wgt,0,len_file);
  /* set right pointer to average image */

  return 0;
}  

int 
sync_binary_file(complex float *d, complex float *im, complex float *wgt, int Nx, int Ny, int imgmode) {
  /* sync to disk */
  long int count;
  if (imgmode==IMG_I0||imgmode==IMG_I) {
   count=Nx*Ny;
  } else {
   count=4*Nx*Ny;
  }
  msync(d, (size_t)count*sizeof(complex float), MS_SYNC );
  count=Nx*Ny;
  msync(im, (size_t)count*sizeof(complex float), MS_SYNC );
  msync(wgt, (size_t)count*sizeof(complex float), MS_SYNC );

  return 0;
}

int 
close_binary_file(const char *uvgridname, int fid0, complex float *d, const char *imname, int fid1, complex float *im,  const char *wtgridname, int fid2, complex float *wgt, int Nx, int Ny, int nx, int ny, int imgmode, int snapshot) {
  /* sync to disk */
  long int count;
  if (imgmode==IMG_I0||imgmode==IMG_I) {
   count=Nx*Ny;
  } else {
   count=4*Nx*Ny;
  }
  msync(d, (size_t)count*sizeof(complex float), MS_SYNC );
  munmap((void*)d, (size_t)count*sizeof(complex float));
  count=Nx*Ny;
  msync(im, (size_t)count*sizeof(complex float), MS_SYNC );
  munmap((void*)im, (size_t)count*sizeof(complex float));
  if (snapshot) {
   msync(wgt, (size_t)count*sizeof(complex float)+nx*ny*sizeof(float), MS_SYNC );
   munmap((void*)wgt, (size_t)count*sizeof(complex float)+nx*ny*sizeof(float));
  } else {
   msync(wgt, (size_t)count*sizeof(complex float), MS_SYNC );
   munmap((void*)wgt, (size_t)count*sizeof(complex float));
  }

  /* close files */
  close(fid0);
  close(fid1);
  close(fid2);

  /* delete files */
  unlink(uvgridname);
  unlink(imname);
  unlink(wtgridname);
  return 0;
}



int
reset_binary_file(complex float *d,  complex float *im,  complex float *wgt, int Nx, int Ny,int imgmode) {


  long int len_file;
  /* file size in bytes */
  if (imgmode==IMG_I0||imgmode==IMG_I) {
   len_file=Nx*Ny*sizeof(complex float);
  } else {
   len_file=4*Nx*Ny*sizeof(complex float);
  }
  /* set to zero */
  memset((void*)d,0,len_file);

  len_file=Nx*Ny*sizeof(complex float);
  /* set to zero */
  memset((void*)im,0,len_file);
  /* set to zero (not the imaging part) */
  memset((void*)wgt,0,len_file);

  return 0;
}
