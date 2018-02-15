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

#include "data.h"
#include "helper.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <largefft.h>
#include <gridder.h>


/****************************** pthgridder.cpp ********************/
/* darr: Nrows x 1  of data from MS 
   uvscale: scale uv coords to fit [-0.5,0.5]
   lambda: wavelength
   deltaU: uv pixel size
   uvgrid: array to store gridded data : Nx*Ny complex float
   psfgrid: array to store PSF (gridded 1) : Nx*Ny complex float
   Nx,Ny : grid size (also image size)
   wparr: W plane values (with guard values) : Nw+2
   Nw: how many W planes
   wpsupport (X,Y): conv. kernel support for each W plane : Nw 
   wkernel: conv. kernels : each kernel NpxNp, Nw planes, complex float
   uvxgrid: array for pixel u,v values, with guard valuesm Np+2
   Np: conv. kernel width in uv pixels
   Nt: no of threads 
   
   Note about weight calculation
   for each data point, and pixel i in the grid that has a contribution
   uvgrid +=conv_kernel*data
   psfgrid +=(conv_kernel)
   
   Finally, when imaging
   wgrid_k <= |psfgrid_k|
   imgmode: grid I or also QUV
*/
extern int
griddata(iodata *darr, int Nrows, float uvscale, float lambda, float deltaU, complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, float *wparr, int Nw, int *wpsupportX, int *wpsupportY, complex float *wkernel, float *uvxgrid, int Np,int Nt, int imgmode);



typedef struct thread_data_weight_t_ {
  long int startrow,Nrows;
  float *wgrid;
  complex float *psfgrid;
  complex float *uvgrid;
  int pix0;
  int Nx,Ny;
  float deltaU;
  int taperfunc;
} thread_data_weight_t;


/* weight the gridded data and PSF with calculated weights
*/
extern int
weightdata(complex float *uvgrid, complex float *psfgrid, float *wgrid, int Nx, int Ny, float robust, int Nt);


/****************************** pthweighter.cpp ********************/
/* weight the uv data */
/* note: only one grid cell is mapped to one uv data point, 
  so convolution kernel support is 1 */
/* sumweights: sum of weights, sumweights2: sum of weights^2 */
extern int
weightuvdata(iodata *darr, int Nrows, float uvscale, float lambda, float robustval, float *wgrid, int Nx, int Ny, int Nt, int weightmode, float *sumweights, float *sumweights2);


/****************************** pthpipe_menon.cpp ********************/
/* convmode : convolution kernel support in pixels */
/* witermax: iterations for the weight update */
/* imscale: scale of Gaussian exp(-x^2/scale^2) (in radians) to match image weights */
/* wkernel: extract NpxNp w=zero plane from this 
   grid values in Np+2 array uvxgrid (for binary search) 
  M : con kernel support is 2M */
extern int
weightuvdata_pipe_menon(iodata *darr, int Nrows, float uvscale, float deltaU, float lambda, float *wgrid, int Nx, int Ny, int Nt, int convmode, int witermax, float *sumweights, float *sumweights2, float imscale, complex float *wkernel, float *uvxgrid, int Np, int M);

extern int
weightuvdata_pipe_menon_gpu(iodata *darr, int Nrows, float uvscale, float deltaU, float lambda, float *wgrid, int Nx, int Ny, int Nt, int convmode, int witermax, float *sumweights, float *wumweights2, float imscale, complex float *wkernel, float *uvxgrid, int Np, int M);


/* taper both uvgrid and psfgrid (if not zero)
  note deltaU is correctly scaled to get back the correct u,v coords
  if psfgrid=0, only taper uvgrid
  tapefunc: CONV_MODE */
extern int
tapergrid(complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, int pix0, float deltaU, int taperfunc, int Nt);


/*************** snapshot.cpp ***********************/
extern int
project_wplane(iodata *darr, int Nrows, int Nt, double *wpa,double *wpb, float *wmax, float *wmin);
