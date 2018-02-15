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

#ifndef __CUGRIDDER_H__
#define __CUGRIDDER_H__
#include "cuda.h"
#include <cuComplex.h>

#include "gridder.h"

extern int
cuda_griddata(iodata *darr, int Nrows, float uvscale, float lambda, float deltaU, complex float *uvgrid, complex float *psfgrid, int Nx, int Ny, float expW, float *wparr, int Nw, int *wpsupportX, int *wpsupportY, int maxsupport, complex float *wkernel, int Np, int Nz, int Nt,int Ngpu, int imgmode);



#endif /* __CUGRIDDER_H__ */
