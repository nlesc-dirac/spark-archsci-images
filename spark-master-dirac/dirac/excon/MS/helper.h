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

#ifndef __HELPER_H__
#define __HELPER_H__
#include <iostream>
#include <fstream>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Vector.h>
#include <casa/Arrays/IPosition.h>

#include <gridder.h> /* for fcomp definition */
using namespace std;

#define ARCSEC_TO_RADIANS M_PI/648000.0f

/* write weight map */
extern void
writeFITSW(char *fname, float *ptr, unsigned int w, unsigned int h, double delta, double *ref, double freq, double deltaf, float bscale, int stokes);
/* stokes: 1, 2, 3, 4 : I,Q,U,V */
extern void 
writeFITSIMG(char *fname, float *ptr, unsigned int w, unsigned int h, double delta, double *ref, double freq, double deltaf, int stokes, double bmaj, double bmin, double bpa, float sumweight);

namespace Cmd
{
    extern bool W_Projection; /* true: enable w projection */
    extern int AvgChannels; /* 0: MFS, 1: average channels, 2: channel imaging */
    extern char *TableName; /* MS name */
    extern const char *File_uvgrid; /* file names for binray files */
    extern const char *File_imgrid;
    extern const char *File_wgrid;
    extern const char *file_qualifier; /* unique qualifier to image files names */
    extern casa::String DataField; /* which column to image DATA/CORRECTED_DATA */
    extern float UVscale; /* scale of uv coords to fit [-0.5,0.5] */
    extern float ImRes; /* pixel size, arcseconds */
    extern float PointingRef[2]; /* phase shifting direction */
    extern bool HasPointing; /* true: different phase center than one in MS */
    extern int Nx; /* image width, pixels */
    extern int Ny; /* image height, pixels */
    extern int Stacks; /* no of W projection planes */
    extern char Stokes; /* which stokes component */
    extern double TimeInterval; //Time interval (in minutes) for the epoch imager.
    extern int Mc; /* conv. kernel bandwidth, in uv pixels */
    extern int variable_M; /* if 1, vary conv. kernel bandwidth with freq */
    extern int Nt; /* no of worker threads */
    extern int Npc; /* conv. kernel size in pixels Npc x Npc */
    extern float Nzc; /* zero padding for conv kernel (how many Npcs) padded size Npc(1+2*Nzc)  */
    extern float padding; /* padding factor > 1 for image */
    extern float robustval; /* robust weighting parameter */
    extern int imgmode; /* 0,1,2,3: 0: I, 1:I+extra, 2:IQUV, 3:IQUV+extra will create weights, PSF FITS files */
    extern float conv_cutoff; /* convlutional kernel truncation, as a fraction of peak value */
    extern float apod_cutoff; /* cutoff for apodization correction */
    extern int wgrid_step; /* how to sample w values, 0: sqrt(), 1: log(), 2 : linear,  3: w^(expW), 4 : user supplied text file */
    extern float expW; /* only used if wgrid_step is 3 */
    extern const char *File_wstep; /* file name for reading no. of w planes  and w steps  */
    extern int use_gpu; /* 0: CPU gridding, 1: GPU gridding */
    extern int create_imaginary; /* if 1, imaginary parts of all images/grids created */
    extern int gpu_threads; /* no. of worker threads per GPU */
    extern int w_snapshot; /* if >0 enable W snapshot */

    extern int weightmode; /* which weighting scheme to use ? */
    extern int pm_maxiter; /* iterations for Pipe-Menon density compensation */
    extern float pm_imscale; /* scale matched image weighting scale, radians FWHM of Gaussian */
    extern int pm_convmode; /* which convolution mode to use 2D or 3D */

    extern int buckets; /* no of buckets to use to sort data */
    void ParseCmdLine(int ac, char**av);

    extern int convkernel; /* if >0, create convolution kernel as a FITS file */
}

namespace Data
{

    //data loaded from the MS.
    extern iodata **data;
    extern float minW; /* min max W */
    extern float maxW;
    extern float maxUV[3]; /* max U, max V, sqrt(max(U^2+V^2)) */
    extern double refDir[2];
    extern int numChannels; 
    extern float chanFreq;
    extern float chanBW;
    extern float waveLen;
    extern unsigned int numRows;
    extern float min_uvcut;
    extern float max_uvcut; /* limits to select data based on uv range (lambda) */
    extern float max_wcut; /* limits max W term (lambda) */

    /* for rescaling data due to leverage */
    extern float rescale_uvcut; /* all baselines <= this value will be rescaled */
    extern float rescale_factor; /* rescaling factor */

    /* for clippling data XX,XY,YX,YY upper limits */
    extern float maxXY[4];
    extern bool clipval;
   
    /* for flagging using sigma XX,XY,YX,YY */
    extern float sigXY[4];
    extern bool flagsigma;

    /* for individual channel imaging */
    extern float *chanFreqs;
    extern float *chanBWs;
}

#endif //__HELPER_H__
