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
#include "pthgridder.h"
#include <fstream>
#include <vector> 
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#ifdef HAVE_CUDA
#include "cugridder.h"
#endif

#include <largefft.h>
using namespace std;

int main(int argc, char **argv) {
    /* reset random seed */
    srand(time(0));
    float ref[2]={0,0}; //reference angle for this table
    Cmd::ParseCmdLine(argc, argv);
    int epoch = 0;
    bool done = false; 
    bool convkernel_created=false;

    float uvLim,L,deltaU;
     //Create a table with the ms so we can directly image it if this is not epoch imaging.
    Table t(Cmd::TableName);
    //stuff required for epoch imaging.
    MeasurementSet ms(Cmd::TableName); 
    //backup MS name
    char *TableName1=strdup(Cmd::TableName);
    char *Base=basename(TableName1);
//    Block<int> sort(3);
//    sort[0] = MS::TIME; /* note: TIME,ANT1,ANT2 for no ms iteration */
//    sort[1] = MS::ANTENNA1;
//    sort[2] = MS::ANTENNA2;
    Block<int> sort(1);
    sort[0] = MS::TIME; /* note: only sort over TIME for ms iterator to work */
    MSIter msIter(ms,sort,Cmd::TimeInterval);
    msIter.origin();
    cout << "Pixel size (arcsec):"<<Cmd::ImRes<<" (rad): " << Cmd::ImRes * ARCSEC_TO_RADIANS << endl;
    L = 1.0f/(2.0f*Cmd::ImRes*ARCSEC_TO_RADIANS); 
    cout << "Longest supportable baseline ~ 1/delta(l) L (lambda) = " << L << endl;
    /* padding */
    int Npad=((int)((Cmd::padding-1.0f)*(float)Cmd::Nx))/2;
    int Nx=Cmd::Nx+2*Npad;
    int Ny=Nx; /* always square images */
    int M=Cmd::Mc; /* convolution support in uv pixels */
    int Np=Cmd::Npc; /* no of points PSWF is sampled in [-1,1] */
    /* zero padded conv. kernel length */
    int Nz=Np*(int)(1.0f+2.0f*Cmd::Nzc); 

    cout<<"Padded image size : "<<Nx<<"x"<<Ny<<endl;
    /* FOV covered by Np conv kernel pixels */
    float delta_l=Cmd::ImRes*ARCSEC_TO_RADIANS*Nx/(float)(Np);
    int Nthreads=Cmd::Nt;
    fftwf_plan fftp0,fftp1,fftp2,fftp3,fftp4;
    init_fftw(Nthreads);


    complex float *din;
    complex float *dout;
    complex float *wgtin;
    const char *uvgridname;
    const char *imgridname;
    const char *wtgridname;
    if (Cmd::File_uvgrid) {
     uvgridname=Cmd::File_uvgrid;
    } else {
     uvgridname="uvgrid.tmp";
    }
    if (Cmd::File_imgrid) {
     imgridname=Cmd::File_imgrid;
    } else {
     imgridname="imgrid.tmp";
    }
    if (Cmd::File_wgrid) {
     wtgridname=Cmd::File_wgrid;
    } else {
     wtgridname="weights.tmp";
    }
    int snapshotmode=0;
    float *avgImage=0;
    if( Cmd::TimeInterval > 0.0 ) {
     snapshotmode=1;
    }

    int uvgridid,imgridid,wtgridid;
    open_binary_file(uvgridname,&uvgridid,&din,imgridname,&imgridid,&dout,wtgridname, &wtgridid, &wgtin, Nx,Ny, Cmd::Nx, Cmd::Ny, Cmd::imgmode, snapshotmode);
    if (snapshotmode) {
     avgImage=(float*)&wgtin[Nx*Ny];
    }
    /* FFT orders 
     dout => wgtin for PSF
     din => dout for image
    */
    do_create_fftw_plan(&fftp0,dout,wgtin,Nx,Ny);
    do_create_fftw_plan(&fftp1,din,dout,Nx,Ny);
    unsigned long int woffset;
    if(Cmd::imgmode==IMG_IQUV0||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
      /* need 3 more plans for Q,U,V FFT */
      woffset=Nx*Ny;
      /* din => dout for image */
      do_create_fftw_plan(&fftp2,&din[woffset],dout,Nx,Ny);
      woffset+=Nx*Ny;
      do_create_fftw_plan(&fftp3,&din[woffset],dout,Nx,Ny);
      woffset+=Nx*Ny;
      do_create_fftw_plan(&fftp4,&din[woffset],dout,Nx,Ny);
    }


    double bmaj,bmin,bpa; /* evaluate bmaj,bmin in pixels, bpa in radians */
    double avbmaj=0.0,avbmin=0.0,avbpax=0.0,avbpay=0.0; /* to store average value */
    double sumweights_avg=0.0f;
    double sumweights2_avg=0.0f;
    char buff[4096];
    /************* snapshot iterations *****************/
    while(!done){
        //we are using TimeInterval to decide if we are in epoch imaging or not.
        if( snapshotmode ) {
            //snapshot
            if( msIter.more() ) {
                loadData(msIter.table(), Cmd::TableName, Cmd::AvgChannels, (Cmd::create_imaginary>0), false);
                msIter++;
            } else {
                //this is the last epoch.
                done = true;
                break;
            }
        cout<<endl<<"======== Making snapshot "<<epoch<<" ========"<<endl;
        } else {
            //we are doing normal imaging so only one iteration.
            done = true; 
            //use the normal table instead of Ms
            // only sort table is image is very big >12000 pix
            loadData(t, Cmd::TableName, Cmd::AvgChannels, (Cmd::create_imaginary>0), false);
        }
        cout<< "Freq (MHz): "<<Data::chanFreq/1e6f<<", BW (kHz): "<<Data::chanBW/1e3<<", Channels: "<<Data::numChannels<<endl;
        cout << "MaxU (m): " << Data::maxUV[0] << endl;
        cout << "MaxV (m): " << Data::maxUV[1] << endl;
        //calculate the UVscale so we can have the longest baseline in the image.
        uvLim=0.50f;
        /* scale to fit _supportable_ coords (m) from maxuv to [-uvLim,uvLim] */
        Cmd::UVscale = uvLim/L; 
        cout << "UVscale (for baselines in lambda): " << Cmd::UVscale << endl;
        cout << "Longest sqrt|u^2+v^2| in data (m): " << Data::maxUV[2] << endl;
        cout << "L/Data:MaxUV (~ 1 for optimum sampling) : " << L/Data::maxUV[2] << endl; //how much imaging parameters match the uv coverage
        deltaU=2.0f*uvLim/((float)Nx);
        cout<< "UV pixel size (lambda) : "<<deltaU/Cmd::UVscale<<" scaled: "<<deltaU<<endl;
        cout<<"W range (m): "<<Data::minW<<" : "<<Data::maxW<<" planes: "<<Cmd::Stacks<<endl;
        float amaxW=(fabs(Data::minW)>fabs(Data::maxW)?fabs(Data::minW):fabs(Data::maxW));
        amaxW/=Data::waveLen;
        ref[0] = Data::refDir[0];
        ref[1] = Data::refDir[1];
        cout << "Phase center";
        if (Cmd::HasPointing) {
         cout << " (shifted) : (" << ref[0] << "," << ref[1] << ")"<<endl;
        } else {
         cout << " : (" << ref[0] << "," << ref[1] << ")"<<endl;
        }

        for(int chan = 0; chan < (Cmd::AvgChannels==2? Data::numChannels: 1); chan++) {
        //for(int chan = 0; chan <1; chan++) 
            iodata *output;
            if (Cmd::AvgChannels==2) {
             output = &(Data::data[0][chan*Data::numRows]);
            } else {
             //only one column, even when not averaging data
             output = Data::data[0];
            }
            /* sort data if buckets>1 */
            if (Cmd::buckets>1) {
             cout<<"Sorting data using "<<Cmd::buckets<<" buckets."<<endl;
             sort_data(output,Data::numRows,Cmd::buckets,Nthreads);
            }
            /* if flagging based on sigma is enabled, flag */
            if (Data::flagsigma) {
             flag_data_sigma(output,Data::numRows,Data::sigXY[0],Data::sigXY[1],Data::sigXY[2],Data::sigXY[3],Nthreads);
            }
        /* w-snapshot, project w coords to best fitting plane */
        /* w= a*u+b*v+delta_w */
        double wpa,wpb; float wmax,wmin;
        if (Cmd::w_snapshot) {
         project_wplane(output,Data::numRows, Nthreads, &wpa, &wpb, &wmax, &wmin);
         cout<<"W projection plane: "<<wpa<<","<<wpb<<endl;
         cout<<"Snapshot W range (m): "<<wmin<<" : "<<wmax<<endl;
         amaxW=(fabs(wmin)>fabs(wmax)?fabs(wmin):fabs(wmax));
        } else {
          wpa=wpb=0.0;
          wmax=wmin=0.0f;
        }

/*************************************************/
    // W: kernel support in image, Nx*delta_l=image FOV
    // delta_u= 1/(Nx*delta_l)
    /*
       -1 |__________| 1 => Np image pixels  = (Nx*delta_l) range
     
      so, image conv. kernel pix width = (Nx delta_l) x 1/Np

      with zero padding, con. kernel support in image = (Nx delta_l) x 1/Np x Nz
      Nz: zero padded length
      
      After the FT, delta_u1 = 1/conv.pix_width= 1/(Nx delta_l) x Np x 1/Nz
      and uv grid pix width delta_u= 1/(Nx delta_l)
      so for 1 original uv grid pixel, we have Nz/Np pixels

      after scaling uv grid to [-0.5,0.5], 
      delta_u1= 1/(Nx delta_l) x Np/Nz x UVscale
              = 1/delta_l1 x 1/Nz x UVscale
    */


    float *wparr; /* w projection kernel W values */
    int Nw=Cmd::Stacks;
    int wgrid_step=Cmd::wgrid_step; /* 0: sqrt(), 1: log() :2 linear, 3: user, 4: w^(expW) */
    /* override this value for user supplied w step option -i 3 */
    if (wgrid_step==4) {
      /* read text file to determine no. of w planes */
      FILE *cfilep=fopen(Cmd::File_wstep,"r");
      fscanf(cfilep,"%d\n",&Nw);
      if ((wparr=(float*)calloc((size_t)Nw,sizeof(float)))==0) {
       fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
       exit(1);
      }
      for (int ci=0; ci<Nw; ci++) {
       fscanf(cfilep,"%f\n",&wparr[ci]);
      }
      fclose(cfilep);
    } else {
     if ((wparr=(float*)calloc((size_t)Nw,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
     }
     wparr[0]=0.0f; 
    }

    //if (Cmd::use_gpu) {  Obsolete ???
    //  wgrid_step=0; /* use sqrt() spacing for GPU gridding */
    //}
    float wp_delta=1.0f;
    if (wgrid_step==0) {
    /* grid in steps of sqrt(w) */ 
      wp_delta=sqrtf(amaxW)/((float)Nw-1.0f);
      Cmd::expW=0.5f;
    } else if (wgrid_step==1) {
     /* grid in steps of log(w) */
     wp_delta=logf(amaxW)/((float)Nw-1.0f);
    } else if (wgrid_step==2) {
     /* steps linear */
     wp_delta=(amaxW)/((float)Nw-1.0f);
     Cmd::expW=1.0f;
    } else if (wgrid_step==3) {
     wp_delta=powf(amaxW,Cmd::expW)/((float)Nw-1.0f);
    } else { 
     /* user supplied w steps, not needed */
    }
    if (wgrid_step<=3) {
     for (int ci=1; ci<Nw; ci++) {
      wparr[ci]=wparr[ci-1]+wp_delta;
     }
     if (wgrid_step==0) {
      /* take back the sqrt(w)^2 */
      for (int ci=1; ci<Nw; ci++) {
        wparr[ci]=wparr[ci]*wparr[ci];
      }
     } else if (wgrid_step==1) {
      /* take back the exp(log(w)) */
      for (int ci=1; ci<Nw; ci++) {
        wparr[ci]=expf(wparr[ci]);
      }
     } else if (wgrid_step==3) {
      /* take back the pow(w,1/expW) */
      for (int ci=1; ci<Nw; ci++) {
        wparr[ci]=powf(wparr[ci],1.0f/Cmd::expW);
      }
     } 
    }
    /* oversampling ratio for different kernels */
    double *oversamp;
    if ((oversamp=(double*)calloc((size_t)Nw,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    double dl1=(double)Cmd::ImRes*ARCSEC_TO_RADIANS;
    /* common term =sqrt(delta_l * Np * 1/pi * |(sqrt(1-(l^2+m^2))-1)|)*/
    /* max(l^2+m^2) = 2 * ((delta_l)*Npix)^2 */
    double comterm=2.0*(dl1*(double)Nx)*(dl1*(double)Nx);
    if (comterm>1.0) { 
     comterm=1.0; 
    } else {
     comterm=1.0-sqrt(1.0-comterm);
    }
    comterm *=dl1*(double)Np/M_PI;
    comterm = sqrt(comterm);
    for (int ci=0; ci<Nw; ci++) {
      oversamp[ci]=comterm*sqrt(wparr[ci]);
      /* check if this fits with GPU, else reduce oversample ratio */
      if ((double)Np*oversamp[ci]>(double)(1024-Np)) {
       oversamp[ci]=(ci>0?oversamp[ci-1]:0.0);;
      }
      //printf("%d k=%lf\n",ci,oversamp[ci]);
    }

    complex float *wkernel;
    if ((wkernel=(complex float*)calloc((size_t)Np*Np*Nw,sizeof(complex float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    int *wpsupportX,*wpsupportY;
    if ((wpsupportX=(int*)calloc((size_t)Nw,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((wpsupportY=(int*)calloc((size_t)Nw,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    float wtol=Cmd::conv_cutoff;

    /* calculate wkernel in 4 stages, depending on oversampling ratio
       [1,1.5) [1.5,2.5) [2.5,3.5) [3.5,...) */
    int *startWplane,*noWplanes;
    int Ovs=12;
    if ((startWplane=(int*)calloc((size_t)Ovs,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((noWplanes=(int*)calloc((size_t)Ovs,sizeof(int)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    int wx=0;
    for (int ci=0; ci<Ovs; ci++) {
      startWplane[ci]=wx;
      while(wx<Nw && oversamp[wx]<(double)ci+1.0) {
        wx++;
        noWplanes[ci]++;
      }
      /* last one should include all other remaining values */
      if (ci==Ovs-1 && wx<Nw) {
       noWplanes[ci]+=Nw-wx;
      }
    }
    //for(int ci=0;ci<Ovs;ci++) {
    // printf("<%lf %d %d\n",(double)ci+1.0,noWplanes[ci],startWplane[ci]);
    //}

    int Msupport;
    //generate_2d_wconv_kernels(M, Np, delta_l, wparr, Nw, Nz, wkernel,wtol,wpsupportX,wpsupportY,&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,Nthreads);
    /* iterate with different oversampling to generate wkernels */
    for (int ci=0; ci<Ovs; ci++) {
      if (noWplanes[ci]>0) {
       /* for very short kernels, lower the cutoff to sample the tail better ?? */
#ifdef HAVE_CUDA
        if (Cmd::use_gpu) {
         generate_2d_wconv_kernels_over_gpu(M, Np, delta_l, &wparr[startWplane[ci]], noWplanes[ci], Nz, &wkernel[Np*Np*startWplane[ci]],wtol,&wpsupportX[startWplane[ci]],&wpsupportY[startWplane[ci]],&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,(double)ci+1.0,Nthreads);
        } else {
         generate_2d_wconv_kernels_over(M, Np, delta_l, &wparr[startWplane[ci]], noWplanes[ci], Nz, &wkernel[Np*Np*startWplane[ci]],wtol,&wpsupportX[startWplane[ci]],&wpsupportY[startWplane[ci]],&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,(double)ci+1.0,Nthreads);
        }
#endif /* HAVE_CUDA */
#ifndef HAVE_CUDA
         generate_2d_wconv_kernels_over(M, Np, delta_l, &wparr[startWplane[ci]], noWplanes[ci], Nz, &wkernel[Np*Np*startWplane[ci]],wtol,&wpsupportX[startWplane[ci]],&wpsupportY[startWplane[ci]],&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,(double)ci+1.0,Nthreads);
#endif /* !HAVE_CUDA */
      }
    }

    free(oversamp);
    free(startWplane);
    free(noWplanes);

/* DEBUG WPROJ */
  //FILE *cfilep=fopen("debug_w_r","w+");
  //fprintf(cfilep,"%d %d\n",Np,Nw);
  //for (int ci=0; ci<Np*Np*Nw; ci++) {
  //  fprintf(cfilep,"%f %f\n",crealf(wkernel[ci]),cimagf(wkernel[ci]));
  //}
  //fclose(cfilep);


    /* array to lookup kernel values in [-1,1] add guard value at beginning
       and at end. so array is [-inf,-max_u,-.,..,0,...,max_u,+inf] so binary search works
     Np/2 is the zero pixel  */
    float *ugrid;
    if ((ugrid=(float*)calloc((size_t)Np+2,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }

    ugrid[0]=-1e6; /* guard values */
    ugrid[Np+1]=1e6;
    /* find zero pixel by finding zero wplane 1.0 location */
    int pix0=my_icamax(Np*Np,wkernel,1);
    pix0=pix0/Np+1;
    /* now PSWF is in [-1,1], scale this to conv. kernel delta_u 
       so we can directly search for right value */
    float scaled_uvdel=1.0f/delta_l*(1.0f/(float)(Nz))*Cmd::UVscale;
    ugrid[pix0]=0.0f; /* zero pixel */
    for (int ci=pix0+1; ci<=Np; ci++) {
      ugrid[ci]=ugrid[ci-1]+scaled_uvdel;
    }
    for (int ci=pix0-1; ci>=1; ci--) {
      ugrid[ci]=ugrid[ci+1]-scaled_uvdel;
    }

    /* regenerate wparray with guard values */
    float *wparr1;
    if ((wparr1=(float*)calloc((size_t)Nw,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    memcpy(wparr1,wparr,sizeof(float)*Nw);
    free(wparr);
    if ((wparr=(float*)calloc((size_t)Nw+2,sizeof(float)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    wparr[0]=-1e6; 
    wparr[Nw+1]=1e6; 
    memcpy(&wparr[1],wparr1,sizeof(float)*Nw);
    free(wparr1);

/*************************************************/
    float *weights=(float*)wgtin;
    /* second half of wgtin is used to store sum of lm transforemed image, when
      w-snapshots are used */
    complex float *psfgrid=dout;

    int weightmode=Cmd::weightmode; /* 0: uniform, 1: robust, 2:natural 3: PM */
    float rval=5.0f*powf(10.0f,-Cmd::robustval); /* robust parameter 5 x 10^-r */
    float sumweights=0.0f; /* add up the sum of weights */
    float sumweights2=0.0f; /* add up the sum of weights^2 (for Parseval) */
    if (weightmode==WMODE_PIPEMENON) { 
     if (!Cmd::use_gpu) {
      weightuvdata_pipe_menon(output,Data::numRows, Cmd::UVscale, deltaU, Data::waveLen, weights,  Nx, Ny, Nthreads, Cmd::pm_convmode, Cmd::pm_maxiter, &sumweights, &sumweights2, Cmd::pm_imscale*M_PI, wkernel, ugrid, Np, wpsupportX[0]); /* imscale scaled by PI for the FT scaling */
     } else {
      /* disable GPU version -- not working properly */
      //weightuvdata_pipe_menon_gpu(output,Data::numRows, Cmd::UVscale, deltaU, Data::waveLen, weights,  Nx, Ny, Nthreads, Cmd::pm_convmode, Cmd::pm_maxiter, &sumweights, &sumweights2, Cmd::pm_imscale*M_PI, wkernel, ugrid, Np, wpsupportX[0]); /* imscale scaled by PI for the FT scaling */
      weightuvdata_pipe_menon(output,Data::numRows, Cmd::UVscale, deltaU, Data::waveLen, weights,  Nx, Ny, Nthreads, Cmd::pm_convmode, Cmd::pm_maxiter, &sumweights, &sumweights2, Cmd::pm_imscale*M_PI, wkernel, ugrid, Np, wpsupportX[0]); /* imscale scaled by PI for the FT scaling */
     }
     sumweights_avg+=sumweights;
     sumweights2_avg+=sumweights2;
    } else {
     /* calculate the weights for the (ungridded) uv data */
     weightuvdata(output,Data::numRows, Cmd::UVscale, Data::waveLen, rval*rval, weights,  Nx, Ny, Nthreads, (weightmode==WMODE_PMTAPER?WMODE_UNIFORM:weightmode), &sumweights, &sumweights2);
     sumweights_avg+=sumweights;
     sumweights2_avg+=sumweights2;
    }
/************* print out weights ******************/
/*    FILE *cfilep=fopen("debug_weights","w+");
    for (unsigned int nr=0; nr<Data::numRows; nr++) {
      if (!output[nr].flag) {
       fprintf(cfilep,"%f %f %e\n",output[nr].u,output[nr].v,output[nr].wt);
      }
    }
    fclose(cfilep);
*/
/************* end print out weights ******************/
    /* recalculate uaxis */
    pix0=(int)roundf(0.5f*(float)Np-(float)Nz/(float)Np);
    ugrid[pix0]=0.0f; /* zero pixel */
    for (int ci=pix0+1; ci<=Np; ci++) {
      ugrid[ci]=ugrid[ci-1]+scaled_uvdel;
    }
    for (int ci=pix0-1; ci>=1; ci--) {
      ugrid[ci]=ugrid[ci+1]-scaled_uvdel;
    }
#ifdef HAVE_CUDA
    if (Cmd::use_gpu) {
     cuda_griddata(output,Data::numRows, Cmd::UVscale, Data::waveLen, scaled_uvdel, din, psfgrid, Nx, Ny, Cmd::expW, wparr, Nw, wpsupportX, wpsupportY, Msupport, (complex float*)wkernel, Np, Nz, Cmd::gpu_threads,MAX_GPU,Cmd::imgmode);
    } else {
     griddata(output,Data::numRows, Cmd::UVscale, Data::waveLen, deltaU, din,  psfgrid, Nx, Ny, wparr, Nw, wpsupportX, wpsupportY, wkernel, ugrid, Np, Nthreads,Cmd::imgmode);
    }
#endif
#ifndef HAVE_CUDA
    /* image grid is din */
    griddata(output,Data::numRows, Cmd::UVscale, Data::waveLen, deltaU, din,  psfgrid, Nx, Ny, wparr, Nw, wpsupportX, wpsupportY, wkernel, ugrid, Np, Nthreads,Cmd::imgmode);
#endif

    free(ugrid);
    free(wparr);
    free(wkernel);
    free(wpsupportX);
    free(wpsupportY);
   if (weightmode==WMODE_PMTAPER) {
    /* apply any tapering to gridded psfgrid and gridded data (din) */
    tapergrid(din,psfgrid,Nx,Ny,Nx/2-1,deltaU/Cmd::UVscale,Cmd::pm_convmode,Nthreads);
   }

   if (Cmd::imgmode==IMG_I||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
    weightdata(din,psfgrid,weights,Nx,Ny,rval,Nthreads);
   /* save weights */
   /* FFTshift weights */
   woffset=Nx*Ny;
   do_weight_fftshift(weights,&weights[woffset],Nx,Ny,Cmd::Nx,Npad);
   char wbuff[4096];
   if( snapshotmode ) {
    if (Cmd::file_qualifier) {
     if (Cmd::AvgChannels==2) { 
      sprintf(wbuff, "!%s_%s_W_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
     } else {
      sprintf(wbuff, "!%s_%s_W_%d.fits", Base,Cmd::file_qualifier,epoch); 
     }
    } else { 
     if (Cmd::AvgChannels==2) { 
      sprintf(wbuff, "!%s_W_%d_%d.fits", Base,chan,epoch); 
     } else {
      sprintf(wbuff, "!%s_W_%d.fits", Base,epoch); 
     }
    }
   } else {
    if (Cmd::file_qualifier) {
     if (Cmd::AvgChannels==2) { 
      sprintf(wbuff, "!%s_%s_W_%d.fits", Base,Cmd::file_qualifier,chan); 
     } else {
      sprintf(wbuff, "!%s_%s_W.fits", Base,Cmd::file_qualifier); 
     }
    } else { 
     if (Cmd::AvgChannels==2) { 
      sprintf(wbuff, "!%s_W_%d.fits", Base,chan); 
     } else {
      sprintf(wbuff, "!%s_W.fits", Base); 
     }
    }
   }
   writeFITSW(wbuff, &weights[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1.0f, 1);
   }


/*************************************************/
    // save conv. kernel as a FITS image
    if (Cmd::convkernel) {
    if(!convkernel_created) {
     float *convBuff;
     if ((convBuff=(float*)calloc((size_t)Cmd::Nx*Cmd::Nx,sizeof(float)))==0) {
      fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
      exit(1);
     }
     sprintf(buff,"!%s_convkernel.fits", Base); 
     create_conv_kernel(convBuff,Cmd::Nx,M);
     writeFITSIMG(buff, convBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1, 1.0, 1.0, 0.0, 1.0);
     free(convBuff);
     convkernel_created=true;
    }
    }
/*************************************************/

    float apodization_tol=Cmd::apod_cutoff; /* cutoff for apodization correction */
    /* FFT dout => wgtin */
    fftwf_execute(fftp0);
    /* find the peak value of the 4 corners of the image to normalize PSF */
    float peakpsf=cabsf(wgtin[0]);
    /* catch if PSF==0 */
    if (peakpsf==0.0f) {
     peakpsf=1.0f;
    }
    float invpsf=1.0f/peakpsf;
    /* scale factor to get back right values for gridded weights,
     when writing set BSCALE to this and then set BSCALE to 1 
     */
    /* grid should be downscaled by the con. kernel width, 
     only when non-adaptive weights are used (because of grid pixels) */
    float gridscale=invpsf;
    if (weightmode!=WMODE_PIPEMENON) {
     gridscale *=(float)M; 
    }
    cout<<"Peak PSF ="<<peakpsf<<" Grid Scale="<<gridscale<<endl;
/********* fit Gaussian PSF ************************/
    /* maxUV/(1/pixel_size) gives approx PSF half width in pixels
       use 5 times this value */
    int P=(int)ceilf(5.0f*Data::maxUV[2]/L);
    double *Pval,*Pgridl,*Pgridm;
    if ((Pval=(double*)calloc((size_t)2*P*2*P,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((Pgridl=(double*)calloc((size_t)2*P*2*P,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    if ((Pgridm=(double*)calloc((size_t)2*P*2*P,sizeof(double)))==0) {
     fprintf(stderr,"%s: %d: no free memory\n",__FILE__,__LINE__);
     exit(1);
    }
    /* read pixel by pixel (need to cast to double) */
    /* 0,0 in grid is at pixel P,P */
    /* top left */
    for (int cq=0; cq<P; cq++) {
     for (int cp=0; cp<P; cp++) {
      Pgridl[P*2*(P+cq)+P+cp]=(double)(cp);
      Pgridm[P*2*(P+cq)+P+cp]=(double)(cq);
      Pval[P*2*(P+cq)+P+cp]=(double)crealf(wgtin[cp+cq*Nx]);
     } 
    }
    /* top right */
    for (int cq=0; cq<P; cq++) {
     for (int cp=0; cp<P; cp++) {
      Pgridl[P*2*(P+cq)+cp]=(double)(cp-P);
      Pgridm[P*2*(P+cq)+cp]=(double)(cq);
      Pval[P*2*(P+cq)+cp]=(double)crealf(wgtin[Nx-P+cp+cq*Nx]);
     } 
    }
    /* bottom left*/
    for (int cq=0; cq<P; cq++) {
     for (int cp=0; cp<P; cp++) {
      Pgridl[P*2*cq+P+cp]=(double)(cp);
      Pgridm[P*2*cq+P+cp]=(double)(cq-P);
      Pval[P*2*(cq)+P+cp]=(double)crealf(wgtin[cp+(Ny-P+cq)*Nx]);
     } 
    }
    /* bottom right */
    for (int cq=0; cq<P; cq++) {
     for (int cp=0; cp<P; cp++) {
      Pgridl[P*2*cq+cp]=(double)(cp-P);
      Pgridm[P*2*cq+cp]=(double)(cq-P);
      Pval[P*2*(cq)+cp]=(double)crealf(wgtin[Nx-P+cp+(Ny-P+cq)*Nx]);
     } 
    }
    my_dscal(4*P*P,1.0/Pval[P*2*P+P],Pval);
    fit_gaussian_psf(Pval,Pgridl,Pgridm,4*P*P,&bmaj,&bmin,&bpa);
    /* convert bmaj,bmin,bpa to degrees */
    bmaj *=(double)(2.0f*Cmd::ImRes/3600.0f); 
    bmin *=(double)(2.0f*Cmd::ImRes/3600.0f); 
    double t1,t2;
    sincos(bpa,&t1,&t2);
    bpa *=-180.0/M_PI;
    avbpax+=t2;
    avbpay+=t1;
    /* if snapshot mode, also store average value of PSF */
    avbmaj+=bmaj;
    avbmin+=bmin;
    free(Pval); free(Pgridl); free(Pgridm);
/***************************************************/
    float *imgBuff = (float*)dout; 
   if (Cmd::imgmode==IMG_I||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
    /* SHIFT wgtin => dout */
    do_data_fftshift(wgtin,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_PSF_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_PSF_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_PSF_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_PSF_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_PSF_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_PSF.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_PSF_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_PSF.fits", Base); 
       }
     }
    }
    writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);

   /*********** also write imaginary part of PSF **********/
   if(Cmd::create_imaginary) {
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_PSFi_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_PSFi_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_PSFi_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_PSFi_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_PSFi_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_PSFi.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_PSFi_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_PSFi.fits", Base); 
      }
     }
    }

    writeFITSIMG(buff, &imgBuff[Nx*Ny], Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);
   }
   /*********** done writing imaginary part of PSF ********/

    /* write Gridded data real+imag as separate FITS files */
    woffset=Nx*Ny;
    do_grid_fftshift(din,(float*)dout,invpsf,&weights[woffset],Nx,Ny,Cmd::Nx,Npad);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GR_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GR_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GR_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GR_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GR_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GR.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GR_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GR.fits", Base); 
       }
     }
    }
    writeFITSW(buff,imgBuff, Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 1);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GI_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GI_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GI_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GI_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GI_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GI.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GI_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GI.fits", Base); 
       }
     }
    }
    woffset=Cmd::Nx*Cmd::Ny;
    writeFITSW(buff,&imgBuff[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 1);
   }


/*************************************************/
    /* FFT din => dout */
    fftwf_execute(fftp1);
    /* dout is the image, write back to din, real and imaginary parts separate */
    /* SHIFT dout => din */
    imgBuff = (float*)din; 
    do_data_fftshift(dout,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
    /* if w-snapshot mode, lm transform from din => wgtin (because wgtin is not used) */
    double delta_ll=(double)Cmd::ImRes*M_PI/648000.0;
    if(Cmd::w_snapshot) {
     do_image_lm_transform(imgBuff,(float*)wgtin,Cmd::Nx,Cmd::Nx,wpa,wpb,delta_ll,ref[0],ref[1],Nthreads);
     imgBuff=(float*)wgtin;
    }
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_I_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_I_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_I_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_I_%d.fits", Base,epoch); 
       }
     }
     /* also add current snapshot to average image, only I image */
     my_faxpy(Cmd::Nx*Cmd::Ny,imgBuff,1,1.0f,avgImage,1);
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_I_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_I.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_I_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_I.fits", Base); 
       }
     }
    }

    writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);

   /*********** also write imaginary part of image **********/
   if(Cmd::create_imaginary) {
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Ii_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_Ii_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Ii_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_Ii_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Ii_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_Ii.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Ii_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_Ii.fits", Base); 
       }
     }
    }

    writeFITSIMG(buff, &imgBuff[Nx*Ny], Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);
   }
   /*********** done writing imaginary part of image ********/

   if (Cmd::imgmode==IMG_IQUV0||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
    /******************* Begin Q Grid *********************************/
    if (Cmd::imgmode==IMG_IQUVF) {
    imgBuff = (float*)dout; 
    woffset=Nx*Ny;
    if (weightmode==WMODE_PMTAPER) {
     /* apply any tapering to gridded data (din) */
     tapergrid(&din[woffset],0,Nx,Ny,Nx/2-1,deltaU/Cmd::UVscale,Cmd::pm_convmode,Nthreads);
    }

    do_grid_fftshift(&din[woffset],(float*)dout,invpsf,&weights[Nx*Ny],Nx,Ny,Cmd::Nx,Npad);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRQ_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GRQ_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRQ_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GRQ_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRQ_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GRQ.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRQ_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GRQ.fits", Base); 
       }
     }
    }
    writeFITSW(buff,imgBuff, Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 2);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIQ_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GIQ_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIQ_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GIQ_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIQ_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GIQ.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIQ_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GIQ.fits", Base); 
       }
     }
    }
    woffset=Cmd::Nx*Cmd::Ny;
    writeFITSW(buff,&imgBuff[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 2);
    }

    /******************* End Q Grid *********************************/
    imgBuff = (float*)din; 
    fftwf_execute(fftp2);
    do_data_fftshift(dout,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
    if(Cmd::w_snapshot) {
     do_image_lm_transform(imgBuff,(float*)wgtin,Cmd::Nx,Cmd::Nx,wpa,wpb,delta_ll,ref[0],ref[1],Nthreads);
     imgBuff=(float*)wgtin;
    }
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Q_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_Q_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Q_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_Q_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Q_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_Q.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Q_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_Q.fits", Base); 
       }
     }
    }
    writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 2,bmaj,bmin,bpa,sumweights);

    if(Cmd::create_imaginary) {
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Qi_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_Qi_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Qi_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_Qi_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Qi_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_Qi.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Qi_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_Qi.fits", Base); 
       }
     }
    }

    writeFITSIMG(buff, &imgBuff[Nx*Ny], Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 2,bmaj,bmin,bpa,sumweights);
   }


   /******************* Begin U Grid *********************************/
    if (Cmd::imgmode==IMG_IQUVF) {
    imgBuff = (float*)dout; 
    woffset=2*Nx*Ny;
    if (weightmode==WMODE_PMTAPER) {
     /* apply any tapering to gridded data (din) */
     tapergrid(&din[woffset],0,Nx,Ny,Nx/2-1,deltaU/Cmd::UVscale,Cmd::pm_convmode,Nthreads);
    }

    do_grid_fftshift(&din[woffset],(float*)dout,invpsf,&weights[Nx*Ny],Nx,Ny,Cmd::Nx,Npad);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRU_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GRU_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRU_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GRU_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRU_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GRU.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRU_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GRU.fits", Base); 
       }
     }
    }
    writeFITSW(buff,imgBuff, Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 3);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIU_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GIU_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIU_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GIU_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIU_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GIU.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIU_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GIU.fits", Base); 
       }
     }
    }
    woffset=Cmd::Nx*Cmd::Ny;
    writeFITSW(buff,&imgBuff[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 3);
    }
    /******************* End U Grid *********************************/

    imgBuff = (float*)din; 
    fftwf_execute(fftp3);
    do_data_fftshift(dout,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
    if(Cmd::w_snapshot) {
     do_image_lm_transform(imgBuff,(float*)wgtin,Cmd::Nx,Cmd::Nx,wpa,wpb,delta_ll,ref[0],ref[1],Nthreads);
     imgBuff=(float*)wgtin;
    }
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_U_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_U_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_U_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_U_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_U_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_U.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_U_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_U.fits", Base); 
       }
     }
    }
    writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 3,bmaj,bmin,bpa,sumweights);

    if(Cmd::create_imaginary) {
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Ui_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_Ui_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Ui_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_Ui_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Ui_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_Ui.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Ui_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_Ui.fits", Base); 
       }
     }
    }

    writeFITSIMG(buff, &imgBuff[Nx*Ny], Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 3,bmaj,bmin,bpa,sumweights);
   }

   /******************* Begin V Grid *********************************/
    if (Cmd::imgmode==IMG_IQUVF) {
    imgBuff = (float*)dout; 
    woffset=3*Nx*Ny;
    if (weightmode==WMODE_PMTAPER) {
     /* apply any tapering to gridded data (din) */
     tapergrid(&din[woffset],0,Nx,Ny,Nx/2-1,deltaU/Cmd::UVscale,Cmd::pm_convmode,Nthreads);
    }

    do_grid_fftshift(&din[woffset],(float*)dout,invpsf,&weights[Nx*Ny],Nx,Ny,Cmd::Nx,Npad);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRV_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GRV_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRV_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GRV_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GRV_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GRV.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GRV_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GRV.fits", Base); 
       }
     }
    }
    writeFITSW(buff,imgBuff, Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, Data::chanFreq, Data::chanBW, gridscale, 4);
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIV_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_GIV_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIV_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_GIV_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_GIV_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_GIV.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_GIV_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_GIV.fits", Base); 
       }
     }
    }
    woffset=Cmd::Nx*Cmd::Ny;
    writeFITSW(buff,&imgBuff[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir,  (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), gridscale, 4);
    }
    /******************* End V Grid *********************************/


    imgBuff = (float*)din; 
    fftwf_execute(fftp4);
    do_data_fftshift(dout,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
    if(Cmd::w_snapshot) {
     do_image_lm_transform(imgBuff,(float*)wgtin,Cmd::Nx,Cmd::Nx,wpa,wpb,delta_ll,ref[0],ref[1],Nthreads);
     imgBuff=(float*)wgtin;
    }
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_V_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_V_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_V_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_V_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_V_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_V.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_V_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_V.fits", Base); 
       }
     }
    }
    writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir,  (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 4,bmaj,bmin,bpa,sumweights);

    if(Cmd::create_imaginary) {
    if(snapshotmode) {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Vi_%d_%d.fits", Base,Cmd::file_qualifier,chan,epoch); 
       } else {
        sprintf(buff, "!%s_%s_Vi_%d.fits", Base,Cmd::file_qualifier,epoch); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Vi_%d_%d.fits", Base,chan,epoch); 
       } else {
        sprintf(buff, "!%s_Vi_%d.fits", Base,epoch); 
       }
     }
    } else {
     if (Cmd::file_qualifier) {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_%s_Vi_%d.fits", Base,Cmd::file_qualifier,chan); 
       } else {
        sprintf(buff, "!%s_%s_Vi.fits", Base,Cmd::file_qualifier); 
       }
     } else {
       if (Cmd::AvgChannels==2) {
        sprintf(buff, "!%s_Vi_%d.fits", Base,chan); 
       } else {
        sprintf(buff, "!%s_Vi.fits", Base); 
       }
     }
    }

    writeFITSIMG(buff, &imgBuff[Nx*Ny], Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir,  (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 4,bmaj,bmin,bpa,sumweights);
   }


   }
/*************************************************/
/*************************************************/
        /* reset binary files to zero, if reusing them */
        if ( Cmd::AvgChannels==2 ) {
         reset_binary_file(din, dout, wgtin, Nx, Ny,Cmd::imgmode);
        }

     }
     freeData();//free the buffer containing data from table.
     epoch++;
     /* reset binary files to zero, if reusing them */
     if ( snapshotmode ) {
         reset_binary_file(din, dout, wgtin, Nx, Ny,Cmd::imgmode);
     }
    }
    /************** end snapshot iterations *********************/
    /* also close the file descriptors, delete the files */
    do_destroy_fftw_plan(fftp0);
    do_destroy_fftw_plan(fftp1);
    if(Cmd::imgmode==IMG_IQUV0||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
     do_destroy_fftw_plan(fftp2);
     do_destroy_fftw_plan(fftp3);
     do_destroy_fftw_plan(fftp4);
    } 
    /* write average image, if snapshots are made */
    if( snapshotmode ) {
     if (Cmd::file_qualifier) {
      sprintf(buff, "!%s_%s_I.fits", Base,Cmd::file_qualifier); 
     } else {
      sprintf(buff, "!%s_I.fits", Base); 
     }
     double avscale=1.0/(double)epoch;
     my_fscal(Cmd::Nx*Cmd::Ny,(float)avscale,avgImage); 
     /* calculate average PSF */
     bmaj=avbmaj*avscale;
     bmin=avbmin*avscale;
     bpa=atan2(avbpay*avscale,avbpax*avscale);
     writeFITSIMG(buff, avgImage, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, Data::chanFreq, Data::chanBW, 1,bmaj,bmin,-bpa*180.0/M_PI,sumweights_avg);
    }

    close_binary_file(uvgridname,uvgridid,din,imgridname,imgridid,dout,wtgridname, wtgridid, wgtin, Nx,Ny, Cmd::Nx, Cmd::Ny, Cmd::imgmode, snapshotmode);
    stop_fftw();

 
   free(TableName1);
   cout<<"Done."<<endl;    
}
