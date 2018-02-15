#include "data.h"
#include "helper.h"
#include "pthgridder.h"
#include <fstream>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <pthread.h>

#include <largefft.h>

#include "Driver.h"

JNIEXPORT jstring JNICALL Java_Driver_stringMethod
  (JNIEnv *env, jobject obj, jstring string) {
    const char *tablename=env->GetStringUTFChars(string,0);

    Table t(tablename);

    env->ReleaseStringUTFChars(string,tablename);

    return string;
}

JNIEXPORT jint JNICALL Java_Driver_readMSAndBack
  (JNIEnv *env, jobject obj, jobjectArray argc) {

  int stringCount = env->GetArrayLength(argc);
  /* alloc memory for char* array */
  char **myArgs;
  if ((myArgs=(char **)malloc((size_t)(stringCount+2)*sizeof(char*)))==0) {
    fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
    exit(1);
  }
  /* first two locations are not used in getopt() */
  if ((myArgs[0]=(char *)calloc((size_t)(1),sizeof(char)))==0) {
         fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
         exit(1);
  }
  if ((myArgs[1]=(char *)calloc((size_t)(1),sizeof(char)))==0) {
         fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
         exit(1);
  }


  for (int ci=0; ci<stringCount; ci++) {
        jstring javastring = (jstring) (env->GetObjectArrayElement(argc, ci));
        const char *rawString = env->GetStringUTFChars(javastring, 0);
        /* alloc memory to store this string, set to zero */
        if ((myArgs[ci+2]=(char *)calloc((size_t)(strlen(rawString)+1),sizeof(char)))==0) {
         fprintf(stderr,"%s: %d: No free memory\n",__FILE__,__LINE__);
         exit(1);
        }
        strcpy(myArgs[ci+2],rawString);
        // Don't forget to call `ReleaseStringUTFChars` when you're done.
        env->ReleaseStringUTFChars(javastring,rawString);
  }

  Cmd::ParseCmdLine(stringCount+2,myArgs);

  float ref[2]={0,0}; //reference angle for this table
  int epoch = 0;
  bool done = false;
  bool convkernel_created=false;

  float uvLim,L,deltaU;

  // FIXME: add hdfs support for reading tables!!!
  // left here
  cout << "  ** MS name ** : " << Cmd::TableName << endl;

  /* now make an image */
  //Create a table with the ms so we can directly image it if this is not epoch imaging.
  Table t(Cmd::TableName);
  //stuff required for epoch imaging.
  MeasurementSet ms(Cmd::TableName);
  //backup MS name
  char *TableName1=strdup(Cmd::TableName);
  char *Base=basename(TableName1);

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

    while(!done){
           //we are doing normal imaging so only one iteration.
            done = true;
            //use the normal table instead of Ms
            // only sort table is image is very big >12000 pix

            // FIXME: add hdfs support here
            loadData(t, Cmd::TableName, Cmd::AvgChannels, (Cmd::create_imaginary>0), false);

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
    int Msupport;
    //generate_2d_wconv_kernels(M, Np, delta_l, wparr, Nw, Nz, wkernel,wtol,wpsupportX,wpsupportY,&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,Nthreads);
    /* iterate with different oversampling to generate wkernels */
    for (int ci=0; ci<Ovs; ci++) {
      if (noWplanes[ci]>0) {
 generate_2d_wconv_kernels_over(M, Np, delta_l, &wparr[startWplane[ci]], noWplanes[ci], Nz, &wkernel[Np*Np*startWplane[ci]],wtol,&wpsupportX[startWplane[ci]],&wpsupportY[startWplane[ci]],&Msupport,ref[0],ref[1],wpa,wpb,Cmd::w_snapshot,(double)ci+1.0,Nthreads);
     }
    }

    free(oversamp);
    free(startWplane);
    free(noWplanes);

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


    float *weights=(float*)wgtin;
    /* second half of wgtin is used to store sum of lm transforemed image, when
      w-snapshots are used */
    complex float *psfgrid=dout;

    int weightmode=Cmd::weightmode; /* 0: uniform, 1: robust, 2:natural 3: PM */
    float rval=5.0f*powf(10.0f,-Cmd::robustval); /* robust parameter 5 x 10^-r */
    float sumweights=0.0f; /* add up the sum of weights */
    float sumweights2=0.0f; /* add up the sum of weights^2 (for Parseval) */
    if (weightmode==WMODE_PIPEMENON) {
weightuvdata_pipe_menon(output,Data::numRows, Cmd::UVscale, deltaU, Data::waveLen, weights,  Nx, Ny, Nthreads, Cmd::pm_convmode, Cmd::pm_maxiter, &sumweights, &sumweights2, Cmd::pm_imscale*M_PI, wkernel, ugrid, Np, wpsupportX[0]);
     sumweights_avg+=sumweights;
     sumweights2_avg+=sumweights2;
    } else {
     /* calculate the weights for the (ungridded) uv data */
     weightuvdata(output,Data::numRows, Cmd::UVscale, Data::waveLen, rval*rval, weights,  Nx, Ny, Nthreads, (weightmode==WMODE_PMTAPER?WMODE_UNIFORM:weightmode), &sumweights, &sumweights2);
     sumweights_avg+=sumweights;
     sumweights2_avg+=sumweights2;
    }


    /* recalculate uaxis */
    pix0=(int)roundf(0.5f*(float)Np-(float)Nz/(float)Np);
    ugrid[pix0]=0.0f; /* zero pixel */
    for (int ci=pix0+1; ci<=Np; ci++) {
      ugrid[ci]=ugrid[ci-1]+scaled_uvdel;
    }
    for (int ci=pix0-1; ci>=1; ci--) {
      ugrid[ci]=ugrid[ci+1]-scaled_uvdel;
    }

  /* image grid is din */
    griddata(output,Data::numRows, Cmd::UVscale, Data::waveLen, deltaU, din,  psfgrid, Nx, Ny, wparr, Nw, wpsupportX, wpsupportY, wkernel, ugrid, Np, Nthreads,Cmd::imgmode);

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
   writeFITSW(wbuff, &weights[woffset], Cmd::Nx, Cmd::Nx, (double)deltaU/Cmd::UVscale*3600.0f, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1.0f, 1);
   }


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

    float *imgBuff = (float*)dout;
    /* SHIFT wgtin => dout */
    do_data_fftshift(wgtin,imgBuff,invpsf,Nx,Ny,M,Cmd::Nx,Npad,apodization_tol);
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
     writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);


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
 writeFITSIMG(buff, imgBuff, Cmd::Nx, Cmd::Nx, (double)Cmd::ImRes, Data::refDir, (Cmd::AvgChannels==2? Data::chanFreqs[chan] : Data::chanFreq), (Cmd::AvgChannels==2? Data::chanBWs[chan] : Data::chanBW), 1,bmaj,bmin,bpa,sumweights);

      /* reset binary files to zero, if reusing them */
        if ( Cmd::AvgChannels==2 ) {
         reset_binary_file(din, dout, wgtin, Nx, Ny,Cmd::imgmode);
        }

   }

     freeData();//free the buffer containing data from table.
     epoch++;

   }

    /* also close the file descriptors, delete the files */
    do_destroy_fftw_plan(fftp0);
    do_destroy_fftw_plan(fftp1);
    if(Cmd::imgmode==IMG_IQUV0||Cmd::imgmode==IMG_IQUV||Cmd::imgmode==IMG_IQUVF) {
     do_destroy_fftw_plan(fftp2);
     do_destroy_fftw_plan(fftp3);
     do_destroy_fftw_plan(fftp4);
    }
   close_binary_file(uvgridname,uvgridid,din,imgridname,imgridid,dout,wtgridname, wtgridid, wgtin, Nx,Ny, Cmd::Nx, Cmd::Ny, Cmd::imgmode, snapshotmode);
    stop_fftw();


   free(TableName1);

  /* free up all memory */
  for (int ci=0; ci<stringCount+2; ci++) {
    free(myArgs[ci]);
  }
  free(myArgs);


  return 0;
}
