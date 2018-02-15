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

#include "helper.h"
#include "data.h"
#include <gridder.h>
#include <fitsio.h>
#include <fitsio2.h>
#include <largefft.h>
using namespace casa;

void 
writeFITSW(char *fname, float *ptr, unsigned int w, unsigned int h, double delta, double *ref, double freq, double deltaf, float bscale, int stokes)
{
    cout << "Writing FITS : " << fname << endl;
    long dims[4] = {w, h, 1, 1};
    int status = 0;
    fitsfile *fptr;
    double ref_ptr[2] = {0.0, 0.0};

    //create the file for fits data.
    fits_create_file(&fptr, fname, &status);
    //write the image header in fits file
    fits_create_img(fptr, FLOAT_IMG, 4, dims, &status);
    //add the BUNIT the keyword and the value JY/BEAM
    char bunit_val[] = "JY";
    fits_write_key(fptr, TSTRING, "BUNIT", bunit_val, "", &status);

    //add the BSCALE value
    fits_write_key(fptr, TFLOAT, "BSCALE", &bscale, "scale to get back data", &status);         
    //add the BZERO value
    float bzero = 0;
    fits_write_key(fptr, TFLOAT, "BZERO", &bzero, "", &status);         
    //add btype value        
    char btype[] = "Intensity";
    fits_write_key(fptr, TSTRING, "BTYPE", btype, "", &status);

    float pos[] = {
    2.000000000000E+03,
    1.800000000000E+02,
    4.821741666667E+01
    };
    
    fits_write_key(fptr, TFLOAT, "EQUINOX", pos, "", &status);
    fits_write_key(fptr, TFLOAT, "LONPOLE", pos+1, "", &status);
    fits_write_key(fptr, TFLOAT, "LATPOLE", pos+2, "", &status);

    //add the Ctype1
    char ctype1[] = "U---WAV";
    fits_write_key(fptr, TSTRING, "CTYPE1", ctype1, "", &status);
    double crpix1 =(double) (w/2+1);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", ref_ptr, "", &status); 
    double dl = delta/3600.0;
    fits_write_key(fptr, TDOUBLE, "CDELT1", &dl, "", &status);
    char cunit[] = "lambda";
    fits_write_key(fptr, TSTRING, "CUNIT1", cunit, "", &status);
    //add the Ctype2         
    char ctype2[] = "V---WAV";
    fits_write_key(fptr, TSTRING, "CTYPE2", ctype2, "", &status);
    crpix1 = (double)(h/2+1); //crpix2 using the same var.        
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", ref_ptr+1, "", &status); 
    double dm =  dl;        
    fits_write_key(fptr, TDOUBLE, "CDELT2", &dm, "", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", cunit, "", &status);

    char ctype3[] = "FREQ";
    fits_write_key(fptr, TSTRING, "CTYPE3", ctype3, "", &status);
    crpix1 = 1; 
    fits_write_key(fptr, TDOUBLE, "CRPIX3", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", &freq, "", &status); 
    fits_write_key(fptr, TDOUBLE, "CDELT3", &deltaf, "", &status);
    char cunit3[]="Hz";
    fits_write_key(fptr, TSTRING, "CUNIT3", cunit3, "", &status);

    char ctype4[] = "STOKES";
    fits_write_key(fptr, TSTRING, "CTYPE4", ctype4, "", &status);
    crpix1 = (double)stokes; 
    double one=1.0;
    fits_write_key(fptr, TDOUBLE, "CRPIX4", &one, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL4", &crpix1, "", &status); 
    fits_write_key(fptr, TDOUBLE, "CDELT4", &one, "", &status);
    char cunit4[]="";
    fits_write_key(fptr, TSTRING, "CUNIT4", cunit4, "", &status);


    // auxilliary  info
    char origin[] = "ExCon (C) Sarod Yatawatta";
    fits_write_key(fptr, TSTRING, "ORIGIN", origin, "", &status);

    //write the pixel data to the file
    long fpixel[4] = {1, 1, 1, 1};
    long nelements = w * h;
    fits_write_pix(fptr, TFLOAT, fpixel, nelements, ptr, &status);

    //now reset BSCALE to 1
    float b0=1.0f;
    fits_update_key(fptr,TFLOAT,"BSCALE",&b0,0,&status);

    //close the file
    fits_close_file(fptr, &status);

    if(status)
        fits_report_error(stderr, status);
}


void 
writeFITSIMG(char *fname, float *ptr, unsigned int w, unsigned int h, double delta, double *ref, double freq, double deltaf, int stokes, double bmaj, double bmin, double bpa,float sumweight)
{
    cout << "Writing FITS : " << fname << endl;
    long dims[4] = {w, h, 1, 1};
    int status = 0;
    fitsfile *fptr;
    double ref_ptr[2] = {ref[0] * 180/M_PI, ref[1] * 180/M_PI};

    //create the file for fits data.
    fits_create_file(&fptr, fname, &status);
    //write the image header in fits file
    fits_create_img(fptr, FLOAT_IMG, 4, dims, &status);
    //add the BUNIT the keyword and the value JY/BEAM
    char bunit_val[] = "JY/BEAM";
    fits_write_key(fptr, TSTRING, "BUNIT", bunit_val, "", &status);

    //add the BSCALE value
    float bscale = 1.0;
    fits_write_key(fptr, TFLOAT, "BSCALE", &bscale, "", &status);         
    //add the BZERO value
    float bzero = 0;
    fits_write_key(fptr, TFLOAT, "BZERO", &bzero, "", &status);         
    //add btype value        
    char btype[] = "Intensity";
    fits_write_key(fptr, TSTRING, "BTYPE", btype, "", &status);

    /* write PSF info */
    float bmajf=(float)bmaj;
    fits_write_key(fptr, TFLOAT, "BMAJ", &bmajf, "", &status); 
    float bminf=(float)bmin;
    fits_write_key(fptr, TFLOAT, "BMIN", &bminf, "", &status); 
    float bpaf=(float)bpa;
    fits_write_key(fptr, TFLOAT, "BPA", &bpaf, "", &status); 

    float pos[] = {
    2.000000000000E+03,
    1.800000000000E+02,
    4.821741666667E+01
    };
    
    fits_write_key(fptr, TFLOAT, "EQUINOX", pos, "", &status);
    fits_write_key(fptr, TFLOAT, "LONPOLE", pos+1, "", &status);
    fits_write_key(fptr, TFLOAT, "LATPOLE", pos+2, "", &status);

    //add the Ctype1
    char ctype1[] = "RA---SIN";
    fits_write_key(fptr, TSTRING, "CTYPE1", ctype1, "", &status);
    double crpix1 =(double) (w/2+1);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL1", ref_ptr, "", &status); 
    double dl = -delta/3600.0;
    fits_write_key(fptr, TDOUBLE, "CDELT1", &dl, "", &status);
    char cunit[] = "deg";
    fits_write_key(fptr, TSTRING, "CUNIT1", cunit, "", &status);
    //add the Ctype2         
    char ctype2[] = "DEC--SIN";
    fits_write_key(fptr, TSTRING, "CTYPE2", ctype2, "", &status);
    crpix1 = (double)(h/2+1); //crpix2 using the same var.        
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", ref_ptr+1, "", &status); 
    double dm =  -dl;        
    fits_write_key(fptr, TDOUBLE, "CDELT2", &dm, "", &status);
    fits_write_key(fptr, TSTRING, "CUNIT2", cunit, "", &status);

    char ctype3[] = "FREQ";
    fits_write_key(fptr, TSTRING, "CTYPE3", ctype3, "", &status);
    crpix1 = 1; 
    fits_write_key(fptr, TDOUBLE, "CRPIX3", &crpix1, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL3", &freq, "", &status); 
    fits_write_key(fptr, TDOUBLE, "CDELT3", &deltaf, "", &status);
    char cunit3[]="Hz";
    fits_write_key(fptr, TSTRING, "CUNIT3", cunit3, "", &status);

    char ctype4[] = "STOKES";
    fits_write_key(fptr, TSTRING, "CTYPE4", ctype4, "", &status);
    double one=1.0;
    crpix1 = (double)stokes; 
    fits_write_key(fptr, TDOUBLE, "CRPIX4", &one, "", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL4", &crpix1, "", &status); 
    fits_write_key(fptr, TDOUBLE, "CDELT4", &one, "", &status);
    char cunit4[]="";
    fits_write_key(fptr, TSTRING, "CUNIT4", cunit4, "", &status);


    // auxilliary  info
    char origin[] = "ExCon (C) Sarod Yatawatta";
    fits_write_key(fptr, TSTRING, "ORIGIN", origin, "", &status);
    // weights
    fits_write_key(fptr, TFLOAT, "WEIGHT", &sumweight, "", &status);

    //write the pixel data to the file
    long fpixel[4] = {1, 1, 1, 1};
    long nelements = w * h;
    fits_write_pix(fptr, TFLOAT, fpixel, nelements, ptr, &status);
    //close the file
    fits_close_file(fptr, &status);

    if(status)
        fits_report_error(stderr, status);
}


bool Cmd::W_Projection = true;
char *Cmd::TableName = NULL;
String Cmd::DataField = "CORRECTED_DATA";
double Cmd::TimeInterval = 0.0; //minutes, converted to seconds
int Cmd::Nx = 1024;
int Cmd::Ny = Cmd::Nx;
int Cmd::Stacks = 128;
bool Cmd::HasPointing = false;
float Cmd::PointingRef[2] = {0.0f, 0.0f};
float Cmd::UVscale = 5;
float Cmd::ImRes = 2.0f;
char Cmd::Stokes = 'I'; //I,Q,U,V
int Cmd::AvgChannels = 0; /* 0:MFS, 1: average, 2: channel imaging */
int Cmd::Mc= 4; 
int Cmd::variable_M= 0; /* if 1, enable variable conv. kernel bandwitdh with freq */ 
int Cmd::Nt= 16;  /* no of threads */
int Cmd::Npc= 256; 
float Cmd::Nzc= 2.0f; /* padded size Npc(1+2*Nzc) */ 
float Cmd::padding=1.25f; /* image padding ratio */
float Cmd::robustval=0.0f; /* robust weight */
int Cmd::imgmode=IMG_I; /* image mode, I|QUV and whether to create weights,PSF FITS files */
float Cmd::conv_cutoff=0.01f;
float Cmd::apod_cutoff=0.001f;
int Cmd::wgrid_step=0; /* how to sample w values, 0: sqrt(), 1: log(), 2: linear,  3: w^(expW) spacing, 4: user supplied text file, */
float Cmd::expW=0.5f; /* w^(expW) spacing */
const char *Cmd::File_wstep=NULL;
int Cmd::use_gpu=1; /* >0 to use GPU */
int Cmd::create_imaginary=0; /* >0 to create imaginary part */
int Cmd::gpu_threads=96; /* no of worker threads per GPU (each will create more threads) */
int Cmd::w_snapshot=0; /* if >0 enable W snapshots */

int Cmd::buckets=1; /* if > 1, sort data using this many buckets */

/* default binary file locations */
const char *Cmd::File_uvgrid=NULL;
const char *Cmd::File_imgrid=NULL;
const char *Cmd::File_wgrid=NULL;
const char *Cmd::file_qualifier=NULL;
int Cmd::weightmode=WMODE_ROBUST;
int Cmd::pm_maxiter=3;
int Cmd::pm_convmode=0;
float Cmd::pm_imscale=0.0f; /* default is uniform for 0.0 */

int Cmd::convkernel=0; /* if >0, write convolution kernel as a FITS file */

void
print_copyright(void) {
  cout<<"ExCon 0.0.31 (C) 2013-2017 Sarod Yatawatta"<<endl;
}


void
print_help(void) {
   cout << "Usage:" << endl;
   cout << "-m MS name" << endl;
   cout << "-c column (DATA/CORRECTED_DATA/...) : default " <<Cmd::DataField<< endl;
   cout << "-p pixel width (arcsec) : default "<<Cmd::ImRes << endl;
   cout << "-d Npix image dimension in pixels : default "<<Cmd::Nx << endl;
   cout << "-w no of Wprojection planes : default "<<Cmd::Stacks << endl;
   cout << "-l lower uv cutoff (lambda) : default "<<Data::min_uvcut<<endl;
   cout << "-u upper uv cutoff (lambda) : default "<<Data::max_uvcut<<endl;
   cout << "-M 0,1,2 : 0: average channels, 1: MFS mode, 2: image each channel: default "<<(!Cmd::AvgChannels?1:(Cmd::AvgChannels==1?0:2))<< endl;
   cout << "-W upper abs(W) cutoff (lambda) : default "<<Data::max_wcut<<endl;
   cout << "-R weighting scheme 0:uniform, 1:robust, 2:natural 3:Adaptive 4:uniform+taper: default "<<Cmd::weightmode<<endl;
   cout << "-r robust weight parameter [-2 (uniform),2 (natural)] (only when R=1) : default "<<Cmd::robustval<<endl;
   cout << "-T snapshot time interval in minutes (0 for full integration): default "<<Cmd::TimeInterval<<endl;
   cout<<endl;
   cout << "Advanced options:"<<endl;
   cout << "-t no of worker threads : default "<<Cmd::Nt << endl;
   cout << "-b convolutional kernel bandwidth in uv pixels : default "<<Cmd::Mc << endl;
   cout << "-f 0,1 : if 1, convolutional kernel bandwidth change with frequency: default "<<Cmd::variable_M << endl;
   cout << "-s convolutional kernel size in pixels : default "<<Cmd::Npc << endl;
   cout << "-z convolutional kernel zero padding size (as a factor 1 + 2 x ?) : default "<<Cmd::Nzc << endl;
   cout << "-a image zero padding ratio (as a factor) : default "<<Cmd::padding<< endl;

   cout << "-x 0,1,2,3,4 : 0: I, 1: I+extra (I Grid), 2:IQUV, 3:IQUV+extra (I Grid), 4:IQUV+extra (IQUV Grid), extra means creating Weights,PSF,GRID as FITS files : default "<<Cmd::imgmode<< endl;
   cout << "-e convolutional kernel truncation, as a fraction of peak value : default "<<Cmd::conv_cutoff<< endl;
   cout << "-o apodization correction cutoff: default "<<Cmd::apod_cutoff<< endl;
   cout << "-g 0,1 : if 1, use GPU for gridding : default "<<Cmd::use_gpu<< endl;
   cout << "-q 0,1 : if 1, produce imaginary parts of all images/grid : default "<<Cmd::create_imaginary<< endl;
   cout << "-y worker threads per GPU for gridding : default "<<Cmd::gpu_threads<< endl;
   cout << "-i howto sample w axis 0: sqrt(), 1: log(),  2: linear,  3: power sampling with user supplied exponent -Z, 4: user supplied with -H, default "<<Cmd::wgrid_step<< endl;
   cout << "-Z sample w axis as w^(expW): only used if -i 3 default "<<Cmd::expW<< endl;
   cout << "-H filename: gives the no. of w planes and the spacing (used with  -i 4)"<< endl;
   cout << "-A | -B | -C : binary file names (on separate disks) : default current directory"<< endl;
   cout << "-Q string: add unique qualifier to image file names : default none"<< endl;
   cout << "-S 0,1 : if 1, enable W-snapshots: default "<<Cmd::w_snapshot<< endl;
   cout << "-P hh,mm,ss,dd,mm,ss : RA,Dec : if given will shift phase center to this: default: none"<< endl;
   cout << "-D PM iterations: default "<<Cmd::pm_maxiter<< endl;
   cout << "-G 0,1,2,3 : PM convolution mode: 0: uniform, 1: NCP (50-250L) use 30-300L input, 2: NCP (0-40kL), 3: 50-250L disk, default "<<Cmd::pm_convmode<< endl;
   cout << "-E uvcut,scale_factor: for rescaling inner baselines: default "<<Data::rescale_uvcut<<","<<Data::rescale_factor<< endl;
   cout << "-F XX,XY,YX,YY: if given (MFS mode), flag data using these as upper absolute limits: default "<<Data::maxXY[0]<<","<<Data::maxXY[1]<<","<<Data::maxXY[2]<<","<<Data::maxXY[3]<< endl;
   cout << "-X XX,XY,YX,YY: if given, find standard deviation, and flag data using these as upper limits in sigma: default "<<Data::sigXY[0]<<","<<Data::sigXY[1]<<","<<Data::sigXY[2]<<","<<Data::sigXY[3]<< endl;
   cout << "-K no. of buckets to sort data, if K>1: default "<<Cmd::buckets<< endl;
   cout << "-V 0,1 : if 1, write convolution kernel as a FITS file: default "<<Cmd::convkernel<< endl;
   cout<<endl;
   cout<<"Output: image_I|Q|U|V.fits image_W.fits (weights) image_PSF.fits (PSF) image_GR|I.fits (gridded data)  *i.fits (imaginary part)"<<endl;
   cout<<endl;
   cout <<"Report bugs to <sarod@users.sf.net>"<<endl;
}

void 
Cmd::ParseCmdLine(int ac, char **av) {
    print_copyright();
    char c;
    if(ac < 2)
    {
        print_help();
        exit(0);
    }
    W_Projection = true;
    char *ptr;
    while((c=getopt(ac, av, "A:B:C:D:E:F:G:H:K:L:M:P:Q:R:S:T:V:W:X:Z:a:b:c:d:e:f:g:i:l:m:o:p:q:r:s:t:u:w:x:y:z:h"))!= -1)
    {
        switch(c)
        {
            case 'p':
                ImRes = atof(optarg);
                break;
            case 'm':
                TableName = optarg;
                break;
            case 'A':
                File_uvgrid= optarg;
                break;
            case 'Q':
                file_qualifier= optarg;
                break;
            case 'B':
                File_imgrid= optarg;
                break;
            case 'C':
                File_wgrid= optarg;
                break;
            case 'S':
                w_snapshot= atoi(optarg);
                break;
            case 'T':
                TimeInterval = atof(optarg);
                break;
            case 'w': 
                Stacks = atoi(optarg);
                break;
            case 'i': 
                wgrid_step= atoi(optarg);
                if (wgrid_step>4) { 
                 print_help(); exit(1);
                }
                break;
            case 'Z':
                expW=atof(optarg);
                break;
            case 'H':
                File_wstep= optarg;
                break;
            case 'c': 
                DataField = optarg;
                break;
            case 'd':
                ptr = strchr(optarg, 'x');
                if(ptr) {
                    *ptr = '\0';
                    Nx = atoi(optarg);
                    Ny = atoi(ptr+1);
                } else {
                    Nx = atoi(optarg);
                    Ny = Nx;
                }
                break;
            case 't': 
                Nt= atoi(optarg);
                break;
            case 'b': 
                Mc= atoi(optarg);
                break;
            case 'f': 
                variable_M= atoi(optarg);
                break;
            case 'x': 
                imgmode= atoi(optarg);
                if (imgmode>IMG_IQUVF) { imgmode=IMG_I; }
                break;
            case 's': 
                Npc= atoi(optarg);
                break;
            case 'z': 
                Nzc= atof(optarg);
                break;
            case 'a': 
                padding= atof(optarg);
                break;
            case 'l': 
                Data::min_uvcut= atof(optarg);
                break;
            case 'u': 
                Data::max_uvcut= atof(optarg);
                break;
            case 'W': 
                Data::max_wcut= atof(optarg);
                break;
            case 'r':
                robustval= atof(optarg);
                break;
            case 'e':
                conv_cutoff= atof(optarg);
                break;
            case 'o':
                apod_cutoff= atof(optarg);
                break;
            case 'g':
                use_gpu= atoi(optarg);
                break;
            case 'q':
                create_imaginary= atoi(optarg);
                break;
            case 'y':
                gpu_threads= atoi(optarg);
                break;
            case 'R': 
                weightmode= atoi(optarg);
                break;
            case 'D': 
                pm_maxiter= atoi(optarg);
                break;
            case 'G': 
                pm_convmode= atoi(optarg);
                break;
            case 'K': 
                buckets= atoi(optarg);
                break;
            case 'V': 
                convkernel= atoi(optarg);
                break;
            case 'L': 
                /* arcminutes */
                pm_imscale=atof(optarg);
                break;
            case 'M': 
            {    
                int tmpachn=atoi(optarg);
                if (tmpachn==0) { AvgChannels=1; }
                if (tmpachn==1) { AvgChannels=0; }
                if (tmpachn==2) { AvgChannels=2; }
            }
                break;
            case 'P': //phase shift to new center
            {
                HasPointing = true;
                char rahr_[128];
                char ramin_[128];
                char rasec_[128];
                char decd_[128];
                char decmin_[128];
                char decsec_[128];
                int matched=sscanf(optarg,"%127[^,],%127[^,],%127[^,],%127[^,],%127[^,],%127[^,]",rahr_,ramin_,rasec_,decd_,decmin_,decsec_);
                if(matched==6) {
                   double rahr,ramin,rasec,decd,decmin,decsec;
                   sscanf(rahr_,"%lf",&rahr);
                   sscanf(ramin_,"%lf",&ramin);
                   sscanf(rasec_,"%lf",&rasec);
                   sscanf(decd_,"%lf",&decd);
                   sscanf(decmin_,"%lf",&decmin);
                   sscanf(decsec_,"%lf",&decsec);
                   if (rahr<0) {
                     rahr=-(-rahr+ramin/60.0+rasec/3600.0)*M_PI/12.0;
                   } else {
                     rahr=(rahr+ramin/60.0+rasec/3600.0)*M_PI/12.0;
                   }
                   if (decd<0) {
                     decd=-(-decd+decmin/60.0+decsec/3600.0)*M_PI/180.0;
                   } else {
                     decd=(decd+decmin/60.0+decsec/3600.0)*M_PI/180.0;
                   }
                   PointingRef[0] = (float)(rahr);
                   PointingRef[1] = (float)(decd);
                } else{
                    cout << "Error in RA,Dec format : hh,mm,ss,dd,mm,ss." << endl;
                    exit(0);
                }
                break;
            }
            case 'E': //find uvcut,scalefactor for rescaling
            {
                char uvc_[128];
                char rescale_[128];
                int matched=sscanf(optarg,"%127[^,],%127[^,]",uvc_,rescale_);
                if(matched==2) {
                   sscanf(uvc_,"%f",&Data::rescale_uvcut);
                   sscanf(rescale_,"%f",&Data::rescale_factor);
                } else {
                    cout << "Error in uvcut,scale_factor format" << endl;
                    exit(0);
                }
                break;
            }
            case 'F': //find XX,XY,YX,YY abs limit values to flag
            {
                char xx_[128];
                char xy_[128];
                char yx_[128];
                char yy_[128];
                int matched=sscanf(optarg,"%127[^,],%127[^,],%127[^,],%127[^,]",xx_,xy_,yx_,yy_);
                if(matched==4) {
                   Data::clipval=true;
                   sscanf(xx_,"%f",&Data::maxXY[0]);
                   sscanf(xy_,"%f",&Data::maxXY[1]);
                   sscanf(yx_,"%f",&Data::maxXY[2]);
                   sscanf(yy_,"%f",&Data::maxXY[3]);
                } else {
                    cout << "Error in data limits XX,XY,YX,YY format" << endl;
                    exit(0);
                }
                break;
            }
            case 'X': //find sigma(XX,XY,YX,YY) limit to flag
            {
                char xx_[128];
                char xy_[128];
                char yx_[128];
                char yy_[128];
                int matched=sscanf(optarg,"%127[^,],%127[^,],%127[^,],%127[^,]",xx_,xy_,yx_,yy_);
                if(matched==4) {
                   Data::flagsigma=true;
                   sscanf(xx_,"%f",&Data::sigXY[0]);
                   sscanf(xy_,"%f",&Data::sigXY[1]);
                   sscanf(yx_,"%f",&Data::sigXY[2]);
                   sscanf(yy_,"%f",&Data::sigXY[3]);
                } else {
                    cout << "Error in sigma limits XX,XY,YX,YY format" << endl;
                    exit(0);
                }
                break;
            }
            case 'h': 
                print_help();
                exit(1);
            default:
                print_help();
                exit(1);
        }
    }

   if (!TableName) {
    cout<<"An MS name must be given."<<endl;
    print_help();
    exit(1);
   }
   if (!File_wstep && wgrid_step==4) {
    cout<<"Text file with w values should be provided for -i 3"<<endl;
    exit(1);
   }
   if (wgrid_step==3) {
    cout<<"Sampling w axis as w^("<<expW<<")"<<endl;
   }
    cout<<"Image Size: "<<Nx<<"x"<<Ny<<" Pixels: "<<ImRes<<"arcsec"<<" Wplanes: "<<Stacks<<" MS: "<<TableName<<endl;
    cout<<"Selecting baselines in ["<<Data::min_uvcut<<","<<Data::max_uvcut<<","<<Data::max_wcut<<"] wavelengths";
 if (weightmode==WMODE_UNIFORM) {
  cout<<", Uniform weight"<<endl;
 } else if (weightmode==WMODE_NATURAL) {
  cout<<", Natural weight"<<endl;
 } else if (weightmode==WMODE_ROBUST) {
  cout<<", Robust weight : "<<robustval<<endl;
 } else if (weightmode==WMODE_PIPEMENON) {
  cout<<", Adaptive weighting iterations : "<<pm_maxiter<<", convolution mode (also for tapering) : "<<pm_convmode<<endl;
 } else {
  cout<<", Uniform weight and tapering with function no. "<<pm_convmode<<"."<<endl;
 }
 if (Data::rescale_uvcut!=0.0f) {
  cout<<"Rescaling baselines shorter than "<<Data::rescale_uvcut<<" by "<<Data::rescale_factor<<endl;
 }
 /* convert pm_imscale from arcminutes to radians */
 pm_imscale*=M_PI/(18000.0f);
 if (TimeInterval>0.0) {
  cout<<"Snapshot imaging with interval "<<TimeInterval<<" minutes."<<endl;
  /* convert interval to seconds */
  TimeInterval *=60.0;
 }
 if (AvgChannels==1) {
  cout<<"Image mode: average channels"<<endl;
 } else if (AvgChannels==2) {
  cout<<"Image mode: individual channel"<<endl;
 } else {
  cout<<"Image mode: MFS"<<endl;
 } 
}
