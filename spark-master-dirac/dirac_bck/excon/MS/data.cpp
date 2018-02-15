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
#include <measures/Measures/MDirection.h>
#include <measures/Measures/UVWMachine.h>
#include <casa/Quanta.h>
#include <casa/Quanta/Quantum.h>
#include <casa/Arrays/MatrixMath.h>
#include <tables/TaQL/TableParse.h> // old tables/Tables/

using namespace casa;

//initialize the data variables.
iodata **Data::data = NULL;
float *Data::chanFreqs= NULL;
float *Data::chanBWs= NULL;
float Data::minW = 0.0f;
float Data::maxW = 0.0f; //wavelengths
float Data::maxUV[3] = {0.0f,0.0f,0.0f};
double Data::refDir[2] = {0, 0};
int Data::numChannels = 1; 
float Data::chanFreq = 160e6f; //MHz
float Data::chanBW = 200e3f; //kHz
float Data::waveLen = 1.8f; //meters
unsigned int Data::numRows = 0;
float Data::min_uvcut=0.0f;
float Data::max_uvcut=1e9f;
float Data::max_wcut=1e9f;

bool Data::clipval = false;
float Data::maxXY[4] = {0.0f,0.0f,0.0f,0.0f};

bool Data::flagsigma= false;
float Data::sigXY[4] = {0.0f,0.0f,0.0f,0.0f};

float Data::rescale_uvcut=0.0f;
float Data::rescale_factor=1.0f;

// get complex number from CASA format
static fcomp 
getfComplex(Complex c) {
    fcomp tmp;
    tmp.x = (float)c.real();
    tmp.y = (float)c.imag();
    return tmp;
}

// get flipped complex number from CASA format
// (X+jY) becomes -Y+jX
static fcomp 
getfComplexConj(Complex c) {
    fcomp tmp;
    tmp.x = -(float)c.imag();
    tmp.y = (float)c.real();
    return tmp;
}

// get absolute value CASA format
static float
getfAbs(Complex c) {
    return sqrtf(c.real()*c.real()+c.imag()*c.imag());
}


using namespace Data;

static void 
fillTransMatrix (Matrix<double>& mat, double ra, double dec) {
      double sinra;
      double cosra;
      double sindec;
      double cosdec;
      sincos(ra,&sinra,&cosra);
      sincos(dec,&sindec,&cosdec);
      mat(0,0) = cosra;
      mat(1,0) = -sinra;
      mat(2,0) = 0.0;
      mat(0,1) = -sinra*sindec;
      mat(1,1) = -cosra*sindec;
      mat(2,1) = cosdec;
      mat(0,2) = sinra*cosdec;
      mat(1,2) = cosra*cosdec;
      mat(2,2) = sindec;
}


// avg: 0: no averaging MFS, 1: average to one channel, 2: channel imaging
static iodata ** 
loadFromMS(Table ti, char *fname, int avg=0,bool imaginarypart=false, bool sort=false) {

    char buff[2048] = {0};
    /* sort input to iterate over baselines:
      to preserve some locality in grid io  */
    Table t;
    if (sort) {
     sprintf(buff, "select from %s orderby sumsqr (UVW[0]),(UVW[1])", fname);
     t=tableCommand(buff,ti);
    } else {
     t=ti;
    }

    ROArrayColumn<Complex> dataCol(t, Cmd::DataField);
    ROArrayColumn<double> uvwCol(t, "UVW"); 
    ROArrayColumn<bool> flagCol(t, "FLAG");
    ROScalarColumn<int> a1(t, "ANTENNA1"), a2(t, "ANTENNA2");
    numRows = t.nrow();
    Data::numChannels = dataCol.shape(0)[1];
    iodata **out = NULL;

    //obtain the field refdir infor
    Table _field = Table(t.keywordSet().asTable("FIELD"));
    ROArrayColumn<double> ref_dir(_field, "PHASE_DIR"); //old code used REFERENCE_DIR
    Array<double> dir = ref_dir(0);
    double *c = dir.data();
    refDir[0] = c[0];
    refDir[1] = c[1];

    //Create MDirection object of the original Phase center.
    if(Cmd::HasPointing){
        refDir[0] = Cmd::PointingRef[0];
        refDir[1] = Cmd::PointingRef[1];
    }

      double newRa  = refDir[0];
      double newDec = refDir[1];
      double oldRa  = c[0];
      double oldDec = c[1];
      Matrix<double> oldUVW(3,3);
      Matrix<double> newUVW(3,3);
      Matrix<double> M(3,3);
      double XYZ[3]={0.0,0.0,0.0};
      const double* mat1 = NULL;

    if(Cmd::HasPointing){
      fillTransMatrix (oldUVW, oldRa, oldDec);
      fillTransMatrix (newUVW, newRa, newDec);
      M.reference(product(transpose(newUVW), oldUVW));
      Matrix<double> wold(oldUVW(IPosition(2,0,2),IPosition(2,2,2)));
      Matrix<double> wnew(newUVW(IPosition(2,0,2),IPosition(2,2,2)));
      Matrix<double> tt= product(transpose(Matrix<double>(wold-wnew)), oldUVW);
      XYZ[0] = tt(0,0);
      XYZ[1] = tt(0,1);
      XYZ[2] = tt(0,2);
      mat1 =M.data();
    }

    /* for channel imaging,setup freq/BW arrays */
    if (avg==2) {
     Data::chanFreqs=new float[numChannels];
     Data::chanBWs=new float[numChannels];
    } 

    //obtain the chanel freq information
    Table _freq = Table(t.keywordSet().asTable("SPECTRAL_WINDOW"));
    ROArrayColumn<double> chan_freq(_freq, "CHAN_FREQ"); 
    /* need channel widths to calculate bandwidth */
    ROArrayColumn<double> chan_width(_freq, "CHAN_WIDTH");
    chanBW=(float)numChannels*(chan_width(0).data()[0]);

    // find average of frequencies
    chanFreq=0.0;
    for (int ci=0; ci<numChannels; ci++) {
     chanFreq += chan_freq(0).data()[ci];
     if (avg==2) {
       chanFreqs[ci]=(float)chan_freq(0).data()[ci];
       chanBWs[ci]=(float)chan_width(0).data()[ci];
     }
    }
    chanFreq/=(double)numChannels;

    double speed_light = C::c;
    double invfreq=1.0/chanFreq;
    waveLen = speed_light*invfreq;
    double invwl=1.0/waveLen;
    double freqTerm = 2.0 * C::pi * chanFreq/C::c; 
 
    if(avg==1) {
        out = new iodata* [1];
        out[0] = new iodata[numRows]();
    } else {
        out = new iodata* [1];
        out[0] = new iodata[numRows*numChannels]();
    }

    float chanscale=1.0f/(float)numChannels;
    
    bool rescale_data=(rescale_uvcut>0.0f && rescale_factor!=1.0f?true:false);

    for(uInt row = 0; row < numRows; row++) {
        uInt i = a1(row); //antenna1 
        uInt j = a2(row); //antenna2 
        Array<Complex> data = dataCol(row);
        Matrix<double> uvw = uvwCol(row);
        Array<bool> flag = flagCol(row);
        Complex cxx(0.0f, 0.0f);
        Complex cxy(0.0f, 0.0f);
        Complex cyx(0.0f, 0.0f);
        Complex cyy(0.0f, 0.0f);
        if(avg!=1) {
/***********************************************************************/
         /* calculate sqrt(u^2+v^2) to select uv cuts */
         double *c = uvw.data();
         bool flag_uvcut=0;
         if (i==j) {
           flag_uvcut=true;
         } 

         if( flag_uvcut ) {
            for(int k = 0; k < numChannels; k++) {
                int offset=(avg==2?k*numRows+row : row*numChannels+k);
                out[0][offset].xx = getfComplex(cxx);
                out[0][offset].xy = getfComplex(cxy);
                out[0][offset].yx = getfComplex(cyx);
                out[0][offset].yy = getfComplex(cyy);
                out[0][offset].u = 0.0f;
                out[0][offset].v = 0.0f;
                out[0][offset].w = 0.0f;
                out[0][offset].flag=1;
            }
         } else {
         //place visibility data suitable for MFS mode, change uvw coordinates
         // for freq0 
           float uvd=(float)sqrt(c[0]*c[0]+c[1]*c[1])*invwl;
           float wd=(float)fabs(c[2])*invwl;
           for(int k = 0; k < numChannels; k++) {
             double frat=(chan_freq(0).data()[k])*invfreq;
             int offset=(avg==2?k*numRows+row : row*numChannels+k);
             /* re-test for uv distance flagging */
             float uvd1=(float)uvd*frat;
             float wd1=(float)wd*frat;
             bool *flgptr=flag[k].data(); /* check flag of each channel */
             bool flag_uvcut1=(flag_uvcut||flgptr[0]||flgptr[1]||flgptr[2]||flgptr[3]);
                if (uvd1<=min_uvcut || uvd1>max_uvcut || wd1>max_wcut ) {
                 flag_uvcut1=true;
                } 
                /* check to see if we need to clip */

                Complex *ptr = data[k].data();
                if (!flag_uvcut1 && clipval) {
                  flag_uvcut1=(maxXY[0]<getfAbs(ptr[0])||maxXY[3]<getfAbs(ptr[3])||maxXY[1]<getfAbs(ptr[1])||maxXY[2]<getfAbs(ptr[2]));
                }
                if( flag_uvcut1 ) {
                    out[0][offset].xx = getfComplex(cxx);
                    out[0][offset].xy = getfComplex(cxy);
                    out[0][offset].yx = getfComplex(cyx);
                    out[0][offset].yy = getfComplex(cyy);
                    out[0][offset].u = 0.0f;
                    out[0][offset].v = 0.0f;
                    out[0][offset].w = 0.0f;
                    out[0][offset].flag=1;
                } else {
                   /* determine if we need to rescale this datapoint (not flagged) */
                   bool scaleuv=false;
                   if (rescale_data && uvd1<=rescale_uvcut) {
                     scaleuv=true;
                   }

                    double ut,vt,wt;
                    if(Cmd::HasPointing){
                      //find the new phase-shifted UVW to use.
                       double u1=c[0]*frat;
                       double v1=c[1]*frat;
                       double w1=c[2]*frat;
                       ut = u1*mat1[0] + v1*mat1[3] + w1*mat1[6];
                       vt = u1*mat1[1] + v1*mat1[4] + w1*mat1[7];
                       wt = u1*mat1[2] + v1*mat1[5] + w1*mat1[8];
                       double phase = XYZ[0]*u1 + XYZ[1]*v1 + XYZ[2]*w1;
                       double phasewvl = phase * freqTerm;
                       double cosph,sinph;
                       sincos(phasewvl,&sinph,&cosph);
                       Complex phasor((float)cosph, (float)sinph);
                      if (!imaginarypart) {
                        out[0][offset].xx = getfComplex(ptr[0]*phasor);
                        out[0][offset].xy = getfComplex(ptr[1]*phasor);
                        out[0][offset].yx = getfComplex(ptr[2]*phasor);
                        out[0][offset].yy = getfComplex(ptr[3]*phasor);
                      } else {
                        out[0][offset].xx = getfComplexConj(ptr[0]*phasor);
                        out[0][offset].xy = getfComplexConj(ptr[1]*phasor);
                        out[0][offset].yx = getfComplexConj(ptr[2]*phasor);
                        out[0][offset].yy = getfComplexConj(ptr[3]*phasor);
                      }
                    } else {
                       ut=c[0]*frat;
                       vt=c[1]*frat;
                       wt=c[2]*frat;
                      if (!imaginarypart) {
                        out[0][offset].xx = getfComplex(ptr[0]);
                        out[0][offset].xy = getfComplex(ptr[1]);
                        out[0][offset].yx = getfComplex(ptr[2]);
                        out[0][offset].yy = getfComplex(ptr[3]);
                       } else {
                        out[0][offset].xx = getfComplexConj(ptr[0]);
                        out[0][offset].xy = getfComplexConj(ptr[1]);
                        out[0][offset].yx = getfComplexConj(ptr[2]);
                        out[0][offset].yy = getfComplexConj(ptr[3]);
                       }
                    }
                    if(maxW < wt)
                        maxW = wt;
                    if(minW > wt)
                        minW = wt;
                    /* always store u,v,w in wavelengths */
                    out[0][offset].u = (float) ut*invwl;
                    out[0][offset].v = (float) vt*invwl;
                    out[0][offset].w = (float) wt*invwl;
                    out[0][offset].flag=0;
                    maxUV[0] = maxUV[0] < fabs(ut)? fabs(ut) : maxUV[0];
                    maxUV[1] = maxUV[1] < fabs(vt)? fabs(vt) : maxUV[1];
                    float dist =sqrtf(ut*ut + vt*vt);
                    maxUV[2] = (dist < maxUV[2]? maxUV[2] : dist);
                    /* if rescaling is needed, rescale data */
                    if (scaleuv) {
                     out[0][offset].xx.x *=rescale_factor;
                     out[0][offset].xx.y *=rescale_factor;
                     out[0][offset].xy.x *=rescale_factor;
                     out[0][offset].xy.y *=rescale_factor;
                     out[0][offset].yx.x *=rescale_factor;
                     out[0][offset].yx.y *=rescale_factor;
                     out[0][offset].yy.x *=rescale_factor;
                     out[0][offset].yy.y *=rescale_factor;
                    }
                }
            }
         }
/***********************************************************************/
        } else { //average all channels 
/***********************************************************************/
            /* calculate sqrt(u^2+v^2) to select uv cuts */
            double *c = uvw.data();
            float uvd=(float) sqrt(c[0]*c[0]+c[1]*c[1])*invwl;
            float wd=(float) fabs(c[2])*invwl;
            bool flag_uvcut=0;
            if (i==j || uvd<=min_uvcut || uvd>max_uvcut || wd>max_wcut ) {
              flag_uvcut=true;
            } 
            /* for averaging, all channels should be unflagged */
            for(int k = 0; k < numChannels; k++) {
             bool *flgptr=flag[k].data(); /* check flag of each channel */
             flag_uvcut=(flag_uvcut||flgptr[0]||flgptr[1]||flgptr[2]||flgptr[3]);
            }

            if( flag_uvcut ) {
                out[0][row].xx = getfComplex(cxx);
                out[0][row].xy = getfComplex(cxy);
                out[0][row].yx = getfComplex(cyx);
                out[0][row].yy = getfComplex(cyy);
                out[0][row].u = 0.0f;
                out[0][row].v = 0.0f;
                out[0][row].w = 0.0f;
                out[0][row].flag = 1;
            } else {
                bool scaleuv=false;
                if (rescale_data && uvd<=rescale_uvcut) {
                     scaleuv=true;
                }


                VectorIterator<double> vecIter(uvw, 0);
                Vector<double> vec_uvw = vecIter.vector();
                if(Cmd::HasPointing){
                    //find the new phase-shifted UVW to use.
                   double ut = vec_uvw[0]*mat1[0] + vec_uvw[1]*mat1[3] + vec_uvw[2]*mat1[6];
                   double vt = vec_uvw[0]*mat1[1] + vec_uvw[1]*mat1[4] + vec_uvw[2]*mat1[7];
                   double wt = vec_uvw[0]*mat1[2] + vec_uvw[1]*mat1[5] + vec_uvw[2]*mat1[8];
                   double phase = XYZ[0]*vec_uvw[0] + XYZ[1]*vec_uvw[1] + XYZ[2]*vec_uvw[2];
                   vec_uvw[0]=ut;
                   vec_uvw[1]=vt;
                   vec_uvw[2]=wt;
                    double phasewvl = phase * freqTerm;
                    double cosph,sinph;
                    sincos(phasewvl,&sinph,&cosph);
                    Complex phasor((float)cosph, (float)sinph);
                    for(int k = 0; k < numChannels; k++) {
                    Complex *ptr = data[k].data();
                    cxx += ptr[0] * phasor;
                    cxy += ptr[1] * phasor;
                    cyx += ptr[2] * phasor;
                    cyy += ptr[3] * phasor;
                    }
                } else {
                    for(int k = 0; k < numChannels; k++) {
                    Complex *ptr = data[k].data();
                    cxx += ptr[0];
                    cxy += ptr[1];
                    cyx += ptr[2];
                    cyy += ptr[3];
                    }
                }
                if (numChannels>1) {
                 if (!imaginarypart) {
                  out[0][row].xx.x = (float)cxx.real()*chanscale;
                  out[0][row].xx.y = (float)cxx.imag()*chanscale;
                  out[0][row].xy.x = (float)cxy.real()*chanscale;
                  out[0][row].xy.y = (float)cxy.imag()*chanscale;
                  out[0][row].yx.x = (float)cyx.real()*chanscale;
                  out[0][row].yx.y = (float)cyx.imag()*chanscale;
                  out[0][row].yy.x = (float)cyy.real()*chanscale;
                  out[0][row].yy.y = (float)cyy.imag()*chanscale;
                 } else {
                  /* flip data real,imag  */
                  out[0][row].xx.x = (float)cxx.imag()*chanscale;
                  out[0][row].xx.y = (float)cxx.real()*chanscale;
                  out[0][row].xy.x = (float)cxy.imag()*chanscale;
                  out[0][row].xy.y = (float)cxy.real()*chanscale;
                  out[0][row].yx.x = (float)cyx.imag()*chanscale;
                  out[0][row].yx.y = (float)cyx.real()*chanscale;
                  out[0][row].yy.x = (float)cyy.imag()*chanscale;
                  out[0][row].yy.y = (float)cyy.real()*chanscale;
                 }
                } else {
                 if (!imaginarypart) {
                  out[0][row].xx.x = (float)cxx.real();
                  out[0][row].xx.y = (float)cxx.imag();
                  out[0][row].xy.x = (float)cxy.real();
                  out[0][row].xy.y = (float)cxy.imag();
                  out[0][row].yx.x = (float)cyx.real();
                  out[0][row].yx.y = (float)cyx.imag();
                  out[0][row].yy.x = (float)cyy.real();
                  out[0][row].yy.y = (float)cyy.imag();
                 } else {
                  /* flip data real,imag */
                  out[0][row].xx.x = (float)cxx.imag();
                  out[0][row].xx.y = (float)cxx.real();
                  out[0][row].xy.x = (float)cxy.imag();
                  out[0][row].xy.y = (float)cxy.real();
                  out[0][row].yx.x = (float)cyx.imag();
                  out[0][row].yx.y = (float)cyx.real();
                  out[0][row].yy.x = (float)cyy.imag();
                  out[0][row].yy.y = (float)cyy.real();
                 }
                }

                maxUV[0] = maxUV[0] < fabs(vec_uvw[0])? fabs(vec_uvw[0]) : maxUV[0];
                maxUV[1] = maxUV[1] < fabs(vec_uvw[1])? fabs(vec_uvw[1]) : maxUV[1];
                float dist =(float)sqrt(vec_uvw[0]*vec_uvw[0] + vec_uvw[1]*vec_uvw[1]);
                maxUV[2] = (dist < maxUV[2]? maxUV[2] : dist);

                maxW = maxW < vec_uvw[2]? vec_uvw[2] : maxW;
                minW = minW > vec_uvw[2]? vec_uvw[2] : minW;

                out[0][row].u = (float)vec_uvw[0]*invwl;
                out[0][row].v = (float)vec_uvw[1]*invwl;
                out[0][row].w = (float)vec_uvw[2]*invwl;
                out[0][row].flag = 0;
                /* if rescaling is needed, rescale data */
                if (scaleuv) {
                   out[0][row].xx.x *=rescale_factor;
                   out[0][row].xx.y *=rescale_factor;
                   out[0][row].xy.x *=rescale_factor;
                   out[0][row].xy.y *=rescale_factor;
                   out[0][row].yx.x *=rescale_factor;
                   out[0][row].yx.y *=rescale_factor;
                   out[0][row].yy.x *=rescale_factor;
                   out[0][row].yy.y *=rescale_factor;
                }
            }
/***********************************************************************/
        }
    }
    /* effective no of rows gets multiplied by channels for MFS (avg==0) */
    if (!avg) numRows=numRows*numChannels;
    return out; 
}

void 
loadData(Table t, char *fname, int avg, bool imaginarypart, bool sort) {
       data = loadFromMS(t, fname, avg); /* do not enable sort or imaginary transform */
}

void 
freeData() {
    if(data)
    {
        delete [] data[0]; /* we have only one column, even for multi channel data */
        delete [] data;
    }
    if (chanFreqs) {
      delete [] chanFreqs;
    }
    if (chanBWs) {
      delete [] chanBWs;
    }

}
