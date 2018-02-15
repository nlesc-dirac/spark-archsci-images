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


#include <wcs.h>

#include "largefft.h"
/* 
  create a default WCS structure to calculate l,m coords
  deltalm: pixel width (deg)
  int N0: reference pixel 
  ra0,dec0: reference coordinates (deg)
  only if gid==0, enable memory allocation for WCS
*/
int
generate_def_wcs(struct wcsprm *wcs, double deltalm, double N0, double ra0, double dec0, int gid) {
 int NAXIS =2;
 double CRPIX[2] =  {N0,N0}; /* reference pixel */
 double PC[2][2] = {{ 1.0,  0.0},
                    { 0.0, 1.0}};
 double CDELT[2] =  {-deltalm, deltalm}; /* pixel widths */
 char CUNIT[2][9] = {"deg", "deg"};
 char CTYPE[2][9] = {"RA---SIN", "DEC--SIN"};
 /* make sure ra0,dec0 is valid */
 if (dec0>90.0) { dec0=90.0; }
 if (dec0<-90.0) { dec0=-90.0; }
 if (ra0>360.0) { ra0=360.0; }
 if (ra0<-360.0) { ra0=-360.0; }
 double CRVAL[2] = {ra0, dec0}; /* phase center */

 double LONPOLE  = 180.0;
 double LATPOLE  = dec0; /* Dec0 of image */
 

 int i,j;
 double *pcij;
 if (!gid) {
  wcs->flag=-1; /* enable memory allocation, only do this once */
 }
 wcsnpv(1); /* wcsnpv() is not thread safe */
 wcsini(1, NAXIS, wcs);
 for (j = 0; j < NAXIS; j++) {
    wcs->crpix[j] = CRPIX[j];
 }

 pcij = wcs->pc;
 for (i = 0; i < NAXIS; i++) {
    for (j = 0; j < NAXIS; j++) {
      *(pcij++) = PC[i][j];
    }
 }

 for (i = 0; i < NAXIS; i++) {
    wcs->cdelt[i] = CDELT[i];
 }

 for (i = 0; i < NAXIS; i++) {
    strcpy(wcs->cunit[i], &CUNIT[i][0]);
 }

 for (i = 0; i < NAXIS; i++) {
    strcpy(wcs->ctype[i], &CTYPE[i][0]);
 }

 for (i = 0; i < NAXIS; i++) {
    wcs->crval[i] = CRVAL[i];
 }

 wcs->pv[0].i = -1;
 wcs->pv[0].m = -1;
 wcs->pv[0].value = -1.0;
 wcs->npv = 0;


 wcs->lonpole = LONPOLE;
 wcs->latpole = LATPOLE;
 wcsset(wcs);

 return 0;
}
