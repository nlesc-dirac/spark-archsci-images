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

#ifndef __DATA_H__
#define __DATA_H__
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <complex>
/* CASA stuff */
#include <ms/MeasurementSets/MSIter.h>
#include <tables/Tables/Table.h>
#include <tables/Tables/TableVector.h>
#include <tables/Tables/TableRecord.h>
#include <tables/Tables/TableIter.h>
#include <tables/Tables/ScalarColumn.h>
#include <tables/Tables/ArrayColumn.h>
#include <casa/Arrays/Array.h>
#include <casa/Arrays/Cube.h>


using namespace casa;


//void loadData(char *fname, bool avg=true);
// if imaginarypart==true, data flipped real,imag -> imag,real to extract imaginary part
// if avg, 0: no averaging MFS, 1: average to one channel, 2: channel imaging
// if sort=true, data sorted to improve gridding speed (not working)
extern void 
loadData(Table t, char *fname, int avg=0, bool imaginarypart=false, bool sort=false);
extern void 
freeData();


#endif //__DATA_H__
