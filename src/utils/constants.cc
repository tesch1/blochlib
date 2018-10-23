

/* constants.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

 /*
 	constants.cc-->several nice constants...
 */


#ifndef __BL_CONSTANTS_CC__
#define __BL_CONSTANTS_CC__ 1

#include "utils/constants.h"
#include <string>

BEGIN_BL_NAMESPACE

const double roundoff = 0.0001;

const double PIx2 = 2.0*PI;             // 2*pi
const double Pi = PI;             // 2*pi
const double pi = PI;             // 2*pi


const double DEG2RAD = PI/180.0;
const double RAD2DEG = 180.0/PI;

const double HZ2RAD = PIx2;
const double RAD2HZ = 1.0/PIx2;

const double HZ2GAUSS = 0.714567e-6;	// h/beta
const double GHZ2GAUSS = 0.714567e3;	// h/beta*1.e9

const double GFREE = 2.00231928;
const double GAMMA1H = 2.67515255e8;
const std::string DEFISO = "1H";

//here 'T'=Telsa
const double No = 6.0221367e23;	// 'items' avagadros number
const double kb = 1.380658e-23;	// bolztman constant (J/K)
const double permVac=12.566370614e-7;		//permetivity of vacume (T^2 m^3/J)
const double hbar = 1.05457266e-34; // h/2Pi (J s)

//if true...to display LU warning if singular matrix is found
// sometimes you wish to use the singular matrix to test for
// something, and do not want the error message everytime
bool LUwarn=true;

END_BL_NAMESPACE
#endif
