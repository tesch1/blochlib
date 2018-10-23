

/* constants.h ********/


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
 	constants.h-->several nice constants...
 */

#ifndef __BL_CONSTANTS_H__
#define __BL_CONSTANTS_H__ 1


#include "blochconfig.h"
#include <string>

BEGIN_BL_NAMESPACE


extern const double roundoff;
//#ifndef HUGE
//#define HUGE HUGE_VAL
//#endif


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#ifndef PI2
#define PI2 6.283185307179586476925286766559
#endif

extern const double pi;
extern const double Pi;
extern const double roundoff;

extern const double PIx2;             // 2*pi


extern const double DEG2RAD;
extern const double RAD2DEG;

extern const double HZ2RAD;
extern const double RAD2HZ ;

extern const double HZ2GAUSS;	// h/beta
extern const double GHZ2GAUSS;	// h/beta*1.e9

extern const double GFREE;
extern const double GAMMA1H ;
extern const std::string DEFISO;

//here 'T'=Telsa
extern const double No;	// 'items' avagadros number
extern const double kb;	// bolztman constant (J/K)
extern const double permVac;		//permetivity of vacume (T^2 m^3/J)
extern const double hbar; // h/2Pi (J s)

END_BL_NAMESPACE

#endif
