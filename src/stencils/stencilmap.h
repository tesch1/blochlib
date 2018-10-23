/* stencilmap.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-29-01
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

 	stencilmap.h-> this class in entirely used as a container to
 	determin the entire stencil layout....
 	 a 1-D Derivative has a stencil like so
 	  r(i)= r(i-1) - 2 r(i) + r(i+1)

	the stencil map would then be (-1,0,1)

	its base class is the 'stencilextent'...the map is acctually stored in
	the template parameters..up to '9' parameters are allowed

	the classes are Macroed to make typing easier

 */


#ifndef _stencil_map_h_
#define _stencil_map_h_ 1

#include "container/grids/coords.h"
#include "stencils/stencilextent.h"

#include "utils/blassert.h"

#ifdef __WIN32__
 #include "blochconfigCW.h"
#else
 #include "blochconfig.h"
#endif

#ifdef HAVE_CLIMITS
 #include <climits>
#else
 #include <limits.h>
#endif


BEGIN_BL_NAMESPACE


template<int p1=INT_MIN,int p2=INT_MIN,int p3=INT_MIN,int p4=INT_MIN,int p5=INT_MIN,
		int p6=INT_MIN,int p7=INT_MIN,int p8=INT_MIN,int p9=INT_MIN,int p10=INT_MIN>
struct StencilMap
{};


template<int p1>
struct StencilMap<p1>
{
		static const int shift1=p1;
};

template<int p1, int p2>
struct StencilMap<p1,p2,p3>
{
		static const int shift1=p1;
		static const int shift2=p2;
};

template<int p1, int p2, int p3>
struct StencilMap<p1,p2,p3>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
};

template<int p1, int p2, int p3, int p4>
struct StencilMap< p1, p2, p3, p4>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;

};

template<int p1, int p2, int p3, int p4, int p5>
struct StencilMap< p1, p2, p3, p4, p5>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;
		static const int shift5=p5;

};

template<int p1, int p2, int p3, int p4, int p5, int p6>
struct StencilMap< p1, p2, p3, p4, p5, p6>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;
		static const int shift5=p5;
		static const int shift6=p6;
};

template<int p1, int p2, int p3, int p4, int p5, int p6>
struct StencilMap<p1,p2,p3,p4,p5, p6,p7>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;
		static const int shift5=p5;
		static const int shift6=p6;
		static const int shift7=p7;
};

template<int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8>
struct StencilMap<p1,p2,p3, p4,p5, p6,p7, p8>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;
		static const int shift5=p5;
		static const int shift6=p6;
		static const int shift7=p7;
		static const int shift8=p8;
};

template<int p1, int p2, int p3, int p4, int p5, int p6, int p7, int p8, int p9>
struct StencilMap< p1, p2, p3, p4, p5, p6, p7, p8, p9>
{
		static const int shift1=p1;
		static const int shift2=p2;
		static const int shift3=p3;
		static const int shift4=p4;
		static const int shift5=p5;
		static const int shift6=p6;
		static const int shift7=p7;
		static const int shift8=p8;
		static const int shift9=p9;
};



END_BL_NAMESPACE



#endif





