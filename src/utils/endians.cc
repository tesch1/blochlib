/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-04-02
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
	big_to_lit.cc--> converts between big and small endian byte ordering
*/

#ifndef _endians_cc_
#define  _endians_cc_ 1



#include "utils/endians.h"

BEGIN_BL_NAMESPACE


namespace BL_endians{

void ByteSwap(long& i)
{
	longchars tmpj, tmpi;
	tmpi.value=i;
	tmpj.chars[0] = tmpi.chars[3];
	tmpj.chars[1] = tmpi.chars[2];
	tmpj.chars[2] = tmpi.chars[1];
	tmpj.chars[3] = tmpi.chars[0];
	i = tmpj.value;
}

void ByteSwap(float& i)
{
	floatchars tmpj, tmpi;
	tmpi.value=i;
	tmpj.chars[0] = tmpi.chars[3];
	tmpj.chars[1] = tmpi.chars[2];
	tmpj.chars[2] = tmpi.chars[1];
	tmpj.chars[3] = tmpi.chars[0];
	i = tmpj.value;
}

void ByteSwap(int& i)
{
	longchars tmpj, tmpi;
	tmpi.value=i;
	tmpj.chars[0] = tmpi.chars[3];
	tmpj.chars[1] = tmpi.chars[2];
	tmpj.chars[2] = tmpi.chars[1];
	tmpj.chars[3] = tmpi.chars[0];
	i = tmpj.value;
}

void ByteSwap(short& i)
{
	shortchars tmpj, tmpi;
	tmpi.value=i;
	tmpj.chars[0] = tmpi.chars[1];
	tmpj.chars[1] = tmpi.chars[0];
	i = tmpj.value;
}

void ByteSwap(double& d)
{
	doublechars tmpj, tmpi;
	tmpi.value=d;
	tmpj.chars[0] = tmpi.chars[7];
	tmpj.chars[1] = tmpi.chars[6];
	tmpj.chars[2] = tmpi.chars[5];
	tmpj.chars[3] = tmpi.chars[4];
	tmpj.chars[4] = tmpi.chars[3];
	tmpj.chars[5] = tmpi.chars[2];
	tmpj.chars[6] = tmpi.chars[1];
	tmpj.chars[7] = tmpi.chars[0];
	d = tmpj.value;
}


int AreWeBigEndian()
{
	int x = 1;
	if(*(char *)&x == 1) return 0;
	else                 return 1;
}

}; //end namespace BL_endians


END_BL_NAMESPACE


#endif


