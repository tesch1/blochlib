/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	endians.h --> converts little endian bit order to big endian bit ordering or vice versa
*/
#ifndef _endians_h_
#define  _endians_h_ 1

#include "utils/utils.h"

BEGIN_BL_NAMESPACE


namespace BL_endians{

union longchars
{
	long value;
	char chars[4];
};

union doublechars
{
	double value;
	char chars[8];
};

union floatchars
{
	float value;
	char chars[4];
};

union shortchars
{
	short value;
	char chars[2];
};

union intchars
{
	int value;
	char chars[2];
};

void ByteSwap(long& i);
void ByteSwap(int& i);
void ByteSwap(double& d);
void ByteSwap(float& d);
void ByteSwap(short& i);
int  AreWeBigEndian();

}; //end namespace BL_endians

END_BL_NAMESPACE


#endif
