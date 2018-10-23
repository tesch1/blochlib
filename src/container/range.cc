 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-25-02
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
	range.cc--> allows for odd acsess to vectors that do not start at one
*/

//examples;;;
//Vector<double> a(6,0), b(3, 3), d(6,0),c;
//c=a(Range(3,6))+b;  //c is a length 3 vector


#ifndef _range_cc_
#define _range_cc_ 1

#include "container/range.h"


BEGIN_BL_NAMESPACE



Range::Range(int f, int l, int s):
	fir_(f), las_(l), str_(s)
{
	if((fir_-las_)%str_ !=0 )
	{
		std::cerr<<std::endl<<"Error: Range() "<<std::endl;
		std::cerr<<" The (first-last) MUST be evenly divisable by the stride"<<std::endl;
		std::cerr<<(*this)<<std::endl;
	}
}


// Operators
Range Range::operator-(int shift) const
{
	return Range(fir_ - shift, las_ - shift, str_);
}

Range Range::operator+(int shift) const
{
	return Range(fir_ + shift, las_ + shift, str_);
}

void Range::operator=(const Range &rhs)
{
	if(&rhs==this)	return;
	fir_=rhs.fir_;
	las_=rhs.las_;
	str_=rhs.str_;
}

Range Range::All(){	return	Range(Start, End); }


std::ostream& operator<<(std::ostream& os, const Range& range)
{
	os << "Range(" << range.first() << "," << range.last() << ","<< range.stride() << ")";

	return os;
}

END_BL_NAMESPACE


#endif



