/* xyzshape.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-28-01
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
 	xyzshape.cc-->A grid in the 'uniform' way (i.e. rectangular x,y,z) that can
 	take a plethera of shape functions....these shape functions CHOP down the
 	large base class (unifom grid) into something that the 'XYZshapeFunc' specifies
 	the iterators over this grid then ONLY work over the validated grid points.

 	The compilable bits...Basicially the "XYZfull" grid type...

 */


#ifndef _XYZshape_cc_
#define _XYZshape_cc_ 1


#include "container/grids/xyzshape.h"


BEGIN_BL_NAMESPACE


#define  __Current_Shape XYZshape<XYZfull,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>

void __Current_Shape::init()
{
	data_.resize(Grid<UniformGrid>::size());

	//Grid<UniformGrid>::reset();
	Griditerator myit(*this);
	int ct=0;
	while(myit)
	{
		data_(ct)=(myit.Point());
		++ct;
		++myit;
	}
	data_.resizeAndPreserve(ct);
	if(data_.size()<=0) EmptyWarn();
}


const void __Current_Shape::EmptyWarn()
{
	std::cerr<<std::endl<<"***Warning: The XYZ Shape Grid is EMPTY.***"<<std::endl;
}

__Current_Shape &__Current_Shape::operator=(const __Current_Shape &rhs)
{
	if(this==&rhs) return *this;
	data_=rhs.data_;
	Grid<UniformGrid>::operator=(rhs);
	return *this;
}


void __Current_Shape::PrintGrid(std::ostream &oo)
{
	Griditerator MyIt(*this);
	while(MyIt)
	{
		oo<<MyIt.Point()<<std::endl;
		++MyIt;
	}
}

void __Current_Shape::PrintShape(std::ostream &oo)
{
	//iterator MyIt(*this);
	int i=0;
	while( i<data_.size() )
	{
		oo<<data_(i)<<std::endl;
		++(i);
	}
}

std::ostream &operator<<(std::ostream &oo, __Current_Shape &out)
{
	out.PrintShape(oo);
	return oo;
}




#undef __Current_Shape

END_BL_NAMESPACE


#endif
