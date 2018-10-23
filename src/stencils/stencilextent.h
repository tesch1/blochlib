/* stencilextent.h ********/


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

 	stencilbasic.h-> extent of stencils...i.e. how many elements are needed
 	from the current element, and which direction.
 */


#ifndef _stencil_extent_h_
#define _stencil_extent_h_ 1

#include "container/grids/coords.h"
#include "utils/blassert.h"


BEGIN_BL_NAMESPACE


//Stencil Extents...
// these hold the extends of any given stencil for arbitrary dimension
// a 1-D Derivative (of O(h^2)) would have and extent of min=-1 and max=1
//  StencilExtent<1>(-1,1)
// as its stencil is r(i-1) - 2 r(i) + r(i+1)

template<int Dims>
class StencilExtent
{
	private:
		coord<int, Dims> minext_; //min extent of the stencil
		coord<int, Dims> maxext_; //max extent of the stencil

	public:
		StencilExtent():
			minext_(0),
			maxext_(0)
		{}

		StencilExtent(const StencilExtent &cp):
			minext_(cp.minext_),
			maxext_(cp.maxext_)
		{}

		StencilExtent(coord<int, Dims> &mins, coord<int, Dims> &maxs):
			minext_(mins),
			maxext_(maxs)
		{}

		StencilExtent(Vector<int> &mins, Vector<int> &maxs):
			minext_(mins),
			maxext_(maxs)
		{}

		StencilExtent(int min, int max):
			minext_(min),
			maxext_(max)
		{}

		StencilExtent &operator=(const StencilExtent &rhs)
		{
			if(this==&rhs) return *this;
			minext_=rhs.minext_;
			maxext_=rhs.maxext_;
			return *this;
		}

		inline int &min(int i){	return minext_(i);	}
		inline int &max(int i){	return maxext_(i);	}

		inline int min(int i) const {	return minext_(i);	}
		inline int max(int i) const {	return maxext_(i);	}

		coord<int, Dims> &min(){	return minext_;	}
		coord<int, Dims> &max(){	return maxext_;	}

		coord<int, Dims> min() const {	return minext_;	}
		coord<int, Dims> max() const {	return maxext_;	}

		//merges 2 extents together
		void merge(const StencilExtent &lhs)
		{
			for(int i=0;i<Dims;++i)
			{
				minext_[i]=min(minext_[i], lhs.min(i));
				maxext_[i]=max(maxext_[i], lhs.max(i));
			}
		}


		void print(std::ostream &oo)
		{
			oo<<"Stencil Extent: dimension->"<<Dims<<" mins->["<<minext_<<"] maxs->["<<maxext<<"]";
		}

};

template<>
class StencilExtent<1>
{
	private:
		int minext_; //min extent of the stencil
		int maxext_; //max extent of the stencil

	public:
		StencilExtent():
			minext_(0),
			maxext_(0)
		{}

		StencilExtent(const StencilExtent &cp):
			minext_(cp.minext_),
			maxext_(cp.maxext_)
		{}

		StencilExtent(const coord<int, 1> &mins,const coord<int, 1> &maxs):
			minext_(mins(0)),
			maxext_(maxs(0))
		{}

		StencilExtent(Vector<int> &mins, Vector<int> &maxs):
			minext_(mins(0)),
			maxext_(maxs(0))
		{}

		StencilExtent(int min, int max):
			minext_(min),
			maxext_(max)
		{}

		StencilExtent &operator=(const StencilExtent &rhs)
		{
			if(this==&rhs) return *this;
			minext_=rhs.minext_;
			maxext_=rhs.maxext_;
			return *this;
		}

		inline int &min(int i){	return minext_;	}
		inline int &max(int i){	return maxext_;	}

		inline int &min(){	return minext_;	}
		inline int &max(){	return maxext_;	}

		inline int min(int i) const {	return minext_;	}
		inline int max(int i) const {	return maxext_;	}

		inline int min()const{	return minext_;	}
		inline int max()const{	return maxext_;	}

		//merges 2 extents together
		void merge(const StencilExtent<1> &lhs)
		{
			minext_=std::min(minext_, lhs.min());
			maxext_=std::max(maxext_, lhs.max());
		}


		void print(std::ostream &oo)
		{
			oo<<"Stencil Extent: dimension->1 mins->["<<minext_<<"] maxs->["<<maxext_<<"]";
		}

};


template<int Dims>
std::ostream operator<<(std::ostream &oo, StencilExtent<Dims> &out)
{
	out.print(oo);
	return oo;
}

//add 2 different Extents together
template<int Dim1, int Dim2>
StencilExtent<Dim1+Dim2> operator+(const StencilExtent<Dim1> &in1, const StencilExtent<Dim2> &in2)
{
	coord<int, Dim1+Dim2> mins;
	mins.put(Range(Range::start, Dim1), in1.min());
	mins.put(Range(Dim1, Dim2), in2.min());
	coord<int, Dim1+Dim2> maxs;
	maxs.put(Range(Range::start, Dim1), in1.max());
	maxs.put(Range(Dim1, Dim2), in2.max());
	return StencilExtent<Dim1+Dim2>(mins, maxs);
}

//external merge of two Extents (dims must be the same)
// allow the '2 dims' so i can control a better error message...
template<int Dim1,int Dim2>
StencilExtent<Dim1> merge(StencilExtent<Dim1> &in1, StencilExtent<Dim2> &in2)
{
	CompTimeAssert(Dim1==Dim2);
	coord<int, Dim1> mins(in1.min());
	coord<int, Dim1> maxs(in1.max());
	for(int i=0;i<Dim1;++i)
	{
		mins[i]=min(mins[i], lhs.min(i));
		maxs[i]=max(maxs[i], lhs.max(i));
	}
	return StencilExtent<Dim1>(mins, maxs);
}

END_BL_NAMESPACE



#endif





