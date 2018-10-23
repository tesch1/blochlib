/* stencilprep_func.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-30-01
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

 	stencilprep_func.h-> the application of stencils for the 'stencilprep'
 	objects...
 		these operators require a 'stencilprep' object (the grid mapping
 		of the unordered shapes)


 	the author would like to thank "blitz++" for the stencil definitions...(i.e. the
 	stencil maping functions)
 */


#ifndef _stencil_prep_func_h_
#define _stencil_prep_func_h_ 1

#include "stencils/stencilextent.h"
#include "stencils/stencilvector.h"
#include "stencils/stencilprep.h"

#include "container/grids/edgedetect.h"
#include "container/grids/xyzshape.h"

BEGIN_BL_NAMESPACE



//the first declaration of this class is in "stencilvector.h"
//this one is a specialization


//The Stencil Prep Iterator
template<class Shape_t,class NumType_t, int N>
class StencilIterator<StencilPrep<Shape_t,N>,NumType_t, N> :
	public StencilExtent<N>,
	public StencilPrep<Shape_t,N>::iterator
{
	private:


//these are NOT allowed
		StencilIterator(int min_, int max_):
			StencilExtent<N>(min_,max_),
			iterator()
		{}

		StencilIterator(const coord<int,N> &min_, const coord<int,N> &max_):
			StencilExtent<N>(min_,max_),
			iterator()
		{}

		StencilIterator(const StencilExtent<1> &cp):
			StencilExtent<N>(cp),
			iterator()
		{}


	public:
		typedef typename StencilPrep<Shape_t,N>::iterator iterator;
		typedef typename StencilPrep<Shape_t,N>::iterator::DataType DataType;

		StencilIterator():
			StencilExtent<1>(0,0),
			iterator()
		{}

		StencilIterator(int min_, int max_, StencilPrep<Shape_t,N> &inv):
			StencilExtent<N>(min_,max_),
			iterator(inv)

		{}

		StencilIterator(const coord<int,N> &min_, const coord<int,N> &max_, StencilPrep<Shape_t,N> &inv):
			StencilExtent<N>(min_,max_),
			iterator(inv)
		{}

		~StencilIterator(){	}

		int begin(int gben)	 const	{	return abs(min(0));	}
		int end(int gend)	 const	{	return size()-(max(0));	}

		int ubound(int dir)	const	{	return min(dir);	}
		int lbound(int dir)	const	{	return max(dir);	}

		DataType operator()(int i) const {	moveTo(i);	return iterator::operator()();	}

		DataType operator[](int i) const {	moveTo(i);	return iterator::operator()();		}

};

/* The 'DUMMY' Data_t data structure...this allows me to write the global iterators
for the data strcutures that can contain up to 10 different elements
however, if only '3' of the data strucutures are used, this optimizes things so that
the 'dummy' is not even called/used */

class DummyDataIterator;

class DummyData{
	public:
		typedef DummyDataIterator iterator;
		typedef double numtype;
		DummyData(){};
		DummyData(const DummyData &cp){};
};

class DummyDataIterator{
	private:

	public:
		DummyDataIterator(){}
		DummyDataIterator(const DummyData &loo){}
		DummyDataIterator( DummyData *loo){}

		typedef double DataType;
		typedef double numtype;

		int begin(int gben)	 const	{	return 1;	}
		int end(int gend)	 const	{	return 1;	}

		int ubound(int dir)	const	{	return 1;	}
		int lbound(int dir)	const	{	return 1;	}

		void moveTo(int) {	};

		DataType operator()(int i) const {	return 	numtype(0);	}

		DataType operator[](int i) const {	return 	numtype(0);	}
		void operator++()	{	};
		void operator++(int){	}
};


struct ObjectTrait<DummyData> {
    typedef double numtype;
};


//this is the function 'class' that holds each stencil type.
//it is contained in a macro for ease
// The 'data_t' class is the container the was MADE from the grid
// that acts as the acctuall operation base....
// i.e. the stencils will be appiled to the 'DATA' NOT to the grid...



const DummyData _DummyData__;

#define BEGIN_STENCILPREP_STENCIL(NAME) \
template<class Shape_t,class NumType_t, int Dim,class Data_t1>	\
class NAME ## _stprep_cl:	\
	public StencilIterator<StencilPrep<Shape_t,Dim>,NumType_t, Dim> ,	\
	public Data_t1::iterator	\
{	\
	private:	\
		NAME ## _stprep_cl(){}	\
	public:	\
		typedef NumType_t numtype;	\
		typedef StencilPrep<Shape_t,Dim> grids;	\
		typedef Data_t1 container;	\
		typedef StencilIterator<StencilPrep<Shape_t,Dim>,NumType_t, Dim> stencil_iterator;	\
		typedef typename stencil_iterator::DataType StencilData;	\
		typedef typename Data_t1::iterator data_iterator;	\
\
		NAME ## _stprep_cl(StencilPrep<Shape_t,Dim> &in, Data_t1 &ind):	\
			 stencil_iterator(-2,2,in),	\
			 Data_t1::iterator(ind)	\
		{}	\
\
		NAME ## _stprep_cl(StencilPrep<Shape_t,Dim> *in, Data_t1 *ind):	\
			 stencil_iterator(-2,2,in),	\
			 Data_t1::iterator(ind)	\
		{}	\
\
		typename ObjectTrait<NumType_t>::numtype operator()(int i, int dir, int slice)	\
		{	\
			stencil_iterator::moveTo(i);	\
			data_iterator::moveTo(i);	\
			return NAME ## _stprep_func(*this, dir,slice);	\
		}	\
\
		NumType_t operator()(int i, int dir)	\
		{	\
			stencil_iterator::moveTo(i);	\
			data_iterator::moveTo(i);	\
			return NAME ## _stprep_func(*this, dir);	\
		}	\
	\
		NumType_t operator()(int i)	\
		{	\
			stencil_iterator::moveTo(i);	\
			data_iterator::moveTo(i);	\
			return NAME ## _stprep_func(*this);	\
		}	\
	\
		NumType_t operator[](int i)	\
		{	\
			stencil_iterator::moveTo(i);	\
			data_iterator::moveTo(i);	\
			return NAME ## _stprep_func(*this);	\
		}	\
	\
		void moveTo(int i)	\
		{	\
			stencil_iterator::moveTo(i);	\
			data_iterator::moveTo(i);	\
		}	\
	\
		NumType_t nearest(int idx, int dim)	\
		{	\
			int i=stencil_iterator::nearest()(idx,dim);	\
			if(i!=-1){	\
				return data_iterator::operator()(i);	\
			}else{	\
				return NumType_t(0);	\
			}	\
		}	\
\
		NumType_t nextnearest(int idx, int dim)	\
		{	\
			int i=stencil_iterator::nextnearest()(idx,dim);	\
			if(i!=-1){	\
				return data_iterator::operator()(i);	\
			}else{	\
				return NumType_t(0);	\
			}	\
		}	\
	\
		NumType_t current()	\
		{	\
			return data_iterator::operator()();	\
		}	\
	\
		inline int size(){	return std::min(data_iterator::size(), stencil_iterator::size());	}	\
	\
		void operator++()	\
		{	stencil_iterator::operator++();	data_iterator::operator++();	}	\
	\
		void operator++(int)	\
		{	stencil_iterator::operator++();	data_iterator::operator++();	}	\
	\
		operator bool()	\
		{	return stencil_iterator() || data_iterator();	}	\
	\
		inline int curpos() \
		{	return stencil_iterator::curpos();	}	\
};	\


#define BEGIN_STPREP_SINGLE_STENCIL(NAME, WHICH) \
template<class Data_t, class Shape_t, int Dim>	\
typename ObjectTrait<Data_t>::numtype NAME ## _stprep_func(	\
	NAME ## _stprep_cl<Shape_t,typename Data_t::numtype, Dim,Data_t> &WHICH,	\
	int dir)	\



#define END_STPREP_SINGLE_STENCIL(NAME) \
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME(Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs, int dir)	\
{	\
	NAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		if(iter2.curpos()>=lhs.size()-1) return;	\
		iter.moveTo(iter2.curpos());	\
		iter2()=iter(iter2.curpos(), dir);	\
		++iter2;	\
	}	\
}


//this macro makes the 'slice' function version of the stencil

#define BEGIN_STPREP_SINGLE_SLICE_STENCIL(NAME) \
template<class Data_t, class Shape_t, int Dim>	\
typename ObjectTrait<typename Data_t::numtype>::numtype NAME ## _stprep_func(	\
	NAME ## _stprep_cl<Shape_t,typename Data_t::numtype, Dim,Data_t> &A,	\
	int dir,	\
	int slice)	\



#define END_STPREP_SINGLE_SLICE_STENCIL(NAME)	\
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME ## _3D(Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs)	\
{	\
	NAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		if(iter2.curpos()>=lhs.size()-1) return;	\
		iter.moveTo(iter2.curpos());	\
\
		iter2()[0]=iter(iter2.curpos(), 0,0);	\
		iter2()[1]=iter(iter2.curpos(), 1,1);	\
		iter2()[2]=iter(iter2.curpos(), 2,2);	\
\
		++iter2;	\
	}	\
}	\
\
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME (Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs, int dir, int slice)	\
{	\
	NAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		if(iter2.curpos()>=lhs.size()-1) return;	\
		iter.moveTo(iter2.curpos());	\
\
		iter2()[slice]=iter(iter2.curpos(), dir,slice);	\
\
		++iter2;	\
	}	\
}	\



/*********************** CENTER POINT DERIVATIVES*****************/

//first derivative to second order
BEGIN_STENCILPREP_STENCIL(Derivative_1_2)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_1_2,A){
	return (A.nearest(1,dir)-A.nearest(-1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_1_2)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_2){
	return (A.nearest(1,dir)[slice]-A.nearest(-1,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_2)

//first derivative to 2nd order normalized
BEGIN_STENCILPREP_STENCIL(Derivative_1_2n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_1_2n,A){
//	cout<<A.nearest(1,dir)<<" "<<-A.nearest(-1,dir)<<endl;
	return stencil_r2*(A.nearest(1,dir)-A.nearest(-1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_1_2n)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_2n){
	return stencil_r2*(A.nearest(1,dir)[slice]-A.nearest(-1,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_2n)

//first derivative to 4th order
BEGIN_STENCILPREP_STENCIL(Derivative_1_4)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_1_4,A){
	return (A.nearest(-2,dir)-A.nearest(2,dir))+8.0*(A.nearest(1,dir)-A.nearest(-1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_1_4)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_4){
	return  (A.nearest(-2,dir)[slice]-A.nearest(2,dir)[slice])+
			8.0*(A.nearest(1,dir)[slice]-A.nearest(-1,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_4)

//first derivative to 4th order normalized
BEGIN_STENCILPREP_STENCIL(Derivative_1_4n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_1_4n,A){
	return (A.nextnearest(-2,dir)-A.nextnearest(2,dir))+
			8.0*(A.nearest(1,dir)-A.nearest(-1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_1_4n)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_4n){
	return stencil_r12*((A.nextnearest(-2,dir)[slice]-A.nextnearest(2,dir)[slice])+
			8.0*(A.nearest(1,dir)[slice]-A.nearest(-1,dir)[slice]));
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_1_4n)

//second derivatative to 2nd order
BEGIN_STENCILPREP_STENCIL(Derivative_2_2)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_2_2,A){
	return (A.nearest(-1,dir)-A.current()+A.nearest(1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_2_2)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_2){
	return  (A.nearest(-1,dir)[slice]-A.current()[slice]+A.nearest(1,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_2)

//second derivatative to 2nd order Normalized
BEGIN_STENCILPREP_STENCIL(Derivative_2_2n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_2_2n,A){
	return (A.nearest(-1,dir)-A.current()+A.nearest(1,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_2_2n)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_2n){
	return  (A.nearest(-1,dir)[slice]-A.current()[slice]+A.nearest(1,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_2n)

//second derivatative to 4th order
BEGIN_STENCILPREP_STENCIL(Derivative_2_4)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_2_4,A){
	return (-30.0*A.current() + 16.0*(A.nearest(-1,dir)+A.nearest(1,dir))
			-(A.nextnearest(-2, dir)+A.nextnearest(2,dir)));
}
END_STPREP_SINGLE_STENCIL(Derivative_2_4)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_4){
	return  (-30.0*A.current()[slice] + 16.0*(A.nearest(-1,dir)[slice]+A.nearest(1,dir)[slice])
			-(A.nextnearest(-2, dir)[slice]+A.nextnearest(2,dir)[slice]));
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_4)

//second derivatative to 4th order Normalized
BEGIN_STENCILPREP_STENCIL(Derivative_2_4n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_2_4n,A){
	return stencil_r12*(-30.0*A.current() + 16.0*(A.nearest(-1,dir)+A.nearest(1,dir))
			-(A.nextnearest(-2, dir)+A.nextnearest(2,dir)));
}
END_STPREP_SINGLE_STENCIL(Derivative_2_4n)


BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_4n){
	return stencil_r12*(-30.0*A.current()[slice] + 16.0*(A.nearest(-1,dir)[slice]+A.nearest(1,dir)[slice])
			-(A.nextnearest(-2, dir)[slice]+A.nextnearest(2,dir)[slice]));
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_2_4n)


//third derivatative to 2nd order
BEGIN_STENCILPREP_STENCIL(Derivative_3_2)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_3_2,A){
	return 2.0*(A.nearest(-1,dir)-A.nearest(1,dir))+
			(A.nextnearest(-2,dir)-A.nextnearest(2,dir));
}
END_STPREP_SINGLE_STENCIL(Derivative_3_2)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_3_2){
	return  2.0*(A.nearest(-1,dir)[slice]-A.nearest(1,dir)[slice])+
			(A.nextnearest(-2,dir)[slice]-A.nextnearest(2,dir)[slice]);
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_3_2)

//third derivatative to 2nd order NORMALIZED
BEGIN_STENCILPREP_STENCIL(Derivative_3_2n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_3_2n,A){
	return stencil_r2*(2.0*(A.nearest(-1,dir)-A.nearest(1,dir))+
					(A.nextnearest(-2,dir)-A.nextnearest(2,dir)));
}
END_STPREP_SINGLE_STENCIL(Derivative_3_2n)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_3_2n){
	return  stencil_r2*(2.0*(A.nearest(-1,dir)[slice]-A.nearest(1,dir)[slice])+
			(A.nextnearest(-2,dir)[slice]-A.nextnearest(2,dir)[slice]));
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_3_2n)


//FORTH ORDER 3rd DERIVATIVES NOT ALLOWED FOR THE 'STENCIL PREP'...we'd need 'nextnextnearest'


//4th derivatative to 2nd order
BEGIN_STENCILPREP_STENCIL(Derivative_4_2)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_4_2,A){
	return -4.0*(A.nearest(1,dir)+A.nearest(-1,dir))+
			(A.nextnearest(-2,dir)+A.nextnearest(2,dir))+
			6.0*A.current();
}
END_STPREP_SINGLE_STENCIL(Derivative_4_2)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_4_2){
	return -4.0*(A.nearest(1,dir)[slice]+A.nearest(-1,dir)[slice])+
			(A.nextnearest(-2,dir)[slice]+A.nextnearest(2,dir)[slice])+
			6.0*A.current()[slice];
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_4_2)

//4th derivatative to 2nd order
BEGIN_STENCILPREP_STENCIL(Derivative_4_2n)
BEGIN_STPREP_SINGLE_STENCIL(Derivative_4_2n,A){
	return -4.0*(A.nearest(1,dir)+A.nearest(-1,dir))+
			(A.nextnearest(-2,dir)+A.nextnearest(2,dir))+
			6.0*A.current();
}
END_STPREP_SINGLE_STENCIL(Derivative_4_2n)

BEGIN_STPREP_SINGLE_SLICE_STENCIL(Derivative_4_2n){
	return -4.0*(A.nearest(1,dir)[slice]+A.nearest(-1,dir)[slice])+
			(A.nextnearest(-2,dir)[slice]+A.nextnearest(2,dir)[slice])+
			6.0*A.current()[slice];
}
END_STPREP_SINGLE_SLICE_STENCIL(Derivative_4_2n)

//FORTH ORDER 4th DERIVATIVES NOT ALLOWED FOR THE 'STENCIL PREP'...we'd need 'nextnextnearest'

/************************ GRADIENTS *******************/

#define END_STPREP_SINGLE_STENCIL_GRAD(NAME, SNAME) \
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME(Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs, int dir)	\
{	\
	SNAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		if(iter2.curpos()>=lhs.size()) return;	\
		iter.moveTo(iter2.curpos());	\
		iter2()=iter(iter2.curpos(), dir);	\
		++iter2;	\
	}	\
}

BEGIN_STENCILPREP_STENCIL(Gradient_1_2)
BEGIN_STPREP_SINGLE_STENCIL(Gradient_1_2,A){
	if(!A.edge().IsFace(dir)){
		return stencil_r2*(A.nearest(1,dir)-A.nearest(-1,dir));
	}else{
		if(A.edge().IsFace(1,dir)){// && !A.edge().IsFace(-1,dir)){
		//	cout<<"MOO: "<<A.current()<<endl;
			return A.current()-A.nearest(-1,dir);
		}else if(A.edge().IsFace(-1,dir)){// && !A.edge().IsFace(1,dir)){
		//	cout<<"PPP: "<<A.current()<<endl;
			return A.nearest(1,dir)-A.current();
		}else{
			return A.current();
		}
	}
}

END_STPREP_SINGLE_STENCIL_GRAD(Gradient, Gradient_1_2)

BEGIN_STENCILPREP_STENCIL(Gradient_1_2n)
BEGIN_STPREP_SINGLE_STENCIL(Gradient_1_2n,A){
	if(!A.edge().IsFace(dir)){
		return stencil_r2*(A.nearest(1,dir)-A.nearest(-1,dir))/A.dr(dir);
	}else{
		if(A.edge().IsFace(1,dir)){// && !A.edge().IsFace(-1,dir)){
		//	cout<<"MOO: "<<A.current()<<endl;
			return (A.current()-A.nearest(-1,dir))/A.dr(dir);
		}else if(A.edge().IsFace(-1,dir)){// && !A.edge().IsFace(1,dir)){
		//	cout<<"PPP: "<<A.current()<<endl;
			return (A.nearest(1,dir)-A.current())/A.dr(dir);
		}else{
			return A.current();
		}
	}
}

END_STPREP_SINGLE_STENCIL_GRAD(Gradientn, Gradient_1_2n)

//Second Order Gradient
BEGIN_STENCILPREP_STENCIL(Gradient_1_4)
BEGIN_STPREP_SINGLE_STENCIL(Gradient_1_4,A){
	if(!A.edge().IsFace(dir)){
		if(A.IsNearEdge()){
			return stencil_r2*(A.nearest(1,dir)-A.nearest(-1,dir));
		}else{
			return (A.nextnearest(-2,dir)-A.nextnearest(2,dir))+
				8.0*(A.nearest(1,dir)-A.nearest(-1,dir));
		}
	}else{
		if(A.edge().IsFace(1,dir)){
			return A.current()-A.nearest(-1,dir);
		}else if(A.edge().IsFace(-1,dir)){
			return A.nearest(1,dir)-A.current();
		}else{
			return A.current();
		}
	}
}

BEGIN_STENCILPREP_STENCIL(Gradient_1_4n)
BEGIN_STPREP_SINGLE_STENCIL(Gradient_1_4n,A){
	if(!A.edge().IsFace(dir)){
		if(A.IsNearEdge()){
			return stencil_r2*(A.nearest(1,dir)-A.nearest(-1,dir))/(A.dr(dir));
		}else{
			return (A.nextnearest(-2,dir)-A.nextnearest(2,dir))+
				8.0*(A.nearest(1,dir)-A.nearest(-1,dir))/(A.dr(dir));
		}
	}else{
		if(A.edge().IsFace(1,dir)){
			return (A.current()-A.nearest(-1,dir))/(A.dr(dir));
		}else if(A.edge().IsFace(-1,dir)){
			return (A.nearest(1,dir)-A.current())/(A.dr(dir));
		}else{
			return A.current();
		}
	}
}

END_STPREP_SINGLE_STENCIL_GRAD(Gradient_4n, Gradient_1_4n)





/************************* Lapalcians *********************/


/********** 2D Slices...***********************/
#define BEGIN_STPREP_LAPLACE_STENCIL2D(NAME, WHICH, SLICE) \
template<class Data_t, class Shape_t, int Dim>	\
typename ObjectTrait<Data_t>::numtype NAME ## _stprep_func(	\
	NAME ## _stprep_cl<Shape_t,typename Data_t::numtype, Dim,Data_t> &WHICH, int SLICE) \


#define LAPLACE_STENCIL2D(NAME)	\
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME (Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs, int slice)	\
{	\
	NAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		iter.moveTo(iter2.curpos());	\
		\
		iter2()=iter(iter2.curpos(),slice);	\
\
		++iter2;	\
	}	\
}	\

/*simple cubic spine function...
performs the following
here i=x_i, and j=x_(i-1)...but it assumes
that both i and j are integers
i=1 (the 'cur' index), j=0 (the prev index), k=-1 (the want index)
f(want)=(f''(prev)/6*(i-j) * (

*/

template<class Num>
Num spline(Num cur, Num prev, Num d2yi, Num d2yj)
{
	Num b,a;
	a=2.0;
	b=2.0;
	return a*cur+b*prev+((a*a*a-a)*d2yj+(b*b*b-b)*d2yi)/6.0;
}


BEGIN_STENCILPREP_STENCIL(Laplace2D)
BEGIN_STPREP_LAPLACE_STENCIL2D(Laplace2D,A, slice){

	if(slice==0){
		if(A.edge().IsFace(1) || A.edge().IsFace(2)){
			return A.current();
		}
		return   A.nearest(-1,1) + A.nearest(1,1)
				+ A.nearest(-1,2) + A.nearest(1,2)
				-4.0*A.current();
	}else if(slice==1){
		if(A.edge().IsFace(0) || A.edge().IsFace(2)){
			return A.current();
		}
		return   A.nearest(-1,0) + A.nearest(1,0)
				+ A.nearest(-1,2) + A.nearest(1,2)
				-4.0*A.current();
	}else{
		if(A.edge().IsFace(1) || A.edge().IsFace(0)){
			return A.current();
		}
		return  A.nearest(-1,0) + A.nearest(1,0)
				+ A.nearest(-1,1) + A.nearest(1,1)
				-4.0*A.current();
	}
}



LAPLACE_STENCIL2D(Laplace2D)

BEGIN_STENCILPREP_STENCIL(Laplace2Dn)
BEGIN_STPREP_LAPLACE_STENCIL2D(Laplace2Dn,A, slice){
	if(slice==0){
		if(A.edge().IsFace(1) || A.edge().IsFace(2)){
			return A.current();
		}
		return   (A.nearest(-1,1) -2.0*A.current()+ A.nearest(1,1))/(A.dy()*A.dy())
				+ (A.nearest(-1,2) -2.0*A.current()+ A.nearest(1,2))/(A.dz()*A.dz());
	}else if(slice==1){
		if(A.edge().IsFace(0) || A.edge().IsFace(2)){
			return A.current();
		}
		return   (A.nearest(-1,0) -2.0*A.current()+ A.nearest(1,0))/(A.dx()*A.dx())
				+ (A.nearest(-1,2) -2.0*A.current()+ A.nearest(1,2))/(A.dz()*A.dz());
	}else{
		if(A.edge().IsFace(1) || A.edge().IsFace(0)){
			return A.current();
		}
		return  (A.nearest(-1,0) -2.0*A.current()+ A.nearest(1,0))/(A.dx()*A.dx())
				+ (A.nearest(-1,1) -2.0*A.current()+ A.nearest(1,1))/(A.dy()*A.dy());
	}
}

LAPLACE_STENCIL2D(Laplace2Dn)



#define BEGIN_STPREP_LAPLACE_STENCIL(NAME, WHICH) \
template<class Data_t, class Shape_t, int Dim>	\
typename ObjectTrait<Data_t>::numtype NAME ## _stprep_func(	\
	NAME ## _stprep_cl<Shape_t,typename Data_t::numtype, Dim,Data_t> &WHICH) \


#define LAPLACE_STENCIL(NAME)	\
template<class Shape_t, int Dim, class Data_t1, class Data_t2>	\
void NAME (Data_t1 &lhs, StencilPrep<Shape_t,Dim> &in,  Data_t2 &rhs)	\
{	\
	NAME ## _stprep_cl<Shape_t,typename Data_t2::numtype, Dim,Data_t2> iter(in, rhs);	\
	typename Data_t1::iterator iter2(lhs);	\
	while(iter2){	\
		iter.moveTo(iter2.curpos());	\
		\
		iter2()=iter(iter2.curpos());	\
\
		++iter2;	\
	}	\
}	\




//2nd order NOT normalized
BEGIN_STENCILPREP_STENCIL(Laplace3D)
BEGIN_STPREP_LAPLACE_STENCIL(Laplace3D,A){
  	if(A.edge().IsFace()){
		return A.current();
	}else{
		return ( A.nearest(-1,0) + A.nearest(1,0)
      + A.nearest(-1,1) + A.nearest(1,1)
      + A.nearest(-1,2) + A.nearest(1,2))
     -6.0*A.current();
	}
}

LAPLACE_STENCIL(Laplace3D)


//4th order noramlized....
BEGIN_STENCILPREP_STENCIL(Laplace3Dn)
BEGIN_STPREP_LAPLACE_STENCIL(Laplace3Dn,A){
  	if(A.edge().IsFace()){
		return A.current();
	}else{
		return (A.nearest(-1,0) -2.0*A.current()+ A.nearest(1,0))/(A.dx()*A.dx())
      + (A.nearest(-1,1) -2.0*A.current()+ A.nearest(1,1))/(A.dy()*A.dy())
      + (A.nearest(-1,2) -2.0*A.current()+ A.nearest(1,2))/(A.dz()*A.dz());
	}
}

LAPLACE_STENCIL(Laplace3Dn)

//4th order NOT normalized
BEGIN_STENCILPREP_STENCIL(Laplace3D_4)
BEGIN_STPREP_LAPLACE_STENCIL(Laplace3D_4,A){
	if(A.edge().IsFace()){
		return A.current();
	}else if(A.IsNearEdge()){
		return ( A.nearest(-1,0) + A.nearest(1,0)
			  + A.nearest(-1,1) + A.nearest(1,1)
			  + A.nearest(-1,2) + A.nearest(1,2))
			  -6.0*A.current();
	}else{
		return stencil_r12*((-A.nextnearest(-2,0)+16.0*A.nearest(-1,0) -30.0*A.current()+ 16.0*A.nearest(1,0)-A.nextnearest(2,0))
						   +(-A.nextnearest(-2,1)+16.0*A.nearest(-1,1) -30.0*A.current()+ 16.0*A.nearest(1,1)-A.nextnearest(2,1))
						   +(-A.nextnearest(-2,2)+16.0*A.nearest(-1,2) -30.0*A.current()+ 16.0*A.nearest(1,2)-A.nextnearest(2,2)));
	}
}

LAPLACE_STENCIL(Laplace3D_4)


//4th order noramlized....
BEGIN_STENCILPREP_STENCIL(Laplace3D_4n)
BEGIN_STPREP_LAPLACE_STENCIL(Laplace3D_4n,A){
	if(A.edge().IsFace()){
		return A.current();
	}else if(A.IsNearEdge()){
		return (A.nearest(-1,0) -2.0*A.current()+ A.nearest(1,0))/(A.dx()*A.dx())
	  + (A.nearest(-1,1) -2.0*A.current()+ A.nearest(1,1))/(A.dy()*A.dy())
	  + (A.nearest(-1,2) -2.0*A.current()+ A.nearest(1,2))/(A.dz()*A.dz());
	}else{
		return stencil_r12*((-A.nextnearest(-2,0)+16.0*A.nearest(-1,0) -30.0*A.current()+ 16.0*A.nearest(1,0)-A.nextnearest(2,0))/(A.dx()*A.dx())
						   +(-A.nextnearest(-2,1)+16.0*A.nearest(-1,1) -30.0*A.current()+ 16.0*A.nearest(1,1)-A.nextnearest(2,1))/(A.dy()*A.dy())
						   +(-A.nextnearest(-2,2)+16.0*A.nearest(-1,2) -30.0*A.current()+ 16.0*A.nearest(1,2)-A.nextnearest(2,2))/(A.dz()*A.dz()));
	}
}

LAPLACE_STENCIL(Laplace3D_4n)


END_BL_NAMESPACE


#endif



