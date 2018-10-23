/* stencilvector.h ********/


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

 	stencilvector.h-> the basic 'application' of stencils
 		for vecortial pieces only...

 		to combine these vector stencils...all one needs to do is
 		'operator' them together
 		A=Derivative_1_2(B)+Derivative_3_2(C)/5

 	the author would like to thank "blitz++" for the stencil definitions...(i.e. the
 	stencil maping functions)
 */


#ifndef _stencil_basic_h_
#define _stencil_basic_h_ 1

#include "container/grids/coords.h"
#include "container/Vector/Vector.h"

#include "utils/blassert.h"

#include "stencils/stencilextent.h"



BEGIN_BL_NAMESPACE


template<class Container, class NumType, int Dim>
class StencilIterator;

//The Stencil Vector Iterator
template<class NumType_t>
class StencilIterator<Vector<NumType_t>,NumType_t, 1> :
	public StencilExtent<1>,
	public Vector<NumType_t>::iterator
{
	private:


//these are NOT allowed
		StencilIterator(int min_, int max_):
			StencilExtent<1>(min_,max_),
			iterator()
		{}

		StencilIterator(const coord<int,1> &min_, const coord<int,1> &max_):
			StencilExtent<1>(cp),
			iterator()
		{}

		StencilIterator(const StencilExtent<1> &cp):
			StencilExtent<1>(cp),
			iterator()
		{}


	public:
		typedef typename Vector<NumType_t>::iterator iterator;

		StencilIterator():
			StencilExtent<1>(0,0),
			iterator()
		{}

		StencilIterator(int min_, int max_, Vector<NumType_t> &inv):
			StencilExtent<1>(min_,max_),
			iterator(inv)

		{}

		StencilIterator(const coord<int,1> &min_, const coord<int,1> &max_, Vector<NumType_t> &inv):
			StencilExtent<1>(min_,max_),
			iterator(inv)
		{}

		~StencilIterator(){		}

		int begin(int gben)	 const	{	return abs(min(0));	}
		int end(int gend)	 const	{	return size()-(max(0));	}

		NumType_t operator()(int i) const {	moveTo(i);	return iterator::operator()(i);	}

		NumType_t operator[](int i) const {	moveTo(i);	return iterator::operator()(i);	}

};




//here we have simple Vectorial derivatives...
//the 'underbelly' of this Macro is where on fills in the correct 'shift'
//sentcil functions...

#define VectorDerivsBEG(NAME, ITEM, MIN, MAX) \
template<class NumType_t>		\
class NAME ## _vec_cl:  	\
	public StencilIterator<Vector<NumType_t>,NumType_t, 1> \
{	\
	private:	\
		NAME ## _vec_cl(){};	\
	public:	\
		typedef NumType_t numtype;		\
		typedef Vector<NumType_t> container; \
\
		NAME ## _vec_cl( Vector<NumType_t> &rhs):	\
			StencilIterator<Vector<NumType_t>, NumType_t, 1>(MIN,MAX,rhs)	\
		{}	\
\
		NumType_t operator()(int i) const	 {	moveTo(i);	return NAME ## _func(*this);	}	\
\
		NumType_t operator[](int i) const	 {	moveTo(i);	return NAME ## _func(*this);	}	\
		\
};	\
\
template<class Container>	\
inline typename Container::numtype NAME ## _func(const Container &ITEM)	\



#define StencilVector(NAME) \
template<class NumType_T>	\
inline VExpr<NAME ## _vec_cl<NumType_T> >	\
	NAME(Vector<NumType_T> &rhs)	\
{	\
	return VExpr<NAME ## _vec_cl<NumType_T> >(NAME ## _vec_cl<NumType_T>(rhs));	\
}	\



const double stencil_r2=0.500000000000000000000000000000000;
const double stencil_r6=0.16666666666666666666666666666667;
const double stencil_r8=0.125000000000000000000000000000000;
const double stencil_r12=0.083333333333333333333333333333333;
const double stencil_r144=0.0069444444444444444444444444444444;

/*********************** FOWARD DERIVATIVES **********************/
// Uses the points 'foward' of the tagged point to calc derivatives
//first derivative to 1st order
VectorDerivsBEG(DerivativeF_1_1, A, 0,1){
	return  A.shift(1)-A.shift(0);
}

//first derivative to 1st order nomalized
VectorDerivsBEG(DerivativeF_1_1n, A, 0,1){
	return  A.shift(1)-A.shift(0);
}

//first derivative to 2nd order
VectorDerivsBEG(DerivativeF_1_2, A, 0,2){
	return  4.0*A.shift(1)-3.0*A.shift(0)-A.shift(2);
}

//first derivative to 2nd order nomalized
VectorDerivsBEG(DerivativeF_1_2n, A, 0,2){
	return  stencil_r2*(4.0*A.shift(1)-3.0*A.shift(0)-A.shift(2));
}

//second derivative to 1st order
VectorDerivsBEG(DerivativeF_2_1, A, 0,2){
	return  A.shift(0)+A.shift(2)-2.0*A.shift(1);
}

//second derivative to 1st order nomalized
VectorDerivsBEG(DerivativeF_2_1n, A, 0,3){
	return  A.shift(0)+A.shift(2)-2.0*A.shift(1);
}

//second derivative to 2nd order
VectorDerivsBEG(DerivativeF_2_2, A, 0,3){
	return  2.0*A.shift(0)+4.0*A.shift(2)-5.0*A.shift(1)-A.shift(3);
}

//second derivative to 2nd order normalized
VectorDerivsBEG(DerivativeF_2_2n, A, 0,3){
	return  2.0*A.shift(0)+4.0*A.shift(2)-5.0*A.shift(1)-A.shift(3);
}

//third derivative to 1st order
VectorDerivsBEG(DerivativeF_3_1, A, 0,3){
	return  3.0*A.shift(1)-3.0*A.shift(2)-A.shift(0)+A.shift(3);
}

//third derivative to 1st order nomalized
VectorDerivsBEG(DerivativeF_3_1n, A, 0,3){
	return  3.0*A.shift(1)-3.0*A.shift(2)-A.shift(0)+A.shift(3);
}

//third derivative to 2nd order
VectorDerivsBEG(DerivativeF_3_2, A, 0,4){
	return  18.0*A.shift(1)-24.0*A.shift(2)-5.0*A.shift(0)+14.0*A.shift(3)-3.0*A.shift(4);
}

//third derivative to 2nd order nomalized
VectorDerivsBEG(DerivativeF_3_2n, A, 0,4){
	return  stencil_r2*(18.0*A.shift(1)-24.0*A.shift(2)-5.0*A.shift(0)+14.0*A.shift(3)-3.0*A.shift(4));
}

//4th derivative to 1st order
VectorDerivsBEG(DerivativeF_4_1, A, 0,4){
	return  6.0*A.shift(2)-4.0*A.shift(1)-A.shift(0)-4.0*A.shift(3)+A.shift(4);
}

//4th derivative to 1st order
VectorDerivsBEG(DerivativeF_4_1n, A, 0,4){
	return  6.0*A.shift(2)-4.0*A.shift(1)-A.shift(0)-4.0*A.shift(3)+A.shift(4);
}

//4th derivative to 2nd order
VectorDerivsBEG(DerivativeF_4_2, A, 0,5){
	return  26.0*A.shift(2)-14.0*A.shift(1)-3.0*A.shift(0)-24.0*A.shift(3)+11.0*A.shift(4)-2.0*A.shift(5);
}

//4th derivative to 2nd order
VectorDerivsBEG(DerivativeF_4_2n, A, 0,5){
	return  26.0*A.shift(2)-14.0*A.shift(1)-3.0*A.shift(0)-24.0*A.shift(3)+11.0*A.shift(4)-2.0*A.shift(5);
}


/*********************** BACKWARDS POINT DERIVATIVES *************/
//first derivative to 1st order
VectorDerivsBEG(DerivativeB_1_1, A, -1,0){
	return  A.shift(0)-A.shift(-1);
}

//first derivative to 1st order nomalized
VectorDerivsBEG(DerivativeB_1_1n, A, -1,0){
	return  A.shift(0)-A.shift(-1);
}

//first derivative to 2nd order
VectorDerivsBEG(DerivativeB_1_2, A, -2,0){
	return  A.shift(-2)+3.0*A.shift(0)-4.0*A.shift(-1);
}

//first derivative to 2nd order nomalized
VectorDerivsBEG(DerivativeB_1_2n, A, -2,0){
	return  stencil_r2*(A.shift(-2)+3.0*A.shift(0)-4.0*A.shift(-1));
}

//second derivative to 1st order
VectorDerivsBEG(DerivativeB_2_1, A, -2,0){
	return  A.shift(0)+A.shift(-2)-2.0*A.shift(-1);
}

//second derivative to 1st order nomalized
VectorDerivsBEG(DerivativeB_2_1n, A, -2,0){
	return  A.shift(0)+A.shift(-2)-2.0*A.shift(-1);
}

//second derivative to 2nd order
VectorDerivsBEG(DerivativeB_2_2, A, -3,0){
	return  2.0*A.shift(0)+4.0*A.shift(-2)-5.0*A.shift(-1)-A.shift(-3);
}

//second derivative to 2nd order normalized
VectorDerivsBEG(DerivativeB_2_2n, A, -3,0){
	return  2.0*A.shift(0)+4.0*A.shift(-2)-5.0*A.shift(-1)-A.shift(-3);
}

//third derivative to 1st order
VectorDerivsBEG(DerivativeB_3_1, A, -3,0){
	return  3.0*A.shift(-2)+A.shift(-1)+A.shift(0)-A.shift(-3);
}

//third derivative to 1st order nomalized
VectorDerivsBEG(DerivativeB_3_1n, A, -3,0){
	return  3.0*A.shift(-2)+A.shift(-1)+A.shift(0)-A.shift(-3);
}

//third derivative to 2nd order
VectorDerivsBEG(DerivativeB_3_2, A, -4,0){
	return  5.*A() - 18.*A.shift(-1) + 24.*A.shift(-2) -14.*A.shift(-3)
    + 3.*A.shift(-4);
}

//third derivative to 2nd order nomalized
VectorDerivsBEG(DerivativeB_3_2n, A, -4,0){
	return  stencil_r2*(5.*A() - 18.*A.shift(-1) + 24.*A.shift(-2) -14.*A.shift(-3)
    + 3.*A.shift(-4));
}

//4th derivative to 1st order
VectorDerivsBEG(DerivativeB_4_1, A, -4,0){
	return  6.0*A.shift(2)-4.0*A.shift(-1)-A.shift(0)-4.0*A.shift(-3)+A.shift(-4);
}

//4th derivative to 1st order
VectorDerivsBEG(DerivativeB_4_1n, A, -4,0){
	return  6.0*A.shift(2)-4.0*A.shift(-1)-A.shift(0)-4.0*A.shift(-3)+A.shift(-4);
}

//4th derivative to 2nd order
VectorDerivsBEG(DerivativeB_4_2, A, -5,0){
	return  26.0*A.shift(-2)-14.0*A.shift(-1)-3.0*A.shift(0)-24.0*A.shift(-3)+11.0*A.shift(-4)-2.0*A.shift(-5);
}

//4th derivative to 2nd order
VectorDerivsBEG(DerivativeB_4_2n, A, -5,0){
	return  26.0*A.shift(-2)-14.0*A.shift(-1)-3.0*A.shift(0)-24.0*A.shift(-3)+11.0*A.shift(-4)-2.0*A.shift(-5);
}

/*********************** CENTER POINT DERIVATIVES*****************/
//first derivative to second order
VectorDerivsBEG(Derivative_1_2, A, -1,1){
	return  A.shift(1)-A.shift(-1);
}

//first derivative to 2nd order normalized
VectorDerivsBEG(Derivative_1_2n, A, -1,1){
	return  stencil_r2*(A.shift(1)-A.shift(-1));
}

//first derivative to 4th order
VectorDerivsBEG(Derivative_1_4, A, -2,2){
	return  A.shift(-2)-8.0*A.shift(-1)+8.0*A.shift(1)-A.shift(2);
}

//first derivative to 4th order normlaized
VectorDerivsBEG(Derivative_1_4n, A, -2,2){
	return  stencil_r12*(A.shift(-2)-8.0*A.shift(-1)+8.0*A.shift(1)-A.shift(2));
}

//second derivatative to 2nd order
VectorDerivsBEG(Derivative_2_2, A, -1,1){
	return  A.shift(-1)-2.0*A.shift(0)+A.shift(1);
}

//second derivatative to 2nd order normalized
VectorDerivsBEG(Derivative_2_2n, A, -1,1){
	return  A.shift(-1)-2.0*A.shift(0)+A.shift(1);
}

//second derivatative to 4th order
VectorDerivsBEG(Derivative_2_4, A, -2,2){
	return  -30.*A.shift(0) + 16.*(A.shift(-1)+A.shift(1))
    - (A.shift(-2)+A.shift(2));
}

//second derivatative to 4th order normalized
VectorDerivsBEG(Derivative_2_4n, A, -2,2){
	return  stencil_r12*(-30.*A.shift(0) + 16.*(A.shift(-1)+A.shift(1))
    - (A.shift(-2)+A.shift(2)));
}

//third derivatative to 2nd order
VectorDerivsBEG(Derivative_3_2, A, -2,2){
	return  2.0*A.shift(-1)-1.0*A.shift(-2)-2.0*A.shift(1)+A.shift(2);
}

//third derivatative to 2nd order
VectorDerivsBEG(Derivative_3_2n, A, -2,2){
	return  stencil_r2*(2.0*A.shift(-1)-1.0*A.shift(-2)-2.0*A.shift(1)+A.shift(2));
}

//third derivatative to 4th order
VectorDerivsBEG(Derivative_3_4, A, -3,3){
	return   A.shift(-3) - 8.*A.shift(-2) +13.*A.shift(-1)
     -13.*A.shift(1)+8.*A.shift(2)-A.shift(3);
}

//third derivatative to 4th order normalized
VectorDerivsBEG(Derivative_3_4n, A, -3,3){
	return  stencil_r8*(A.shift(-3) - 8.*A.shift(-2) +13.*A.shift(-1)
     -13.*A.shift(1)+8.*A.shift(2)-A.shift(3));
}

//forth derivatative to 2nd order
VectorDerivsBEG(Derivative_4_2, A, -2,2){
	return  A.shift(-2)-4.0*A.shift(-1)+6.0*A.shift(0)-4.0*A.shift(1)+A.shift(2);
}

//forth derivatative to 2nd order normalized
VectorDerivsBEG(Derivative_4_2n, A, -2,2){
	return  A.shift(-2)-4.0*A.shift(-1)+6.0*A.shift(0)-4.0*A.shift(1)+A.shift(2);
}

//forth derivatative to 4th order
VectorDerivsBEG(Derivative_4_4, A, -3,3){
	return  -1.*A.shift(-3)+12.*A.shift(-2)-39.*A.shift(-1)
    +56.*A-39.*A.shift(1)+12.*A.shift(2)-A.shift(3);
}

//forth derivatative to 4th order normalized
VectorDerivsBEG(Derivative_4_4n, A, -3,3){
	return  stencil_r6*(-1.*A.shift(-3)+12.*A.shift(-2)-39.*A.shift(-1)
    +56.*A-39.*A.shift(1)+12.*A.shift(2)-A.shift(3));
}

StencilVector(Derivative_1_2n);
StencilVector(Derivative_1_4n);
StencilVector(Derivative_1_2);
StencilVector(Derivative_1_4);
StencilVector(Derivative_2_2n);
StencilVector(Derivative_2_4n);
StencilVector(Derivative_2_2);
StencilVector(Derivative_2_4);
StencilVector(Derivative_3_2n);
StencilVector(Derivative_3_4n);
StencilVector(Derivative_3_2);
StencilVector(Derivative_3_4);
StencilVector(Derivative_4_2n);
StencilVector(Derivative_4_4n);
StencilVector(Derivative_4_2);
StencilVector(Derivative_4_4);

StencilVector(DerivativeF_1_1n);
StencilVector(DerivativeF_1_2n);
StencilVector(DerivativeF_2_1n);
StencilVector(DerivativeF_2_2n);
StencilVector(DerivativeF_3_1n);
StencilVector(DerivativeF_3_2n);
StencilVector(DerivativeF_4_1n);
StencilVector(DerivativeF_4_2n);
StencilVector(DerivativeF_1_1);
StencilVector(DerivativeF_1_2);
StencilVector(DerivativeF_2_1);
StencilVector(DerivativeF_2_2);
StencilVector(DerivativeF_3_1);
StencilVector(DerivativeF_3_2);
StencilVector(DerivativeF_4_1);
StencilVector(DerivativeF_4_2);


StencilVector(DerivativeB_1_1n);
StencilVector(DerivativeB_1_2n);
StencilVector(DerivativeB_2_1n);
StencilVector(DerivativeB_2_2n);
StencilVector(DerivativeB_3_1n);
StencilVector(DerivativeB_3_2n);
StencilVector(DerivativeB_4_1n);
StencilVector(DerivativeB_4_2n);
StencilVector(DerivativeB_1_1);
StencilVector(DerivativeB_1_2);
StencilVector(DerivativeB_2_1);
StencilVector(DerivativeB_2_2);
StencilVector(DerivativeB_3_1);
StencilVector(DerivativeB_3_2);
StencilVector(DerivativeB_4_1);
StencilVector(DerivativeB_4_2);



//the application function...for 1 input...for Foward and Backward derivatives
/*  OLD (but perhaps still usefull)
#define Vector1DsFB(NAME,CLASS) \
\
template<class NumType_t1, class NumType_t2> \
void NAME(Vector<NumType_t1> &lhs, Vector<NumType_t2> &rhs) \
{ \
	static const StencilEng<CLASS<NumType_t2>, NumType_t2  > dr;			\
	typename Vector<NumType_t1>::iterator myit(lhs);		\
	typename Vector<NumType_t2>::iterator myit2(rhs);		\
	myit.begin(abs(dr.min(0)));		\
	myit2.begin(abs(dr.min(0)));		\
	myit.end(lhs.size()-dr.max(0));		\
	myit2.end(lhs.size()-dr.max(0));		\
\
	while(myit)	\
	{	\
		myit()=dr.apply(myit2);	\
		++myit; ++myit2;		\
	}	\
}		\
\
\
template<class NumType_t1, class NumType_t2, class NumType_t3> \
void NAME(Vector<NumType_t1> &lhs, Vector<NumType_t2> &rhs, const NumType_t3 &fact) \
{ \
	static const StencilEng<CLASS<NumType_t2>, NumType_t2 > dr;			\
	typename Vector<NumType_t1>::iterator myit(lhs);		\
	typename Vector<NumType_t2>::iterator myit2(rhs);		\
	myit.begin(abs(dr.min(0)));		\
	myit2.begin(abs(dr.min(0)));		\
	myit.end(lhs.size()-dr.max(0));		\
	myit2.end(lhs.size()-dr.max(0));		\
\
	while(myit)	\
	{	\
		myit()=dr.apply(myit2)*fact;	\
		++myit; ++myit2;		\
	}	\
}		\


//macros for CENTRAL stencil functions...uses the foward and backward
// stencils for begining and end points

#define Vector1Ds(NAME,CLASS, FF, BB) \
\
template<class NumType_t1, class NumType_t2> \
void NAME(Vector<NumType_t1> &lhs, Vector<NumType_t2> &rhs) \
{ \
	static const StencilEng<CLASS<NumType_t2>, NumType_t2  >dr;			\
	static const StencilEng<FF<NumType_t2>, NumType_t2 > drf;			\
	static const StencilEng<BB<NumType_t2>, NumType_t2 > drb;			\
	typename Vector<NumType_t1>::iterator myit(lhs);		\
	typename Vector<NumType_t2>::iterator myit2(rhs);		\
	if(myit2.size()>drf.max(0))		\
	{	\
		myit.end(abs(dr.min(0))+1);	\
		myit2.end(abs(dr.min(0))+1);	\
		while(myit)	\
		{	\
			cout<<myit()<<endl;	\
			myit()=drf.apply(myit2);	\
			cout<<myit()<<endl;	\
			++myit; ++myit2;		\
		}	\
	}	\
	myit.reset();	myit2.reset();	\
	if(myit2.size()>abs(drb.min(0))	)	\
	{	\
		myit.begin(myit.size()-dr.max(0)-1);	\
		myit2.begin(myit2.size()-dr.max(0)-1);	\
		myit.end(myit.size());	\
		myit2.end(myit2.size());	\
		while(myit)	\
		{	\
			cout<<myit()<<endl;	\
			myit()=drb.apply(myit2);	\
			cout<<myit()<<endl;	\
			++myit; ++myit2;		\
		}	\
	}	\
	myit.reset();	myit2.reset();	\
	myit.begin(abs(dr.min(0)));		\
	myit2.begin(abs(dr.min(0)));		\
	myit.end(lhs.size()-dr.max(0));		\
	myit2.end(lhs.size()-dr.max(0));		\
\
	while(myit)	\
	{	\
		myit()=dr.apply(myit2);	\
		++myit; ++myit2;		\
	}	\
}		\
\
\
template<class NumType_t1, class NumType_t2, class NumType_t3> \
void NAME(Vector<NumType_t1> &lhs, Vector<NumType_t2> &rhs, const NumType_t3 &fact) \
{ \
	static const StencilEng<CLASS<NumType_t2>, NumType_t2  >dr;			\
	typename Vector<NumType_t1>::iterator myit(lhs);		\
	typename Vector<NumType_t2>::iterator myit2(rhs);		\
	myit.begin(abs(dr.min(0)));		\
	myit2.begin(abs(dr.min(0)));		\
	myit.end(lhs.size()-dr.max(0));		\
	myit2.end(lhs.size()-dr.max(0));		\
\
	while(myit)	\
	{	\
		myit()=dr.apply(myit2)*fact;	\
		++myit; ++myit2;		\
	}	\
}		\

*/
//foward stencil decs
/*
Vector1DsFB(DerivativeF_1_1,	DerivativeF_1_1_cl);
Vector1DsFB(DerivativeF_1_1n,	DerivativeF_1_1n_cl);
Vector1DsFB(DerivativeF_1_2,	DerivativeF_1_2_cl);
Vector1DsFB(DerivativeF_1_2n,	DerivativeF_1_2n_cl);
Vector1DsFB(DerivativeF_2_1,	DerivativeF_2_1_cl);
Vector1DsFB(DerivativeF_2_1n,	DerivativeF_2_1n_cl);
Vector1DsFB(DerivativeF_2_2,	DerivativeF_2_2_cl);
Vector1DsFB(DerivativeF_2_2n,	DerivativeF_2_2n_cl);
Vector1DsFB(DerivativeF_3_1,	DerivativeF_3_1_cl);
Vector1DsFB(DerivativeF_3_1n,	DerivativeF_3_1n_cl);
Vector1DsFB(DerivativeF_3_2,	DerivativeF_3_2_cl);
Vector1DsFB(DerivativeF_3_2n,	DerivativeF_3_2n_cl);
Vector1DsFB(DerivativeF_4_1,	DerivativeF_4_1_cl);
Vector1DsFB(DerivativeF_4_1n,	DerivativeF_4_1n_cl);
Vector1DsFB(DerivativeF_4_2,	DerivativeF_4_2_cl);
Vector1DsFB(DerivativeF_4_2n,	DerivativeF_4_2n_cl);

//backwards stencil decs

Vector1DsFB(DerivativeB_1_1,	DerivativeB_1_1_cl);
Vector1DsFB(DerivativeB_1_1n,	DerivativeB_1_1n_cl);
Vector1DsFB(DerivativeB_1_2,	DerivativeB_1_2_cl);
Vector1DsFB(DerivativeB_1_2n,	DerivativeB_1_2n_cl);
Vector1DsFB(DerivativeB_2_1,	DerivativeB_2_1_cl);
Vector1DsFB(DerivativeB_2_1n,	DerivativeB_2_1n_cl);
Vector1DsFB(DerivativeB_2_2,	DerivativeB_2_2_cl);
Vector1DsFB(DerivativeB_2_2n,	DerivativeB_2_2n_cl);
Vector1DsFB(DerivativeB_3_1,	DerivativeB_3_1_cl);
Vector1DsFB(DerivativeB_3_1n,	DerivativeB_3_1n_cl);
Vector1DsFB(DerivativeB_3_2,	DerivativeB_3_2_cl);
Vector1DsFB(DerivativeB_3_2n,	DerivativeB_3_2n_cl);
Vector1DsFB(DerivativeB_4_1,	DerivativeB_4_1_cl);
Vector1DsFB(DerivativeB_4_1n,	DerivativeB_4_1n_cl);
Vector1DsFB(DerivativeB_4_2,	DerivativeB_4_2_cl);
Vector1DsFB(DerivativeB_4_2n,	DerivativeB_4_2n_cl);

Vector1Ds(Derivative_1_2,	Derivative_1_2_cl,	DerivativeF_1_2_cl,	DerivativeB_1_2_cl);
Vector1Ds(Derivative_1_2n,	Derivative_1_2n_cl,	DerivativeF_1_2n_cl,	DerivativeB_1_2n_cl);
Vector1Ds(Derivative_1_4,	Derivative_1_4_cl,	DerivativeF_1_2_cl,	DerivativeB_1_2_cl);
Vector1Ds(Derivative_1_4n,	Derivative_1_4n_cl,	DerivativeF_1_2n_cl,	DerivativeB_1_2n_cl);
Vector1Ds(Derivative_2_2,	Derivative_2_2_cl,	DerivativeF_2_2_cl,	DerivativeB_2_2_cl);
Vector1Ds(Derivative_2_2n,	Derivative_2_2n_cl,	DerivativeF_2_2n_cl,	DerivativeB_2_2n_cl);
Vector1Ds(Derivative_2_4,	Derivative_2_4_cl,	DerivativeF_2_2_cl,	DerivativeB_2_2_cl);
Vector1Ds(Derivative_2_4n,	Derivative_2_4n_cl,	DerivativeF_2_2n_cl,	DerivativeB_2_2n_cl);
Vector1Ds(Derivative_3_2,	Derivative_3_2_cl,	DerivativeF_3_2_cl,	DerivativeB_3_2_cl);
Vector1Ds(Derivative_3_2n,	Derivative_3_2n_cl,	DerivativeF_3_2n_cl,	DerivativeB_3_2n_cl);
Vector1Ds(Derivative_3_4,	Derivative_3_4_cl,	DerivativeF_3_2_cl,	DerivativeB_3_2_cl);
Vector1Ds(Derivative_3_4n,	Derivative_3_4n_cl,	DerivativeF_3_2n_cl,	DerivativeB_3_2n_cl);
Vector1Ds(Derivative_4_2,	Derivative_4_2_cl,	DerivativeF_4_2_cl,	DerivativeB_4_2_cl);
Vector1Ds(Derivative_4_2n,	Derivative_4_2n_cl,	DerivativeF_4_2n_cl,	DerivativeB_4_2n_cl);
Vector1Ds(Derivative_4_4,	Derivative_4_4_cl,	DerivativeF_4_2_cl,	DerivativeB_4_2_cl);
Vector1Ds(Derivative_4_4n,	Derivative_4_4n_cl,	DerivativeF_4_2n_cl,	DerivativeB_4_2n_cl);
*/

END_BL_NAMESPACE


#endif





