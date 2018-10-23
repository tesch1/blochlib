

/* bloch_int_bfield.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-15-01
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
 	bloch_int_bfield.h-->An offset distribution that depends on another 'Bfield'
 	this requires an 'aux' class that determins what the Bfield is....
 	******************************
 	this class acts as the interface between a custom class that has the function

 	'coord<> Bfield(double t, int indx)'

 	the 'indx' will be index from the 'ListBlochParams' HOWEVER, this
 	index will match the index from the Grid Iterator and thus if you use the
 	same Grid for both, this index will be the same...

 	It is ASSUMED that the custom class (i.e the one you write) uses the same grid.

 	the return value above is the (Bx(t), By(t), Bz(t))
 	******************************

 	There are several options available to you to match the situation you have

 	--The 'normal' situation is that you have
 	a) a large static field in the 'nhat' direction (Box, Boy, Boz)..this large field defines the
 	   Larmor fequency of the system..and thus the natural frequecy of the system and
 	   will be 'factored out'
 	b) The OFFSETS of the system will them be calculated as
 	   gammaGauss()*(Bx(t)-Box, By(t)-Boy, Bz(t)-Boz)


	***NOTE::***************
	* you can only use this Interaction IF you choose the
	* 'Offset_T=coord<>' for the ListBlochParams
	************************

 */


#ifndef _bloch_interatction_Bfield_h_
#define _bloch_interatction_Bfield_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "bloch/bloch_int_offset.h"

BEGIN_BL_NAMESPACE

template<class BCalculator>
class  Offset<BCalculator,typename BCalculator::Offset_T >;


/***************************************************************
***************************************************************
************** The double base Offset Claculator ************
***************************************************************
***************************************************************/

template<class BCalculator>
class  Offset<BCalculator,double >
{
	private:
		BCalculator *CalcB;
		double Bo_;			//static field
		Vector<double > offs; //the list of 'precalced' offsets
		bool i_on;

		template<class GridEngine_t, int BPops>
		void Init(ListBlochParams<GridEngine_t, BPops, double> &bc);

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;

		//need to pre-calc the bfield
		//if the Field in 'BCalculator' varies in time...
		//NOTE:: THUS you need to have a variable 'BCalculator::Dynamic=1(dynaimc), 0(not dynamic)'
		static const int PreCalc=BCalculator::Dynamic;

		Offset():
			CalcB(NULL),Bo_(0), i_on(true)
		{}

		template<class GridEngine_t, int BPops>
		Offset(BCalculator &bc,
		       ListBlochParams<GridEngine_t, BPops, double> &lp,
		       double Bo=ZeroType<double >::zero()):
			CalcB(&bc), Bo_(Bo),offs(lp.size()), i_on(true)
		{
			Init(lp);
		}

		template<class GridEngine_t, int BPops>
		Offset(  BCalculator &bc,
		         ListBlochParams<GridEngine_t, BPops,double> &lc,
		         const Vector<double> &initOffs,
		         double Bo=ZeroType<double>::zero()):
			CalcB(&bc), Bo_(Bo),offs(initOffs), i_on(true)
		{}


		~Offset(){	CalcB=NULL;	}

		inline double &Bo()							{	return Bo_;	}
		inline void Bo(const double &in)				{	Bo_=in;	}

		inline void on()							{	i_on=true;	}
		inline void off()							{	i_on=false;	}

		inline double &offset(int i){	return offs(i);	}
		inline double offset(int i) const {	return offs(i);	}

		inline double &spinOffset(int i){	return offset(i);	}
		inline double spinOffset(int i) const {	return offset(i);	}


		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class PList>
		void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		//derivative of the interaction vs time...
		template<class ParamIter>
		inline void dFdt(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

		template<class PList>
		inline void dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);


};
/**********
IMPLIMENTATION
***********/
template<class BCalculator>
template<class GridEngine_t, int BPops>
void
  Offset<BCalculator, double >::
Init(ListBlochParams<GridEngine_t, BPops, double> &bc)
{
	typename ListBlochParams<GridEngine_t, BPops, double>::iterator myIt(bc);
	while(myIt)
	{
		myIt.Bo(CalcB->Bfield(0.0, myIt.curpos()));
		offs[myIt.curpos()]=myIt.gammaGauss()*(myIt.Bo()-Bo_);
		++myIt;
	}
}


template<class BCalculator>
template<class ParamIter>
inline rmatrix
  Offset<BCalculator, double>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0);

								out(0,1)=offs[pars->curpos()];
	out(1,0)= -offs[pars->curpos()];
	return out;
}

template<class BCalculator>
template<class ParamIter>
inline void
  Offset<BCalculator, double >::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	out(0,1)+=offs[pars->curpos()];
	out(1,0)-=offs[pars->curpos()];
}

//Precalculate the Bfield if need be...
template<class BCalculator>
template<class PList>
void
  Offset<BCalculator, double>::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(PreCalc){
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			myIt.Bo(CalcB->Bfield(t, myIt.curpos()));
			offs[myIt.curpos()]=myIt.gammaGauss()*(myIt.Bo()-Bo_);
			++myIt;
		}
		return;
	}
}



//here we calculate M x B_s (M 'cross' B_s)
template<class BCalculator>
template<class ParamIter>
inline coord<>
  Offset<BCalculator, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		//cout<<"a=["<<offs[pars->curpos()]<<"]; b=["<<M<<"]; c=["<<cross(offs[pars->curpos()], M)<<"]"<<endl;
		return coord<>(
			    offs[pars->curpos()]*M.y(),
			    -offs[pars->curpos()]*M.x(),
			    ZeroType<double>::zero()
			   );
	}
	return ZeroType<coord<> >::zero();
}

template<class BCalculator>
template<class PList>
inline void
  Offset<BCalculator, double >::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(i_on)
	{
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			dMdt[myIt->curpos()].x()+=offs[myIt->curpos()]*M[myIt->curpos()].y();
			dMdt[myIt->curpos()].y()-=offs[myIt->curpos()]*M[myIt->curpos()].x();
			++myIt;
		}
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class BCalculator>
template<class ParamIter>
inline void
  Offset<BCalculator, double>::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class BCalculator>
template<class PList>
inline void
  Offset<BCalculator, double>::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class BCalculator>
template<class PList>
rmatrix
  Offset<BCalculator, double>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct+2,ct+2)=pars->Bo();
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class BCalculator>
template<class PList>
rmatrix
  Offset<BCalculator, double>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct+2, ct+2)=offs[myit.curpos()]/myit.gamma();
		ct+=3;
		++myIt;
	}
	return tmp;
}




/***************************************************************
***************************************************************
************** The coord<> base Offset Claculator ************
***************************************************************
***************************************************************/
template<class BCalculator>
class  Offset<BCalculator,coord<> >
{
	private:
		BCalculator *CalcB;
		coord<> Bo_;			//static field
		Vector<coord<> > offs; //the list of 'precalced' offsets
		bool i_on;

		template<class GridEngine_t, int BPops>
		void Init(ListBlochParams<GridEngine_t, BPops, coord<> > &bc);

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;

		//need to pre-calc the bfield
		//if the Field in 'BCalculator' varies in time...
		//NOTE:: THUS you need to have a variable 'BCalculator::Dynamic=1(dynaimc), 0(not dynamic)'
		static const int PreCalc=BCalculator::Dynamic;

		Offset():
			CalcB(NULL),Bo_(0), i_on(true)
		{}

		template<class GridEngine_t, int BPops>
		Offset(BCalculator &bc,
		       ListBlochParams<GridEngine_t, BPops, coord<> > &lp,
		       coord<> Bo=ZeroType<coord<> >::zero()):
			CalcB(&bc), Bo_(Bo),offs(lp.size()), i_on(true)
		{
			Init(lp);
		}

		template<class GridEngine_t, int BPops>
		Offset(  BCalculator &bc,
		         ListBlochParams<GridEngine_t, BPops, coord<> > &lc,
		         const Vector<coord<> > &initOffs,
		         coord<> Bo=ZeroType<coord<> >::zero()):
			CalcB(&bc), Bo_(Bo),offs(initOffs), i_on(true)
		{}


		~Offset(){	CalcB=NULL;	}

		inline coord<> &Bo()							{	return Bo_;	}
		inline void Bo(const double &in)				{	Bo_=in;	}
		inline void Bo(const coord<> &in)				{	Bo_=in;	}

		inline void on()							{	i_on=true;	}
		inline void off()							{	i_on=false;	}

		inline coord<> &offset(int i){	return offs(i);	}
		inline coord<> offset(int i) const {	return offs(i);	}

		inline coord<> &spinOffset(int i){	return offset(i);	}
		inline coord<> spinOffset(int i) const {	return offset(i);	}


		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class PList>
		void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		//derivative of the interaction vs time...
		template<class ParamIter>
		inline void dFdt(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

		template<class PList>
		inline void dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);


};

/**********
IMPLIMENTATION
***********/
template<class BCalculator>
template<class GridEngine_t, int BPops >
void
  Offset<BCalculator, coord<>  >::
Init(ListBlochParams<GridEngine_t, BPops, coord<> > &bc)
{
	typename ListBlochParams<GridEngine_t, BPops,coord<> >::iterator myIt(bc);
	while(myIt)
	{
		myIt.Bo(CalcB->Bfield(0.0, myIt.curpos()));
		offs[myIt.curpos()]=myIt.gammaGauss()*(myIt.Bo()-Bo_);
		++myIt;
	}
}


template<class BCalculator>
template<class ParamIter>
inline rmatrix
  Offset<BCalculator, coord<>  >::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0);
								out(0,1)=offs[pars->curpos()].z();		out(0,2)= -offs[pars->curpos()].y();
	out(1,0)= -offs[pars->curpos()].z();		out(1,2)=offs[pars->curpos()].x();
	out(2,0)=offs[pars->curpos()].y();									out(2,1)= -offs[pars->curpos()].x();
	return out;
}

template<class BCalculator>
template<class ParamIter>
inline void
  Offset<BCalculator, coord<>  >::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
								out(0,1)=offs[pars->curpos()].z();		out(0,2)= -offs[pars->curpos()].y();
	out(1,0)= -offs[pars->curpos()].z();		out(1,2)=offs[pars->curpos()].x();
	out(2,0)=offs[pars->curpos()].y();									out(2,1)= -offs[pars->curpos()].x();
}

//Precalculate the Bfield if need be...
template<class BCalculator>
template<class PList>
void
  Offset<BCalculator, coord<>  >::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(PreCalc){
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			myIt.Bo(CalcB->Bfield(t, myIt.curpos()));
			offs[myIt.curpos()]=myIt.gammaGauss()*(myIt.Bo()-Bo_);
			++myIt;
		}
		return;
	}
}



//here we calculate M x B_s (M 'cross' B_s)
template<class BCalculator>
template<class ParamIter>
inline coord<>
  Offset<BCalculator, coord<>  >::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		//cout<<"a=["<<offs[pars->curpos()]<<"]; b=["<<M<<"]; c=["<<cross(offs[pars->curpos()], M)<<"]"<<endl;
		return cross(offs[pars->curpos()], M);
	}
	return ZeroType<coord<> >::zero();
}

template<class BCalculator>
template<class PList>
inline void
  Offset<BCalculator, coord<>  >::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(i_on)
	{
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			dMdt[i]=cross(offs[myIt.curpos()], M[pars->curpos()]);
			++myIt;
		}
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class BCalculator>
template<class ParamIter>
inline void
  Offset<BCalculator, coord<> >::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class BCalculator>
template<class PList>
inline void
  Offset<BCalculator, coord<>  >::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class BCalculator>
template<class PList>
rmatrix
  Offset<BCalculator, coord<>  >::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct, ct)=pars->Bo().x();
		tmp(ct+1, ct+1)=pars->Bo().y();
		tmp(ct+2,ct+2)=pars->Bo().z();
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class BCalculator>
template<class PList>
rmatrix
  Offset<BCalculator, coord<> >::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct, ct)=offs[myIt.curpos()].z()*M[pars->curpos()].y()-offs[myIt.curpos()].y()*M[pars->curpos()].z();
		tmp(ct+1, ct+1)=offs[myIt.curpos()].x()*M[pars->curpos()].z()-offs[myIt.curpos()].z()*M[pars->curpos()].x();
		tmp(ct+2, ct+2)=offs[myIt.curpos()].y()*M[pars->curpos()].x()-offs[myIt.curpos()].x()*M[pars->curpos()].y();
		ct+=3;
		++myIt;
	}
	return tmp;
}

END_BL_NAMESPACE

#endif



