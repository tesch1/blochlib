

/* bloch_int_offset.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-13-01
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
 	bloch_int_offset.h-->Offsets... this class is templated to take 2 types of
 	data..either double or coord<>s...and also Templated for 'offset generation'
 	meaning that the magnetic fields could change (thus effecting to offset)
 	or a gradient could be applied, thus altering the offset..etc
 */


#ifndef _bloch_interatction_offset_h_
#define _bloch_interatction_offset_h_ 1

#include "utils/constants.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE

class NullBFcalc{
	public:
		typedef double Offset_T;
	};	//the 'static' offset  (no chaning B field)

//predec so you can specailize
template<class BFcalc=NullBFcalc, class Offset_T=typename BFcalc::Offset_T>
class Offset;


template<>
class  Offset< NullBFcalc, double>
{
	private:
		Vector<double>  offsets;
		bool i_on;

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=0;

		Offset():
			offsets(0,0), i_on(true)
		{}

		template<class Params>
		Offset(Params &bc, double offset=0.0):
			offsets(bc.size(), offset), i_on(true)
		{}

		Offset(int &bc, double offset=0.0):
			offsets(bc, offset), i_on(true)
		{}

		Offset(const Vector<double> &offset):
			offsets(offset), i_on(true)
		{}


		~Offset(){}

		inline double &offset(int i){	return offsets(i);	}
		inline double offset(int i) const {	return offsets(i);	}

		inline double &spinOffset(int i){	return offsets(i);	}
		inline double spinOffset(int i) const {	return offsets(i);	}

		inline void on(){	i_on=true;	}
		inline void off(){	i_on=false;	}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		//a 'dummy' function as no pre-calcualtion is required
		template<class PList>
		void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
		{
			return;
		}

		//here we calculate M x B_s (M 'cross' B_s)
		//template<class ParamIter>
		//inline void function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

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
IMPLIMENTATION for 'NullBFcalc' and 'double' offset type
***********/

template<class ParamIter>
inline rmatrix
  Offset<NullBFcalc, double>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0.0);
	if(!i_on) return out;

								out(0,1)=offsets[pars->curpos()];
	out(1,0)= -offsets[pars->curpos()];
	return out;
}

template<class ParamIter>
inline void
  Offset<NullBFcalc, double>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	if(!i_on) return;

	out(0,1)+=offsets[pars->curpos()];
	out(1,0)-=offsets[pars->curpos()];
}


/*
//here we calculate M x B_s (M 'cross' B_s)
template<class Params>
template<class ParamIter>
inline void
  Offset<Params, NullBFcalc, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		dMdt.x()+=offsets[pars->curpos()]*M.y();
		dMdt.y()-=offsets[pars->curpos()]*M.x();
	}
}
*/

//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<>
  Offset<NullBFcalc, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		return coord<>(
			    offsets[pars->curpos()]*M.y(),
			    -offsets[pars->curpos()]*M.x(),
			    ZeroType<double>::zero()
			   );
	}
	return  ZeroType<coord<> >::zero();
}

template<class PList>
inline void
  Offset<NullBFcalc, double>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(i_on)
	{
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			dMdt[myIt->curpos()].x()+=offsets[myIt->curpos()]*M[myIt->curpos()].y();
			dMdt[myIt->curpos()].y()-=offsets[myIt->curpos()]*M[myIt->curpos()].x();
			++myIt;
		}
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class ParamIter>
inline void
  Offset< NullBFcalc, double>::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class PList>
inline void
  Offset< NullBFcalc, double>::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class PList>
rmatrix
  Offset< NullBFcalc, double>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!i_on) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct+2,ct+2)=offsets[myit.curpos()]/myit.gamma();
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
template<class PList>
rmatrix
  Offset< NullBFcalc, double>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!i_on) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct,ct+1)=offsets[myit.curpos()]*M[myit.curpos()].y();
		tmp(ct+1,ct )=-offsets[myit.curpos()]*M[myit.curpos()].x();
		ct+=3;
		++myIt;
	}
	return tmp;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* coord<> No BF calc offsets !*/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
template<>
class  Offset<NullBFcalc,coord<> >
{
	private:
		Vector<coord<> > offsets; //the list of 'precalced' offsets
		bool i_on;

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;

		static const int PreCalc=0;

		Offset():
			i_on(true)
		{}

		template<class Params>
		Offset(Params &bc, coord<> offset=0.0):
			offsets(bc.size(), offset), i_on(true)
		{}

		Offset(int &bc, coord<> offset=0.0):
			offsets(bc, offset), i_on(true)
		{}

		Offset(const Vector<coord<> > &offset):
			offsets(offset), i_on(true)
		{}


		~Offset(){		}

		inline void on()							{	i_on=true;	}
		inline void off()							{	i_on=false;	}

		inline coord<> &offset(int i){	return offsets(i);	}
		inline coord<> offset(int i) const {	return offsets(i);	}

		inline coord<> &spinOffset(int i){	return offsets(i);	}
		inline coord<> spinOffset(int i) const {	return offsets(i);	}


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


template<class ParamIter>
inline rmatrix
  Offset<NullBFcalc, coord<>  >::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0);
								out(0,1)=offsets[pars->curpos()].z();		out(0,2)= -offsets[pars->curpos()].y();
	out(1,0)= -offsets[pars->curpos()].z();		out(1,2)=offsets[pars->curpos()].x();
	out(2,0)=offsets[pars->curpos()].y();									out(2,1)= -offsets[pars->curpos()].x();
	return out;
}

template<class ParamIter>
inline void
  Offset<NullBFcalc, coord<>  >::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
								out(0,1)=offsets[pars->curpos()].z();		out(0,2)= -offsets[pars->curpos()].y();
	out(1,0)= -offsets[pars->curpos()].z();		out(1,2)=offs[pars->curpos()].x();
	out(2,0)=offsets[pars->curpos()].y();									out(2,1)= -offsets[pars->curpos()].x();
}

//Precalculate the Bfield if need be...
template<class PList>
void
  Offset<NullBFcalc, coord<>  >::
preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	return;
}



//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<>
  Offset<NullBFcalc, coord<>  >::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		//cout<<"a=["<<offs[pars->curpos()]<<"]; b=["<<M<<"]; c=["<<cross(offs[pars->curpos()], M)<<"]"<<endl;
		return cross(offsets[pars->curpos()], M);
	}
	return ZeroType<coord<> >::zero();
}

template<class PList>
inline void
  Offset<NullBFcalc, coord<>  >::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(i_on)
	{
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			dMdt[i]=cross(offsets[myIt.curpos()], M[pars->curpos()]);
			++myIt;
		}
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class ParamIter>
inline void
  Offset<NullBFcalc, coord<> >::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class PList>
inline void
  Offset<NullBFcalc, coord<>  >::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class PList>
rmatrix
  Offset<NullBFcalc, coord<>  >::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct, ct)=offsets[myIt.curpos()].x()/myIt.gammaGauss();;
		tmp(ct+1, ct+1)=offsets[myIt.curpos()].y()/myIt.gammaGauss();;
		tmp(ct+2,ct+2)=offsets[myIt.curpos()].z()/myIt.gammaGauss();
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class PList>
rmatrix
  Offset<NullBFcalc, coord<> >::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct, ct)=offsets[myIt.curpos()].z()*M[pars->curpos()].y()-offsets[myIt.curpos()].y()*M[pars->curpos()].z();
		tmp(ct+1, ct+1)=offsets[myIt.curpos()].x()*M[pars->curpos()].z()-offsets[myIt.curpos()].z()*M[pars->curpos()].x();
		tmp(ct+2, ct+2)=offsets[myIt.curpos()].y()*M[pars->curpos()].x()-offsets[myIt.curpos()].x()*M[pars->curpos()].y();
		ct+=3;
		++myIt;
	}
	return tmp;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* GHRADIENT GRID OFFSETS !*/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

template<class GridEngine_t, int BPops>
class  Offset<  ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,
				double>
{
	private:
		ListBlochParams<GradientGrid<GridEngine_t>, BPops, double > *grads_;
		Vector<double>  offsets;
		bool i_on;
		bool grad_apply;

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=0;

		Offset():
			offsets(0,0), i_on(true), grad_apply(true)
		{}

		Offset( ListBlochParams<GradientGrid<GridEngine_t>, BPops, double > &bc, double offset=0.0):
			grads_(&bc), offsets(bc.size(), offset), i_on(true), grad_apply(true)
		{}

		~Offset(){ grads_=NULL;	}

		inline double offset(int i) const
		{
			if(grad_apply){	return offsets(i)+grads_->offset(i);	}
			return offsets(i);
		}


		inline double &spin_offset(int i){	return offsets(i);	}
		inline double spin_offset(int i) const {	return offsets(i);	}

		inline double &spinOffset(int i){	return offsets(i);	}
		inline double spinOffset(int i) const {	return offsets(i);	}

		inline void on(){	grad_apply=true;	}
		inline void off(){	grad_apply=false;	}

		inline void GradOn(){	grad_apply=true;	}
		inline void GradOff(){	grad_apply=false;	}

		inline void gradOn(){	grad_apply=true;	}
		inline void gradOff(){	grad_apply=false;	}


		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM);

		//a 'dummy' function as no pre-calcualtion is required
		template<class PList>
		void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
		{
			return;
		}

		//here we calculate M x B_s (M 'cross' B_s)
		//template<class ParamIter>
		//inline void function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

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
IMPLIMENTATION for 'NullBFcalc' and 'double' offset type
***********/

template<class GridEngine_t, int BPops>
template<class ParamIter>
inline rmatrix
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double >::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0.0);
	if(!i_on) return out;

								out(0,1)=offset(pars->curpos());
	out(1,0)= -offset(pars->curpos());
	return out;
}

template<class GridEngine_t, int BPops>
template<class ParamIter>
inline void
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	if(!i_on) return;

	out(0,1)+=offset(pars->curpos());
	out(1,0)-=offset(pars->curpos());
}


/*
//here we calculate M x B_s (M 'cross' B_s)
template<class Params>
template<class ParamIter>
inline void
  Offset<Params, NullBFcalc, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		dMdt.x()+=offsets[pars->curpos()]*M.y();
		dMdt.y()-=offsets[pars->curpos()]*M.x();
	}
}
*/

//here we calculate M x B_s (M 'cross' B_s)
template<class GridEngine_t, int BPops>
template<class ParamIter>
inline coord<>
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	if(i_on)
	{
		return coord<>(
			    offset(pars->curpos())*M.y(),
			    -offset(pars->curpos())*M.x(),
			    ZeroType<double>::zero()
			   );
	}
	return  ZeroType<coord<> >::zero();
}

template<class GridEngine_t, int BPops>
template<class PList>
inline void
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	if(i_on)
	{
		typename PList::iterator myIt((*pars));
		while(myIt)
		{
			dMdt[myIt->curpos()].x()+=offset(pars->curpos())*M[myIt->curpos()].y();
			dMdt[myIt->curpos()].y()-=offset(pars->curpos())*M[myIt->curpos()].x();
			++myIt;
		}
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class GridEngine_t, int BPops>
template<class ParamIter>
inline void
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class GridEngine_t, int BPops>
template<class PList>
inline void
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class GridEngine_t, int BPops>
template<class PList>
rmatrix
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!i_on) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct+2,ct+2)=offset(pars->curpos())/myit.gamma();
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
template<class GridEngine_t, int BPops>
template<class PList>
rmatrix
  Offset< ListBlochParams<GradientGrid<GridEngine_t>, BPops, double >,double>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!i_on) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		tmp(ct,ct+1)=offset(pars->curpos())*M[myit.curpos()].y();
		tmp(ct+1,ct )=-offset(pars->curpos())*M[myit.curpos()].x();
		ct+=3;
		++myIt;
	}
	return tmp;
}


END_BL_NAMESPACE

#endif


