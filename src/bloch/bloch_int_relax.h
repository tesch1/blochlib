
/* bloch_int_relax.h ********/


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
 	bloch_int_relax.h-->simply maintains the T1 and T2 relaxation bits...
 */


#ifndef _bloch_interatction_relax_h_
#define _bloch_interatction_relax_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"


BEGIN_BL_NAMESPACE



template< class Offset_T=double>
class  Relax;

template<>
class  Relax<double>{
	friend class Relax<coord<> >;
	private:
		Vector<double> T1s; //list of T1 relaxation rates which if parellel to the 'Bo'
		Vector<double> T2s;	//list of T2 relaxation rates perpendicular to Bo
	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=0;

		Relax(){}

		Relax(int bc, double inT2=double(0),  double inT1=double(0) ):
			T1s(bc,inT1), T2s(bc, inT2)
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		Relax(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, double inT2=double(0),  double inT1=double(0) ):
			T1s(bc.size(),inT1), T2s(bc.size(), inT2)
		{}

		Relax(const Vector<double> &inT2,const Vector<double> &inT1 ):
			T1s(inT1), T2s( inT2)
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		Relax(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, const Vector<double> &inT2,const Vector<double> &inT1 ):
			T1s(inT1), T2s( inT2)
		{}

		~Relax(){}

		inline double &T2(int i)			{	return T2s(i);	}
		inline double &T1(int i)			{	return T1s(i);	}
		inline double T2(int i) const		{	return T2s(i);	}
		inline double T1(int i) const		{	return T1s(i);	}

		inline Vector<double> &T2() {	return T2s;	}
		inline Vector<double> &T1() {	return T1s;	}

		inline Vector<double> T2() const {	return T2s;	}
		inline Vector<double> T1() const {	return T1s;	}

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
IMPLIMENTATION
***********/
template<class ParamIter>
inline rmatrix
  Relax< double>::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix oo(3,3,0);
	oo(0,0)= -T2s(pars->curpos());
								oo(1,1)= -T2s(pars->curpos());
															oo(2,2)= -T1s(pars->curpos());
	return oo;
}

template<class ParamIter>
inline void
  Relax< double>::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	out(0,0)-=T2s[pars->curpos()];
								out(1,1)-=T2s[pars->curpos()];
															out(2,2)-=T1s[pars->curpos()];
}

//here we calculate M x B_s (M 'cross' B_s)

/*
template<class Params>
template<class ParamIter>
inline void
  Relax<Params, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	dMdt.x()-=M.x()*T2s[pars->curpos()];
	dMdt.y()-=M.y()*T2s[pars->curpos()];
	dMdt.z()+=(pars->Mo()-M.z())*T1s[pars->curpos()];
}
*/

//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<>
  Relax< double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	return coord<>(-M.x()*T2s[pars->curpos()],
				   -M.y()*T2s[pars->curpos()],
				   (pars->Mo()-M.z())*T1s[pars->curpos()]
				  );
}

template<class PList>
inline void
  Relax< double>::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	typename PList::iterator myIt((*pars));
	int i=0;
	while(myIt)
	{
		dMdt[myit.curpos()].x()-=M[myit.curpos()].x()*T2s[myit.curpos()];
		dMdt[myit.curpos()].y()-=M[myit.curpos()].y()*T2s[myit.curpos()];
		dMdt[myit.curpos()].z()+=(myit.Mo()-M[myit.curpos()].z())*T1s[myit.curpos()];
		++myIt;
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class ParamIter>
inline void
  Relax< double>::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class PList>
inline void
  Relax< double>::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class PList>
rmatrix
  Relax< double>::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return rmatrix(M.size()*3, M.size()*3,0);
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class PList>
rmatrix
  Relax< double>::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myit((*pars));
	int ct=0;
	while(myit)
	{
		tmp(ct,ct)=M[myit.curpos()].x()*T2s[myit.curpos()];
		tmp(ct+1, ct+1)=M[myit.curpos()].y()*T2s[myit.curpos()];
		tmp(ct+2, ct+2)=(myit.Mo()-M[myit.curpos()].z())*T1s[myit.curpos()];
		ct+=3;
		++myit;
	}
	return tmp;
}

/****************COORD VERSION *********************/
/****************COORD VERSION *********************/
/****************COORD VERSION *********************/

template<>
class  Relax<coord<> > :
	public Relax<double>
{
	private:
		Vector<coord<> > rhat_; //the 'rhat' direction axis where Bo lies
		Vector<double > Momag; //the 'Mo' magnitude along the rhat axis
		Vector<double> theta_; //the theta angel from z-axis to rhat
		Vector<double> phi_; //the phi angle to rotate the x-yplane into rhat

	//this performs the rotation of the relaxation interaction
	// into the (x,y,z) frame from the rhat frame
		inline coord<> makeInt(int i, coord<> &M);

	public:

		Relax(){}

		Relax(int bc, double inT2=double(0),  double inT1=double(0) ):
			Relax<double>(bc, inT2, inT1),
			rhat_(bc, coord<>(0,0,1)),
			Momag(bc, 1),
			theta_(bc, 0), phi_(bc,0)
		{}

		template<class GridEngine_t, int BPops>
		Relax(ListBlochParams<GridEngine_t, BPops, double> &bc, double inT2=double(0),  double inT1=double(0) ):
			Relax<double>(bc, inT2, inT1),
			rhat_(bc.size(), coord<>(0,0,1)),
			Momag(bc.size(), 1),
			theta_(bc.size(), 0), phi_(bc.size(),0)

		{
			//set the rhat vector from the ListBlochPars
			for(int i=0;i<bc.size();++i)
			{
				Momag[i]=bc.Mo(i);
				if(Momag[i]==0) Momag[i]=1;
			}
		}

		template<class GridEngine_t, int BPops>
		Relax(ListBlochParams<GridEngine_t, BPops, coord<> > &bc, double inT2=double(0),  double inT1=double(0) ):
			Relax<double>(bc, inT2, inT1),
			rhat_(bc.size(), coord<>(0,0,1)),
			Momag(bc.size(), 1),
			theta_(bc.size(), 0), phi_(bc.size(),0)
		{
			//set the rhat vector from the ListBlochPars
			for(int i=0;i<bc.size();++i)
			{
				Momag[i]=norm(bc.Mo(i));
				if(Momag[i]==0) Momag[i]=1;
				rhat_[i]=(1.0/Momag[i])*bc.Mo(i);
				theta_[i]=-acos(rhat_[i].z());
				if(rhat_[i].y()==0) phi_[i]=0;
				else phi_[i]=-atan(rhat_[i].x()/rhat_[i].y());
			}
		}

		Relax(const Vector<double> &inT2,const Vector<double> &inT1 ):
			Relax<double>(inT2, inT1),
			rhat_(inT1.size(), coord<>(0,0,1)),
			Momag(inT1.size(), 1),
			theta_(inT1.size(), 0), phi_(inT1.size(),0)

		{}


		~Relax(){}


		inline coord<> Mo(int i) const		{	return rhat_(i);	}

		inline Vector<coord<> > Mo() const {	return rhat_;	}

		inline void setMo(const Vector<coord<> > &in);
		inline void setMo(const coord<>  &in);
		inline void setMo(int i, coord<>  &in);

		template<class GridEngine_t, int BPops>
		inline void setMo(ListBlochParams<GridEngine_t, BPops, double > &in);

		template<class GridEngine_t, int BPops>
		inline void setMo(ListBlochParams<GridEngine_t, BPops, coord<> > &in);


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
IMPLIMENTATION
***********/

//THe Mo setters
inline void  Relax< coord<> >::setMo(const Vector<coord<> > &in)
{
	rhat_.resize(in.size());
	Momag.resize(in.size());
	theta_.resize(in.size());
	phi_.resize(in.size());
	for(int i=0;i<in.size();++i)
	{
		Momag[i]=norm(in[i]);
		if(Momag[i]==0) Momag[i]=1.0;
		rhat_[i]=(1.0/Momag[i]) *in[i];
		theta_[i]=-acos(rhat_[i].z());
		if(rhat_[i].y()==0) phi_[i]=0;
		else phi_[i]=-atan(rhat_[i].x()/rhat_[i].y());
	}
}

inline void  Relax< coord<> >::setMo(const coord<>  &in)
{
	for(int i=0;i<rhat_.size();++i)
	{
		Momag[i]=norm(in);
		if(Momag[i]==0) Momag[i]=1.0;
		rhat_[i]=(1.0/Momag[i])*in;
		theta_[i]=acos(rhat_[i].z());

		if(rhat_[i].y()==0) phi_[i]=0;
		else phi_[i]=-atan(rhat_[i].x()/rhat_[i].y());
	}
}

inline void  Relax< coord<> >::setMo(int i, coord<>  &in)
{
	Momag[i]=norm(in);
	if(Momag[i]==0) Momag[i]=1.0;
	rhat_[i]=(1.0/Momag[i])*in;
	theta_[i]=-acos(rhat_[i].z());
	if(rhat_[i].y()==0) phi_[i]=0;
	else phi_[i]=-atan(rhat_[i].x()/rhat_[i].y());
}

template<class GridEngine_t, int BPops>
inline void  Relax< coord<> >::setMo(ListBlochParams<GridEngine_t, BPops, double > &in)
{
	rhat_.resize(in.size());
	Momag.resize(in.size());
	theta_.resize(in.size());
	phi_.resize(in.size());

	//set the rhat vector from the ListBlochPars
	for(int i=0;i<in.size();++i)
	{
		Momag[i]=in.Mo(i);
		if(Momag[i]==0) Momag[i]=1;
		rhat_[i]=coord<>(0,0,1);
		theta_[i]=0;
		phi_[i]=0;
	}
}

template<class GridEngine_t, int BPops>
inline void  Relax< coord<> >::setMo(ListBlochParams<GridEngine_t, BPops, coord<> > &in)
{
	rhat_.resize(in.size());
	Momag.resize(in.size());
	theta_.resize(in.size());
	phi_.resize(in.size());

	//set the rhat vector from the ListBlochPars
	for(int i=0;i<in.size();++i)
	{
		Momag[i]=norm(in.Mo(i));
		if(Momag[i]==0) Momag[i]=1;
		rhat_[i]=(1.0/Momag[i])*in.Mo(i);
		theta_[i]=-acos(rhat_[i].z());
		if(rhat_[i].y()==0) phi_[i]=0;
		else phi_[i]=-atan(rhat_[i].x()/rhat_[i].y());
	}
}

inline coord<> Relax< coord<> >::makeInt(int i, coord<> &M)
{
	double cp=cos(phi_[i]), sp=sin(phi_[i]);
	double ct=cos(theta_[i]), st=sin(theta_[i]);
	return
	coord<>(
		-T2s[i]*cp*M.x()+sp*(T2s[i]*ct*M.y()+T1s[i]*st*(Momag[i]-M.z())),
		-T2s[i]*cp*ct*M.y()-T2s[i]*M.x()*sp+T1s[i]*cp*st*(M.z()-Momag[i]),
		T1s[i]*ct*(Momag[i]-M.z())-T2s[i]*M.y()*st
	);
}


template<class ParamIter>
inline rmatrix
  Relax< coord<> >::
jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0);
	double cp=cos(phi_[pars->curpos()]), sp=sin(phi_[pars->curpos()]);
	double ct=cos(theta_[pars->curpos()]), st=sin(theta_[pars->curpos()]);

	out(0,0)=-T2s[pars->curpos()]*cp;
	out(0,1)=T2s[pars->curpos()]*ct*sp;
	out(0,2)=-T1s[pars->curpos()]*sp*st;

	out(1,0)=-T2s[pars->curpos()]*sp;
	out(1,1)=-T2s[pars->curpos()]*cp*ct;
	out(1,2)=T1s[pars->curpos()]*cp*st;

	out(2,0)=0.0;
	out(2,1)=-T2s[pars->curpos()]*st;
	out(2,2)=-T1s[pars->curpos()]*ct;

	return out;
}

template<class ParamIter>
inline void
  Relax< coord<> >::
jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	double cp=cos(phi_[pars->curpos()]), sp=sin(phi_[pars->curpos()]);
	double ct=cos(theta_[pars->curpos()]), st=sin(theta_[pars->curpos()]);
	out(0,0)-=T2s[pars->curpos()]*cp;
	out(0,1)+=T2s[pars->curpos()]*ct*sp;
	out(0,2)-=T1s[pars->curpos()]*sp*st;

	out(1,0)-=T2s[pars->curpos()]*sp;
	out(1,1)-=T2s[pars->curpos()]*cp*ct;
	out(1,2)+=T1s[pars->curpos()]*cp*st;

	out(2,1)-=T2s[pars->curpos()]*st;
	out(2,2)-=T1s[pars->curpos()]*ct;
}

//here we calculate M x B_s (M 'cross' B_s)

/*
template<class Params>
template<class ParamIter>
inline void
  Relax<Params, double>::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	dMdt.x()-=M.x()*T2s[pars->curpos()];
	dMdt.y()-=M.y()*T2s[pars->curpos()];
	dMdt.z()+=(pars->Mo()-M.z())*T1s[pars->curpos()];
}
*/

//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<>
  Relax< coord<> >::
function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	//this is the main rotation equation for the interaction
	return makeInt(pars->curpos(), M);
}

template<class PList>
inline void
  Relax< coord<> >::
function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	typename PList::iterator myIt((*pars));
	int i=0;
	while(myIt)
	{
		dMdt[myit.curpos()]= makeInt(myit.curpos(), M[myit.curpos()]);
		++myIt;
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class ParamIter>
inline void
  Relax< coord<> >::
dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class PList>
inline void
  Relax< coord<> >::
dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class PList>
rmatrix
  Relax< coord<> >::
magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	return rmatrix(M.size()*3, M.size()*3,0);
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class PList>
rmatrix
  Relax< coord<> >::
evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	typename PList::iterator myit((*pars));
	int ct=0;
	coord<> Int;
	while(myit)
	{
		Int=makeInt(myit.curpos(), M[myit.curpos()]);
		tmp(ct,ct)=Int.x();
		tmp(ct+1, ct+1)=Int.y();
		tmp(ct+2, ct+2)=Int.z();
		ct+=3;
		++myit;
	}
	return tmp;
}

END_BL_NAMESPACE

#endif


