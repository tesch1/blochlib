


/* bloch_int_bulks.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-6-01
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
 	bloch_int_bulks.h-->Bulk Susceptibility interactions for the Bloch equations
 */


#ifndef _bloch_interatction_bulk_sus_h_
#define _bloch_interatction_bulk_sus_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"

BEGIN_BL_NAMESPACE


/* this class conatins the lovely pieces for including 'Bulk Susceptibility' into
  the Block equation solver
  	the bulk susceptibility simply adds an extra Mz offset of the form
  		B_s=gamma*mu_o*D*Mz
  	where gamma is the spin gamma factor
  	mu_o--> permitivity in a vacume
  	D--> An important "shape factor" that depends on the sample shape
  		D=0 for a sphere (i.e. no suscep)
  		D=1 for a infinitly long cylinder (a typical NMR tube)
  		D=1/3 for a flat disk
  	This factor depends on a nonlocal integral of the magnitization around the sample shape
  	and is hard to calculate..hense here it is only included as a number
  	Mz--> the magnitization along the Z axis

	this is a simple class and acts more as a function storage class...most of the
	spin information is stored in the BlochParams

Things to note about the format presented below
	1) The only yhting that really matters for the Bloch solver is the 'function' function
	   This function takes in dM/dt, M, t, and a parameter list
	   this function format should be consistant across all similar classes and it should
	   integrate seemlessly with the 'Bloch' classes main function

  */

class BulkSus {
	private:
		double D_;			//shape factor (0--1)
		std::string active_;
		const void DError() const 
		{
			BLEXCEPTION(" The shape factor 'D' cannot be more then 1 or less then 0")
		}
	public:

		static const int TotalMag=1;
		static const int TotalVector=0;
		static const int PreCalc=0;

		BulkSus():
			D_(0), active_(""){}

		BulkSus(double D):
			D_(D), active_("")
		{
			//if(D_>1 || D_<0) DError();
		}

		inline double &D()							{	return D_;	}
		inline void D(const double &in)				{	D_=in;	}

		inline std::string &active()				{	return active_;	}
		inline void active(const std::string &in)	{	active_=in;	}


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


/************** BULK SUSCEPTIBILITY *************/

/* this class conatins the lovely pieces for including 'Bulk Susceptibility' into
  the Block equation solver
  	the bulk susceptibility simply adds an extra Mz offset of the form
  		B_s=gamma*mu_o*D*Mz
  	where gamma is the spin gamma factor
  	mu_o--> permitivity in a vacume
  	D--> An important "shape factor" that depends on the sample shape
  		D=0 for a sphere (i.e. no suscep)
  		D=1 for a infinitly long cylinder (a typical NMR tube)
  		D=1/3 for a flat disk
  	This factor depends on a nonlocal integral of the magnitization around the sample shape
  	and is hard to calculate..hense here it is only included as a number
  	Mz--> the magnitization along the Z axis

	this is a simple class and acts more as a function storage class...most of the
	spin information is stored in the BlochParams

Things to note about the format presented below
	1) The only yhting that really matters for the Bloch solver is the 'function' function
	   This function takes in dM/dt, M, t, and a parameter list
	   this function format should be consistant across all similar classes and it should
	   integrate seemlessly with the 'Bloch' classes main function

  */

template<class ParamIter>
inline rmatrix BulkSus::jacobian(double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static rmatrix out(3,3,0);
	static double fact=0.0;
	if(!D_) return rmatrix(3,3,0);
	fact=D_*pars->gamma()*permVac*totM.z();
	//fact=D_*totM.z();
	out(0,0)=0.0;		out(0,1)=-fact;		out(0,2)=0.0;
	out(1,0)=fact;		out(1,1)=0.0;		out(1,2)=0.0;
	out(2,0)=0.0;		out(2,1)=0.0;		out(2,2)=0.0;
	return out;
}

template<class ParamIter>
inline void BulkSus::jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,   coord<> &totM)
{
	static double fact=0.0;
	if(!D_) return;
	fact=D_*pars->gamma()*permVac*totM.z();
	//fact=D_*totM.z();
						out(0,1)-=fact;
	out(1,0)+=fact;
}

//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<> BulkSus::function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM)
{
	static double fact=0;
	if(!D_) return ZeroType<coord<> >::zero();
	fact=D_*pars->gamma()*permVac*totM.z();
	//fact=D_*totM.z();
	return coord<>(fact*M.y(),
	 		      -fact*M.x(),
	 		      ZeroType<double>::zero()
	 		     );
}

template<class PList>
inline void BulkSus::function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
{
	static double fact=0;
	if(!D_) return;
	typename PList::iterator myIt((*pars));
	int i=0;
	while(myIt)
	{
		fact=D_*myIt.gamma()*permVac*totM.z();
		fact=D_*totM.z();
		//cout<<"t: "<<t<<" "<<fact<<endl;
		dMdt[i].x()+=fact*M[i].y();
		dMdt[i].y()-=fact*M[i].x();
		++i;
		++myIt;
	}
}

//the derivative of the interaction with time..(here it is 0)
template<class ParamIter>
inline void BulkSus::dFdt(double t, coord<>  &M, coord<>  &dFdt, ParamIter *pars, coord<> &totM)
{
	return;
}

template<class PList>
inline void BulkSus::dFdt(double t, Vector<coord<> >  &M, Vector<coord<> >  &dFdt, PList *pars, coord<> &totM)
{
	return;
}

//the magentic field matrix here is simple
// there is only Bz and it on the diagonal...
template<class PList>
rmatrix  BulkSus::magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!D_) return tmp;
	double fact;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		fact=D_*myIt.gamma()*permVac*totM.z();
		//fact=D_*totM.z();
		tmp(ct+2,ct+2)=fact;
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
// there is only Mx and My and they are on the diagonal..
template<class PList>
rmatrix  BulkSus::evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(M.size()*3, M.size()*3,0);
	if(!D_) return tmp;
	double fact;
	typename PList::iterator myIt((*pars));
	int ct=0, i=0;
	while(myIt)
	{
		fact=D_*myIt.gamma()*permVac*totM.z();
		//fact=D_*totM.z();
		tmp(ct,ct)=fact*M[i].y();
		tmp(ct+1, ct+1)=-fact*M[i].x();
		++i;
		ct+=3;
		++myIt;
	}
	return tmp;
}

END_BL_NAMESPACE

#endif



