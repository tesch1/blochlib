


/* bloch_int_demag.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-8-01
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
 	bloch_int_demagmodulated.h-->Demagnitizing Field interactions for the Bloch equations
 */


#ifndef _bloch_interatction_Modulated_demag_field_h_
#define _bloch_interatction_Modulated_demag_field_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "container/containers.h"

BEGIN_BL_NAMESPACE

/* dipole-dipole **Demagnetizing Field** coupling in the classical sence
	requires only the total magenization, the 'gradient' angle
	and the 'time factor' (based typically on td=1/muo*gamma*Mo)

		dm/dt=gamma (M x Bd)
		where Bd_i = 1/td*Del*[ mz-<Mz> + 1/3(M-<M>) ]

		where Del=[3*(s.z)^2 -1]

		where s is the direction of the magnetization modulation

		(i.e. the direction the gradient has been applied)

		This interaction assumes that



*/



class ModulatedDemagField {

	private:

							//once upone every 'Bloch' function call
		bool i_on;	// flag to turn the interaction off

		double td_; // the time constant
		coord<> s_; //the modulation direction
		double delFact;

	public:

		static const int TotalMag=1; //just need the total mag
		static const int TotalVector=0;
		static const int PreCalc=0;

		ModulatedDemagField();

		ModulatedDemagField(double td, coord<> s=coord<>(0,0,1));

		ModulatedDemagField(const ModulatedDemagField &cp);

		ModulatedDemagField &operator=(const ModulatedDemagField &cp);

		~ModulatedDemagField();

		void setTd(double in);
		void setDirection(const coord<> &in);

		double td() const;
		coord<> direction() const;

		void off();
		void on();

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline coord<> function(double t,const coord<>  &M,const  coord<>  &dMdt, ParamIter *pars,const  coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);
};



/*
template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, DimLessDemagField<GridEng_t, ScaleFunc_t> &out);

template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, DimLessDemagField<GridEng_t, NoScaleFunc> &out);
*/

/************* DIPOLE_DIPOLE DEMAGNETIZING FIELD**************************/
//iMPLIMENTATION
/************* DIPOLE_DIPOLE DEMAGNETIZING FIELD**************************/



template<class PList>
inline void ModulatedDemagField::preCalculate(
	double t, Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM)
{}

template<class ParamIter>
inline rmatrix ModulatedDemagField::jacobian(
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	static rmatrix out(3,3,0);
	//static double fact=0.0;
	if(!i_on || !td_) return rmatrix(3,3,0);
		coord<> Bd=delFact/td_/3.0*
					coord<>(totM.x()/pars->TotMo()-M.x(),
			        totM.y()/pars->TotMo()-M.y(),
	                2.0*(M.z()-totM.z())/pars->TotMo()
	               );
			out(0,1)=-Bd.z();		out(0,2)=Bd.y();
	out(1,0)=Bd.z();				out(1,2)=-Bd.x();
	out(2,0)=-Bd.y();		out(2,1)=Bd.x();
	return out;
}

template<class ParamIter>
inline void ModulatedDemagField::jacobian(
	rmatrix &out,
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	if(!i_on|| !td_) return;
		coord<> Bd=delFact/td_/3.0*
					coord<>(totM.x()/pars->TotMo()-M.x(),
			        totM.y()/pars->TotMo()-M.y(),
	                2.0*(M.z()-totM.z())/pars->TotMo()
	               );
							out(0,1)-=Bd.z();		out(0,2)+=Bd.y();
	out(1,0)+=Bd.z();								out(1,2)-=Bd.x();
	out(2,0)-=Bd.y();		out(2,1)+=Bd.x();
}



//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<> ModulatedDemagField::function(
	double t,
	const coord<>  &M,
	const coord<>  &dMdt,
	ParamIter *pars,
	const coord<> &totM)
{
	//  dM/dt=M x (gamma Bd)
	if(i_on && td_)
	{

		coord<> Bd(totM.x()/(pars->TotMo())-M.x(),
			        totM.y()/(pars->TotMo())-M.y(),
	                2.0*(M.z()-totM.z()/(pars->TotMo())));
	    return cross(M, Bd)*(delFact/td_/3.0);
	}
	return ZeroType<coord<> >::zero();
}

template<class PList>
inline void ModulatedDemagField::function(
	double t,
	Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM)
{
	if(i_on && td_)
	{
		coord<> Bd;
		for(int i=0;i<M.size();++i)
		{

		coord<> Bd(totM.x()/pars->TotMo()-M.x(),
			        totM.y()/pars->TotMo()-M.y(),
	                2.0*(M.z()-totM.z()/pars->TotMo()));
			dMdt[i]+=cross(M[i], Bd)*delFact/td_/3.0;
		}
	}
}


std::ostream &operator<<(std::ostream &oo, ModulatedDemagField &out);




END_BL_NAMESPACE

#endif

