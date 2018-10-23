


/* bloch_int_raddamp.h ********/


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
 	bloch_int_raddamp.h-->Radiation Damping interactions for the Bloch equations
 */


#ifndef _bloch_interatction_radiation_damping_h_
#define _bloch_interatction_radiation_damping_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "bloch/scalefunc.h"
//#include "stencils/stencilprep_func.h"


BEGIN_BL_NAMESPACE

/* radiation damping
	This functions requires the total magnitization persent in the sample ans simply
	adds the bit
		dMx/dt = -Mx M'z/ (t_r Mo)
		dMy/dt = -My M'z/(t_r Mo)
		dMz/dt = (M'x*Mx+M'y*My)/(t_r Mo)

	where  M'z, M'x, M'y are the TOTAL MAGNITIZATION present
	and t_r is the damping constant
*/

class RadDamp {
	private:
		double t_r;			//shape factor (0--1)

	public:

		static const int TotalMag=1;
		static const int TotalVector=0;
		static const int PreCalc=0;  //no pre calc required

		RadDamp():
			t_r(0){}

		RadDamp(double tr):
			t_r(abs(tr)) {}

		inline void Tr(double tr){	t_r=abs(tr);	}
		inline double Tr()		{	return t_r;		}
		inline void tr(double tr){	t_r=abs(tr);	}
		inline double tr()		{	return t_r;		}

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,  coord<> &totM);
		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		//a 'dummy' function as no pre-calcualtion is required
		template<class PList>
		void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM)
		{
			return;
		}

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline coord<> function(double t,const coord<>  &M,const  coord<>  &dMdt, ParamIter *pars,const  coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};



/*************** RADIATION DAMPING ****************/

/* radiation damping
	This functions requires the total magnitization persent in the sample ans simply
	adds the bit
		dMx/dt = -Mx M'z/ (t_r Mo)
		dMy/dt = -My M'z/(t_r Mo)
		dMz/dt = (M'x*Mx+M'y*My)/(t_r Mo)

	where  M'z, M'x, M'y are the TOTAL MAGNITIZATION present
	and t_r is the damping constant
*/

template<class ParamIter>
inline rmatrix RadDamp::jacobian(double t,coord<>  &M, ParamIter *pars,  coord<> &totM)
{
	static rmatrix out(3,3,0);
	static double fact=0.0;
	if(!t_r) return rmatrix(3,3,0);
	fact=1./(t_r*pars->TotMo() );
										out(0,2)=totM.x()*fact;
										out(1,2)=totM.y()*fact;
	out(2,0)=-totM.x()*fact;		out(2,1)=-totM.y()*fact;
//	cout<<out<<endl;
	return out;
}

template<class ParamIter>
inline void RadDamp::jacobian(
	rmatrix &out,
	double t,
	coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	static double fact=0.0;
	if(!t_r) return;
	fact=1./(t_r*pars->TotMo() );
											out(0,2)+=totM.x()*fact;
											out(1,2)+=totM.y()*fact;
	out(2,0)-=totM.x()*fact;		out(2,1)-=totM.y()*fact;
}


//here we calculate M x B_s (M 'cross' B_s)
template<class ParamIter>
inline coord<> RadDamp::function(
	double t,
	const coord<>  &M,
	const coord<>  &dMdt,
	ParamIter *pars,
	const coord<> &totM)
{
	if(!t_r) return ZeroType<coord<> >::zero();
	double fact=1./(t_r*pars->TotMo());
	//	cout<<endl<<"rad:" <<(M.x()*totM.x()+totM.y()*M.y())<<endl;
	return coord<>(
			-totM.x()*M.z()*fact,
			-totM.y()*M.z()*fact,
			(M.x()*totM.x()+totM.y()*M.y())*fact
		);
	/*return coord<>(
			-totM.x()*M.z()*fact,
			-totM.y()*M.z()*fact,
			(M.x()*totM.x()+totM.y()*M.y())*fact
		);*/
}

template<class PList>
inline void RadDamp::function(
	double t,
	Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM)
{
	if(t_r==0.0) return;
	for(int i=0;i<pars->size();++i)
	{
		//cout<<"t: "<<t<<" "<<fact<<endl;
		double fact=1./(t_r* (pars(i)->TotMo()) );
		dMdt[i].x()-=totM.x()*M[i].z()*fact;
		dMdt[i].y()-=totM.y()*M[i].z()*fact;
		dMdt[i].z()+=(totM.x()*M[i].x()+totM.y()*M[i].y())*fact;
	}
}

//the magentic field matrix here is simply along the diagonal...
template<class PList>
rmatrix  RadDamp::magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(t_r==0.0) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0;
	while(myIt)
	{
		double fact=1./(t_r* (myIt.TotMo()) );
		tmp(ct,ct)=fact*totM.y();
		tmp(ct+1,ct+1)=-fact*totM.x();
		tmp(ct+2, ct+2)=0;
		ct+=3;
		++myIt;
	}
	return tmp;
}

//the evolution matrix (M x B) here is pretty easy
// they are on the diagonal..
template<class PList>
rmatrix  RadDamp::evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM)
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(t_r==0.0) return tmp;
	typename PList::iterator myIt((*pars));
	int ct=0, i=0;
	while(myIt)
	{
		double fact=1./(t_r* (myIt.TotMo()) );
		tmp(ct,ct)=-fact*totM.y()*M[i].z();
		tmp(ct+1,ct+1)=-fact*totM.x()*M[i].z();
		tmp(ct+2, ct+2)=(totM.x()*M[i].x()+totM.y()*M[i].y())*fact;
		ct+=3;
		++i;
		++myIt;
	}
	return tmp;
}


END_BL_NAMESPACE

#endif

