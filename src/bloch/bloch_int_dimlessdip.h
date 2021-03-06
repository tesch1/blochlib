

/* bloch_int_dipdip.h ********/


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
 	bloch_int_dimlessdip.h-->DIMENSIONLESS Dipole-Dipole interactions
 	for the Bloch equations..

 	allows you to set the maximal dipole allowed based on the minimum distance

 	you can set next-nearest neighbors to '1' then all the other
 	dipol-dipole interaction will be smaller....(falling off cubically).. it
 	maintains the possible angular scaling as well...
 */


#ifndef _bloch_interatction_dimless_dipole_dipole_h_
#define _bloch_interatction_dimless_dipole_dipole_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"

BEGIN_BL_NAMESPACE

/* dipole-dipole  coupling in the classical sence
	This functions requires ALL the locations of the spins and thus the GRID!!

		dm/dt= 1.0* (M x Bd)
		where Bd_i = FACTOR*(dr^3/|ri-rj|^3) Sum_j[ (1-3cos(th_ri_rj))/2
		              * [3Mz(rj)*zhat-Mx xhat - My yhat - Mz zhat]]


	FACTOR--> is the parameter you specify
	dr--> the distance to the next nearest neighbor
	Bd_i is the effective magnetic field at point 'ri' generated by all the little
	 dipoles around it.  Each Bd_i thus must be calculated every
	 major function call NOT within the class function call.  It is calculating the Bd_i that
	 are very time consuming as it involves looping over the entire spin list
	 'number of spins' times.

*/

/* This class also has the option of a 'scalefunction'
 * --the scale function 'scales' the caclulated Bd by some ammount as determined by the user
 * this is quite usefull to take into account 'other averaging effects' that are hard
 * to calcuate, or are too time consuming.
 *
 * for instance diffusion can roughly be acounted for by scaling the Bd by some
 * factor proportional to the distance away from the target point
 * the closer it is too the target point the LESS effect Bd has as diffusion will
 * certainly average the interaction out, but points far away will be almost
 * entirely inclded.
 *
 * few examples classes are found in 'scalefunc.h'::THEY MUST HAVE A FUNCTION CALLED 'FUNCTION'!!!
 * and take in the coord I the coord J, and the current magnitization coord
 *
 * function(coord<> &r, coord<> &rp, coord<> &Bd[i])
*/



template<class GridEng_t, class ScaleFunc_t=NoScaleFunc>
class DimLessDipole {

	private:

		GridEng_t *gridm;	//ptr to the grid
		Vector<coord<> > Bd_; //the vector of B_d...this only needs to calcualted
							//once upone every 'Bloch' function call

		Vector<coord<> > factors_;	//the acctual dipolar coupling value for each grid point
		Vector<coord<> > rhat_;	//these are used for calculating the rotated trunction axis (if nessesary)
		bool i_on;	// flag to turn the interaction off
		ScaleFunc_t *func;

		double fact_;
	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=1; //must perform a precalcuation before the main function can be used

		//do i need to recalc 'angles' and 'factors' for each precalculate?
		//i only need to do this if the Main Magnetic field field
		// changes in time (thus the 'Dynamic' flag)

		bool Dynamic;

		DimLessDipole():
			gridm(0), i_on(true), func(0), fact_(PI2), Dynamic(false)
		{}

		DimLessDipole(GridEng_t &gr, double fact=PI2):
			gridm(&gr), Bd_(gr.size()), factors_(gr.size()),
			i_on(true), func(0), fact_(fact), Dynamic(false)
		{}

		DimLessDipole(GridEng_t &gr, ScaleFunc_t &sc, double fact=PI2):
			gridm(&gr), Bd_(gr.size()), factors_(gr.size()),
			i_on(true), func(&sc), fact_(fact), Dynamic(false)
		{}

		DimLessDipole(const DimLessDipole &cp):
			gridm(cp.gridm), Bd_(cp.Bd_), factors_(cp.factors_),
			i_on(cp.i_on), func(cp.func), fact_(cp.fact_), Dynamic(cp.Dynamic)
		{}

		inline DimLessDipole &operator=(const DimLessDipole &cp);

		~DimLessDipole()
		{
			gridm=NULL;
			func=NULL;
		}

		void FactorsInit(double t); //used for the ;static fields' approx (or where Bz only amters)

		template<class PList>
		void FactorsInit(double t,PList *pars); //used for 'off axis' B fields

		inline GridEng_t *grid(){	return gridm;	}
		inline void setGrid(GridEng_t &in){	gridm=&in; Bd_.resize(in.size());	}

		inline ScaleFunc_t *ScaleFunc(){	return func;	}
		inline void SetScaleFunc(ScaleFunc_t &in){	func=&in; 	}

		inline double &factor(){	return fact_;	}
		inline void setFactor(double inf){	fact_=inf;	}


		void Off()	{	i_on=false;	}
		void On()	{	i_on=true;	}

		void off()	{	i_on=false;	}
		void on()	{	i_on=true;	}

		template<class GridEngine_t, int BPops>
		inline void preCalculate(double t,
								Vector<coord<> >  &M,
								Vector<coord<> >  &dMdt,
								ListBlochParams<GridEngine_t, BPops,double>  *pars,
								coord<> &totM);

//this one is for Bfields NOT on the 'Z' axis in the lab frame..
// here the dipolar interaction is trunctaed BUT along a different
//axis...leaving a painful task for dermining the axis of
// truncation for each spin (as each spins Bfield can be different)
		template<class GridEngine_t, int BPops>
		inline void preCalculate(double t,
								Vector<coord<> >  &M,
								Vector<coord<> >  &dMdt,
								ListBlochParams<GridEngine_t, BPops,coord<> >   *pars,
								coord<> &totM);


		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline coord<> function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

		template<class PList>
		inline void function(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix magneticField(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

		template<class PList>
		rmatrix evolutionMatrix(double t, Vector<coord<> > &M, PList *pars, coord<> &totM);

};

/*******************************************/
//IMPLIMETNATION
/*******************************************/
template<class GridEng_t, class ScaleFunc_t>
void
  DimLessDipole<GridEng_t,ScaleFunc_t>::FactorsInit(double t)
{

	if(!gridm)
	{
		std::cerr<<std::endl<<"Warning: DimLessDipole::FactorsInit(...)"<<std::endl;
		std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
		return;
	}

	typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
	typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2

	double dr3=min(GiterO.dr()*GiterO.dr()*GiterO.dr());
	double cosTH=0.0, diffLen=0.0, Z=0.0;
	coord<> ptJ, ptI;
	int i=0, j=0, ct=0;
	factors_.resize(gridm->size()*gridm->size());
	while(GiterO)
	{
		ptI=GiterO.Point(t); //ri
		GiterT.reset();
		j=0;
		dr3=min(GiterO.dr()*GiterO.dr()*GiterO.dr());
		while(GiterT)
		{
			factors_[ct]=0.0;
			if(j!=i)
			{
				/*calculating theta is a bit harder then it seems at first
				 the inital thing one would think would be to take a dot product
				 however, if our spin is on (0,0,0) then we have a problem (i.e. we always
				 get 0 for cos(theta)), so we need to break the problem into
				 a 2-D problem.  We are only interested in the z axis angle

						 o
					(zj) | \
						 |  \
						 Z   R
						 |    \
					(zi) --r--o\

				 here r=sqrt((xi-xj)^2+(yi-yj)^2)
				 we have taken our 2D x-y plane and flattened it to a new basis
				 where the combination basis is the R vector between spin one and spin two
				 The 'Z' hieght is then (zj-zi)
				 big 'R' is the distance btween the two pts
				 THe theta we disire STARTS at the z axis and continues to 'R' in  a clockwise sence
				 SOOO cos(theta)=Z/R
				*/
				ptJ=GiterT.Point(t); //rj
				Z=ptI.z()-ptJ.z();
				diffLen=norm(ptI-ptJ);	//|ri-rj|
				cosTH=Z/diffLen;
				if(diffLen)
				{
					//all three directions have the same factor...
					factors_[ct]=dr3*fact_*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3
					if(func) func->function(ptI, ptJ, diffLen, factors_[ct]);
				}

			}
			++ct;
			++j;
			++GiterT;
		}
		++i;
		++GiterO;
	}
}


template<class GridEng_t, class ScaleFunc_t>
template<class PList>
void
  DimLessDipole<GridEng_t,ScaleFunc_t>::
FactorsInit(double t,PList *pars) //used for 'off axis' B fields
{
	if(!gridm)
	{
		std::cerr<<std::endl<<"Warning: DimLessDipole::FactorsInit(...)"<<std::endl;
		std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
		return;
	}

	typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
	typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2

	double dr3;//=min(GiterO.dr()*GiterO.dr()*GiterO.dr());
	double diffLen;//cosTH=0.0, diffLen=0.0, Z=0.0, ZHATtheta, Phi;
	//coord<> ptJ, ptI;//, diffC;
	int i=0, j=0, ct=0;
	factors_.resize(gridm->size()*gridm->size());
	rhat_.resize(gridm->size()*gridm->size());
	while(GiterO)
	{
		//ptI=GiterO.Point(t); //ri
		GiterT.reset();
		j=0;
		dr3=cube(min(GiterO.dr()));
		while(GiterT)
		{
			factors_[ct]=0.0;
			if(j!=i )
			{

				diffLen=cube(norm(GiterO.Point(t)-GiterT.Point(t)));	//|ri-rj|

				if(diffLen)
				{
					factors_[ct]=dr3*fact_/(diffLen);
					//need the main 'rhat' vector...i.e. our truncation axis...
					rhat_[ct]=pars->Bo(GiterT.curpos())/norm(pars->Bo(GiterT.curpos()));
					if(func) func->function(GiterO.Point(t), GiterT.Point(t), diffLen, factors_[ct]);
				}

			}
			++ct;
			++j;
			++GiterT;
		}
		++i;
		++GiterO;
	}
}



//assignments....
template<class GridEng_t, class ScaleFunc_t>
inline DimLessDipole<GridEng_t,ScaleFunc_t> &
  DimLessDipole<GridEng_t,ScaleFunc_t>::operator=(const DimLessDipole<GridEng_t,ScaleFunc_t> &cp)
{
	if(this==&cp) return *this;
	gridm=(cp.gridm);
	Bd_=(cp.Bd_);
	factors_=cp.factors_;
	i_on=(cp.i_on);
	func=(cp.func);
	Dynamic=cp.Dynamic;
	return *this;
}

template<class GridEng_t, class ScaleFunc_t>
template<class GridEngine_t, int BPops>
inline void DimLessDipole<GridEng_t,ScaleFunc_t>::preCalculate(
	double t, Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	ListBlochParams<GridEngine_t, BPops,double>  *pars,
	coord<> &totM) 
{
	static double oldt=-1e30;
	if(i_on && fact_)
	{
		if(oldt==-1e30 || Dynamic) FactorsInit(t);
		if(oldt==t){ return;}
		else{	oldt=t;	}
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: DimLessDipole::preCalculate(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid..." )
		}
		int ct=0, i, j;
		for(i=0;i<gridm->size();++i)
		{
			Bd_[i]=0.0;
			for(j=0;j<gridm->size();++j)
			{
				if(i!=j){
					Bd_[i].x()-=factors_[ct].x()*M[j].x(); //-factor*Mx*(dr^3)+=Bd_x
					Bd_[i].y()-=factors_[ct].y()*M[j].y(); //-factor*My*(dr^3)+=Bd_y
					Bd_[i].z()+=2.0*factors_[ct].z()*M[j].z(); //2*factor*Mz*(dr^3)+=Bd_z
				}
				++ct;
			}
		}
	}
}

/***** OFFSET=Coord<> *****/
//this one is for Bfields NOT on the 'Z' axis in the lab frame..
// here the dipolar interaction is trunctaed BUT along a different
//axis...leaving a painful task for dermining the axis of
// truncation for each spin (as each spins Bfield can be different)

template<class GridEng_t, class ScaleFunc_t>
template<class GridEngine_t, int BPops>
inline void DimLessDipole<GridEng_t,ScaleFunc_t>::preCalculate(
	double t, Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	ListBlochParams<GridEngine_t, BPops,coord<> >  *pars,
	coord<> &totM) 
{
	static double oldt=-1e30;
	if(i_on && fact_)
	{
		if(oldt==-1e30 || Dynamic){ FactorsInit(t, pars);	}
		if(oldt==t){ return;}
		else{	oldt=t;	}

	//	coord<> tmpB;
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid..." )
		}
		int ct=0, i, j;
		for(i=0;i<gridm->size();++i)
		{
			Bd_[i]=0.0;
			for(j=0;j<gridm->size();++j)
			{
				if(i!=j){

					Bd_[i]+=factors_[ct]*(3.0*dot(M[j], rhat_[ct])*rhat_[ct]-M[j]);
					//cout<<factors_[ct]<<"||"<<rhat_[ct]<<" || "<<norm(rhat_[ct])<<" ||" <<Bd_[i]<<endl;
				}
				++ct;
			}
		}
	}
}


template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline rmatrix DimLessDipole<GridEng_t,ScaleFunc_t>::jacobian(
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	static rmatrix out(3,3,0);
	//static double fact=0.0;
	if(!i_on || !fact_) return rmatrix(3,3,0);
	out(0,0)=0.0;		out(0,1)=-Bd_[pars->Position()].z();		out(0,2)=Bd_[pars->Position()].y();
	out(1,0)=Bd_[pars->Position()].z();		out(1,1)=0.0;		out(1,2)=-Bd_[pars->Position()].x();
	out(2,0)=-Bd_[pars->Position()].y();		out(2,1)=Bd_[pars->Position()].x();		out(2,2)=0.0;
	return out;
}

template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline void DimLessDipole<GridEng_t,ScaleFunc_t>::jacobian(
	rmatrix &out,
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	//static double fact=0.0;
	if(!i_on || !fact_) return;
							out(0,1)-=Bd_[pars->Position()].z();		out(0,2)+=Bd_[pars->Position()].y();
	out(1,0)+=Bd_[pars->Position()].z();								out(1,2)-=Bd_[pars->Position()].x();
	out(2,0)-=Bd_[pars->Position()].y();		out(2,1)+=Bd_[pars->Position()].x();
}



//here we calculate M x B_s (M 'cross' B_s)
template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline coord<> DimLessDipole<GridEng_t,ScaleFunc_t>::function(
	double t,
	coord<>  &M,
	coord<>  &dMdt,
	ParamIter *pars,
	coord<> &totM)
{
	//  dM/dt=M x (gamma Bd)
	if(i_on && fact_)
	{
		//cout<<cross(M,Bd_[pars->Position()])<<endl;
		return cross(M,Bd_[pars->Position()]);


	}
	return ZeroType<coord<> >::zero();
}


/*template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline void
  DimLessDipole<GridEng_t,ScaleFunc_t>::function(
	double t,
	coord<>  &M,
	coord<>  &dMdt,
	ParamIter *pars,
	coord<> &totM)
{
	//  dM/dt=M x (gamma Bd)
	if(i_on)
	{
		dMdt.x()-=(Bd_[pars->Position()].z()*M.y()-Bd_[pars->Position()].y()*M.z());
		dMdt.y()-=(Bd_[pars->Position()].x()*M.z()-Bd_[pars->Position()].z()*M.x());
		dMdt.z()-=(Bd_[pars->Position()].y()*M.x()-Bd_[pars->Position()].x()*M.y());
	}
}
*/

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
inline void DimLessDipole<GridEng_t,ScaleFunc_t>::function(
	double t,
	Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM)
{
	if(i_on && fact_)
	{
		for(int i=0;i<M.size();++i)
		{
			dMdt[i].x()-=(Bd_[i].z()*M[i].y()-Bd_[i].y()*M[i].z());
			dMdt[i].y()-=(Bd_[i].x()*M[i].z()-Bd_[i].z()*M[i].x());
			dMdt[i].z()-=(Bd_[i].y()*M[i].x()-Bd_[i].x()*M[i].y());
		}
	}
}

//to place the Magnetic fields in the correct places in the matrix
// we must perform the 'precalc' calulation here as well....
// NOTE:: the SCALE FUNCTION IS NOT APPLIED HERE!!!!

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
rmatrix  DimLessDipole<GridEng_t,ScaleFunc_t>::magneticField(
	double t,
	Vector<coord<> > &M,
	PList *pars,
	coord<> &totM) 
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(i_on && fact_)
	{
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: DimLessDipole::magneticField(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return tmp;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid...")
		}
		typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
		typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2
	//	double fact;
	//	if(pars->SolveUnits()!=Normal){	fact=1.0;	}
		double cosTH=0.0, diffLen=0.0, factor=0.0, Z=0.0, dr3;
		coord<> ptJ, ptI;
		int i=0, j=0,cti=0, ctj=0;
		while(GiterO)
		{
			ptI=GiterO.Point(t); //ri
			GiterT.reset();
			j=0; ctj=0;
			dr3=min(GiterO.dr()*GiterO.dr()*GiterO.dr());
			while(GiterT)
			{
				if(j!=i)
				{
					ptJ=GiterT.Point(t); //rj
					Z=ptI.z()-ptJ.z();
					diffLen=norm(ptI-ptJ);	//|ri-rj|
					cosTH=Z/diffLen;
					if(diffLen)
					{
						factor=fact_*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3
						tmp(i,j)=-factor*dr3*M[ctj].x(); //-factor*Mx*(dr^3)+=Bd_x
						tmp(i+1, j+1)=-factor*dr3*M[ctj].y(); //-factor*My*(dr^3)+=Bd_y
						tmp(i+2, j+2)=2.0*dr3*factor*M[ctj].z(); //2*factor*Mz*(dr^3)+=Bd_z
	//					cout<<i<<" "<<j<<" ["<<ptI<<"] ["<<ptJ<<"] ("<<acos(cosTH)*180./Pi<<") ("<<M[ctj].z()<<") ["<<tmp(i,j)<<" "<<tmp(i+1,j+1)<<" "<<tmp(i+2,j+2)<<" ]"<<endl;
					}
				}
				j+=3;
				++ctj;
				++GiterT;
			}
			i+=3;
			++cti;
			++GiterO;
		}
	}
	return tmp;
}

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
rmatrix  DimLessDipole<GridEng_t,ScaleFunc_t>::evolutionMatrix(
	double t,
	Vector<coord<> > &M,
	PList *pars,
	coord<> &totM) 
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(i_on && fact_)
	{
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: DimLessDipole::evolutionMatrix(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return tmp;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid...")
		}
		typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
		typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2
	//	static double fact=permVac/(4.0*Pi);
	//	if(pars->SolveUnits()!=Normal){	fact=1.0;	}
		double cosTH=0.0, diffLen=0.0, factor=0.0, Z=0.0, dr3;
		coord<> ptJ, ptI;
		int i=0, j=0,cti=0, ctj=0;
		while(GiterO)
		{
			ptI=GiterO.Point(t); //ri
			GiterT.reset();
			j=0,ctj=0;
			dr3=min(GiterO.dr()*GiterO.dr()*GiterO.dr());
			while(GiterT)
			{
				if(j!=i)
				{
					ptJ=GiterT.Point(t); //rj
					Z=ptI.z()-ptJ.z();
					diffLen=norm(ptI-ptJ);	//|ri-rj|
					cosTH=Z/diffLen;
					if(diffLen)
					{
						factor=fact_*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3

						double bz=2.0*dr3*factor*M[ctj].z();
						double by=factor*dr3*M[ctj].y();
						double bx=-factor*dr3*M[ctj].x();
						tmp(i,j)= -(bz*M[cti].y()-by*M[cti].z());
						tmp(i+1, j+1)=-(bx*M[cti].z()-bz*M[cti].x());
						tmp(i+2, j+2)=-(by*M[cti].x()-bx*M[cti].y());
					}
				}
				j+=3;
				++ctj;
				++GiterT;
			}
			i+=3;
			++cti;
			++GiterO;
		}
	}
	return tmp;
}


template<class GridEng_t, class ScaleFunc_t>
std::ostream &operator<<(std::ostream &oo, DimLessDipole<GridEng_t, ScaleFunc_t> &out)
{
	oo<<"Dipole-Dipole Interactions"<<std::endl;
	oo<<*(out.Grid());
	oo<<*(out.ScaleFunc());
	return oo;
}

template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, DimLessDipole<GridEng_t, NoScaleFunc> &out)
{
	oo<<"Dipole-Dipole Interactions"<<std::endl;
	oo<<*(out.Grid());
	oo<<"No Scaleing Function";
	return oo;
}


template<class GridEng_t, class ScaleFunc_t>
std::ostream &operator<<(std::ostream &oo, DimLessDipole<GridEng_t, ScaleFunc_t> &out);

template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, DimLessDipole<GridEng_t, NoScaleFunc> &out);


END_BL_NAMESPACE


#endif


