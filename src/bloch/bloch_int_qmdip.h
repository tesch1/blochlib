


/* bloch_int_qmdip.h ********/


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
 	bloch_int_dipdip.h-->'Quantum Mechanincal' Dipole-Dipole Interactions
 	for the Bloch equations

 	NOTE:: THIS INTERACTION IS Physically strange...for now it was just to see some numerical
 	scaling problems....not really a valid interaction at all
 */


#ifndef _bloch_interatction_QM_dipole_dipole_h_
#define _bloch_interatction_QM_dipole_dipole_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "bloch/scalefunc.h"

BEGIN_BL_NAMESPACE

/* Quantum Mechanical dipole-dipole  coupling

	the only difference between this and the 'QMDipoleDipole' class
	is the way the Magnetic field is caluclated...it is caluclates AS IF
	the two spins where 'atomic'/quantum mechanical according to the formulae below

	This functions requires ALL the locations of the spins and thus the GRID!!

		dm/dt=gamma*hbar (M x Bd)
		where Bd_i = gamma_i*hbar Sum_j[ (1-3cos(th_ri_rj))/2*|ri-rj|^3
		              * [3Mz(rj)*zhat-Mx xhat - My yhat - Mz zhat]]


	Bd_i is the effective magnetic field at point 'ri' generated by all the little
	 dipoles around it.  Each Bd_i thus must be calculated every
	 major function call NOT within the class function call.  It is calculating the Bd_i that
	 are very time consuming as it involves looping over the entire spin list
	 'number of spins' times.

	This is VERY close to the 'DemagField' class except that here we are assuming individual spins

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
class QMDipoleDipole {

	private:

		GridEng_t *gridm;	//ptr to the grid
		Vector<coord<> > Bd_; //the vector of B_d...this only needs to calcualted
							//once upone every 'Bloch' function call
		bool i_on;	// flag to turn the interaction off
		ScaleFunc_t *func;

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=1; //must perform a precalcuation before the main function can be used

		QMDipoleDipole():
			gridm(0), i_on(true), func(0)
		{}

		QMDipoleDipole(GridEng_t &gr):
			gridm(&gr), i_on(true), func(0)
		{
			Bd_.resize(gr.size());
		}

		QMDipoleDipole(GridEng_t &gr, ScaleFunc_t &sc):
			gridm(&gr), i_on(true), func(&sc)
		{
			Bd_.resize(gr.size());
		}

		QMDipoleDipole(const QMDipoleDipole &cp):
			gridm(cp.gridm), Bd_(cp.Bd_), i_on(cp.i_on), func(cp.func)
		{}

		inline QMDipoleDipole &operator=(const QMDipoleDipole &cp);

		~QMDipoleDipole()
		{
			gridm=NULL;
			func=NULL;
		}

		inline GridEng_t *Grid(){	return gridm;	}
		inline void SetGrid(GridEng_t &in){	gridm=&in; Bd_.resize(in.size());	}

		inline ScaleFunc_t *ScaleFunc(){	return func;	}
		inline void SetScaleFunc(ScaleFunc_t &in){	func=&in; 	}



		void Off()	{	i_on=false;	}
		void On()	{	i_on=true;	}

		template<class PList>
		inline void preCalculate(double t, Vector<coord<> >  &M, Vector<coord<> >  &dMdt, PList *pars, coord<> &totM);

		template<class ParamIter>
		inline rmatrix jacobian(double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		template<class ParamIter>
		inline void jacobian(rmatrix &out, double t,coord<>  &M, ParamIter *pars,  coord<> &totM);

		//here we calculate M x B_s (M 'cross' B_s)
		template<class ParamIter>
		inline void function(double t, coord<>  &M, coord<>  &dMdt, ParamIter *pars, coord<> &totM);

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


//assignments....
template<class GridEng_t, class ScaleFunc_t>
inline QMDipoleDipole<GridEng_t,ScaleFunc_t> &
  QMDipoleDipole<GridEng_t,ScaleFunc_t>::operator=(const QMDipoleDipole<GridEng_t,ScaleFunc_t> &cp)
{
	if(this==&cp) return *this;
	gridm=(cp.gridm);
	Bd_=(cp.Bd_);
	i_on=(cp.i_on);
	func=(cp.func);
	return *this;
}

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
inline void QMDipoleDipole<GridEng_t,ScaleFunc_t>::preCalculate(
	double t, Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM) 
{
	if(i_on)
	{
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: QMDipoleDipole::preCalculate(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid...")
		}
		typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
		typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2
		double fact;
		double cosTH=0.0, diffLen=0.0, factor=0.0, Z=0.0;
		coord<> ptJ, ptI;
		int i=0, j=0;
		while(GiterO)
		{
			ptI=GiterO.Point(); //ri
			GiterT.reset();
			Bd_[i]=0.0;
			j=0;
			typename PList::iterator ParamIT(*pars); //parameter iterator
			while(GiterT)
			{
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

					fact=ParamIT.gamma()*hbar;
					ptJ=GiterT.Point(); //rj
					Z=ptI.z()-ptJ.z();
					diffLen=norm(ptI-ptJ);	//|ri-rj|
					cosTH=Z/diffLen;
					if(diffLen)
					{
						//if(pars->SolveUnits()==Normal)
						//{
							factor=fact*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3
							Bd_[i].x()-=factor*M[j].x(); //-factor*Mx*(dr^3)+=Bd_x
							Bd_[i].y()-=factor*M[j].y(); //-factor*My*(dr^3)+=Bd_y
							Bd_[i].z()+=2.0*factor*M[j].z(); //2*factor*Mz*(dr^3)+=Bd_z
						/*}
						else //Dimmensionless UNITS....r-->r/dr
						{
							factor=fact*(1.0-3.0*cosTH*cosTH)/(2.0); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3
							double difx=diffLen/GiterT.dx(), dify=diffLen/GiterT.dy(), difz=diffLen/GiterT.dz();
							Bd_[i].x()-=factor*M[j].x()/(difx*difx*difx); //-factor*Mx*(DeltaX)+=Bd_x
							Bd_[i].y()-=factor*M[j].y()/(dify*dify*dify); //-factor*My*(DeltaY)+=Bd_y
							Bd_[i].z()+=2.0*factor*M[j].z()/(difz*difz*difz); //2*factor*Mz*(DeltaZ)+=Bd_z
						}*/
						if(func) func->function(ptI, ptJ, diffLen, Bd_[i]);
					}

				}
				++j;
				++GiterT;
				++ParamIT;
			}
			++i;
			++GiterO;
		}
	}
}

template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline rmatrix QMDipoleDipole<GridEng_t,ScaleFunc_t>::jacobian(
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	static rmatrix out(3,3,0);
	static double fact=0.0;
	if(!i_on) return rmatrix(3,3,0);
	fact=pars->gamma()*hbar;
	out(0,0)=0.0;		out(0,1)=-fact*Bd_[pars->Position()].z();		out(0,2)=fact*Bd_[pars->Position()].y();
	out(1,0)=fact*Bd_[pars->Position()].z();		out(1,1)=0.0;		out(1,2)=-fact*Bd_[pars->Position()].x();
	out(2,0)=-fact*Bd_[pars->Position()].y();		out(2,1)=fact*Bd_[pars->Position()].x();		out(2,2)=0.0;
	return out;
}

template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline void QMDipoleDipole<GridEng_t,ScaleFunc_t>::jacobian(
	rmatrix &out,
	double t,coord<>  &M,
	ParamIter *pars,
	coord<> &totM)
{
	static double fact=0.0;
	if(!i_on) return;
	fact=pars->gamma()*hbar;
							out(0,1)-=fact*Bd_[pars->Position()].z();		out(0,2)+=fact*Bd_[pars->Position()].y();
	out(1,0)+=fact*Bd_[pars->Position()].z();								out(1,2)-=fact*Bd_[pars->Position()].x();
	out(2,0)-=fact*Bd_[pars->Position()].y();		out(2,1)+=fact*Bd_[pars->Position()].x();
}



//here we calculate M x B_s (M 'cross' B_s)
template<class GridEng_t, class ScaleFunc_t>
template<class ParamIter>
inline void QMDipoleDipole<GridEng_t,ScaleFunc_t>::function(
	double t,
	coord<>  &M,
	coord<>  &dMdt,
	ParamIter *pars,
	coord<> &totM)
{
	//  dM/dt=M x (gamma Bd)
	if(i_on)
	{
		dMdt.x()-=pars->gamma()*(Bd_[pars->Position()].z()*M.y()-Bd_[pars->Position()].y()*M.z());
		dMdt.y()-=pars->gamma()*(Bd_[pars->Position()].x()*M.z()-Bd_[pars->Position()].z()*M.x());
		dMdt.z()-=pars->gamma()*(Bd_[pars->Position()].y()*M.x()-Bd_[pars->Position()].x()*M.y());
	}
}

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
inline void QMDipoleDipole<GridEng_t,ScaleFunc_t>::function(
	double t,
	Vector<coord<> >  &M,
	Vector<coord<> >  &dMdt,
	PList *pars,
	coord<> &totM)
{
	if(i_on)
	{
		for(int i=0;i<M.size();++i)
		{
			dMdt[i].x()-=pars->gamma(i)*(Bd_[i].z()*M[i].y()-Bd_[i].y()*M[i].z());
			dMdt[i].y()-=pars->gamma(i)*(Bd_[i].x()*M[i].z()-Bd_[i].z()*M[i].x());
			dMdt[i].z()-=pars->gamma(i)*(Bd_[i].y()*M[i].x()-Bd_[i].x()*M[i].y());
		}
	}
}

//to place the Magnetic fields in the correct places in the matrix
// we must perform the 'precalc' calulation here as well....
// NOTE:: the SCALE FUNCTION IS NOT APPLIED HERE!!!!

template<class GridEng_t, class ScaleFunc_t>
template<class PList>
rmatrix  QMDipoleDipole<GridEng_t,ScaleFunc_t>::magneticField(
	double t,
	Vector<coord<> > &M,
	PList *pars,
	coord<> &totM) 
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(i_on)
	{
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: QMDipoleDipole::magneticField(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return tmp;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid...")
		}
		typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
		typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2
		double fact;
		if(pars->SolveUnits()!=Normal){	fact=1.0;	}
		double cosTH=0.0, diffLen=0.0, factor=0.0, Z=0.0;
		coord<> ptJ, ptI;
		int i=0, j=0,cti=0, ctj=0;
		while(GiterO)
		{
			ptI=GiterO.Point(); //ri
			GiterT.reset();
			j=0; ctj=0;
			typename PList::iterator ParamIT(*pars); //parameter iterator
			while(GiterT)
			{
				if(j!=i)
				{
					fact=ParamIT.gamma()*hbar;
					ptJ=GiterT.Point(); //rj
					Z=ptI.z()-ptJ.z();
					diffLen=norm(ptI-ptJ);	//|ri-rj|
					cosTH=Z/diffLen;
					if(diffLen)
					{
						factor=fact*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3
						tmp(i,j)=-factor*M[ctj].x(); //-factor*Mx*(dr^3)+=Bd_x
						tmp(i+1, j+1)=-factor*M[ctj].y(); //-factor*My*(dr^3)+=Bd_y
						tmp(i+2, j+2)=2.0*factor*M[ctj].z(); //2*factor*Mz*(dr^3)+=Bd_z
						cout<<i<<" "<<j<<" ["<<ptI<<"] ["<<ptJ<<"] ("<<acos(cosTH)*180./Pi<<") ("<<M[ctj].z()<<") ["<<tmp(i,j)<<" "<<tmp(i+1,j+1)<<" "<<tmp(i+2,j+2)<<" ]"<<endl;
					}
				}
				j+=3;
				++ctj;
				++GiterT;
				++ParamIT;
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
rmatrix  QMDipoleDipole<GridEng_t,ScaleFunc_t>::evolutionMatrix(
	double t,
	Vector<coord<> > &M,
	PList *pars,
	coord<> &totM) 
{
	rmatrix tmp(pars->size()*3, pars->size()*3,0);
	if(i_on)
	{
		if(!gridm)
		{
			std::cerr<<std::endl<<"Warning: QMDipoleDipole::evolutionMatrix(...)"<<std::endl;
			std::cerr<<" No grid defined, cannot calculate interaction..."<<std::endl;
			return tmp;
		}
		if(gridm->size() != pars->size())
		{
			BLEXCEPTION(" The parameter list is not the same size as the grid...")
		}
		typename GridEng_t::iterator GiterO((*gridm)); //grid iterator 1
		typename GridEng_t::iterator GiterT((*gridm)); //grid iterator 2
		double fact;
		if(pars->SolveUnits()!=Normal){	fact=1.0;	}
		double cosTH=0.0, diffLen=0.0, factor=0.0, Z=0.0;
		coord<> ptJ, ptI;
		int i=0, j=0,cti=0, ctj=0;
		while(GiterO)
		{
			ptI=GiterO.Point(); //ri
			GiterT.reset();
			j=0,ctj=0;
			typename PList::iterator ParamIT(*pars); //parameter iterator
			while(GiterT)
			{
				if(j!=i)
				{
					fact=ParamIT.gamma()*hbar;
					ptJ=GiterT.Point(); //rj
					Z=ptI.z()-ptJ.z();
					diffLen=norm(ptI-ptJ);	//|ri-rj|
					cosTH=Z/diffLen;
					if(diffLen)
					{
						factor=fact*(1.0-3.0*cosTH*cosTH)/(2.0*diffLen*diffLen*diffLen); //permVac/4Pi *(1-3cos^2(th))/2*|ri-rj|^3

						double bz=2.0*factor*M[ctj].z();
						double by=factor*M[ctj].y();
						double bx=-factor*M[ctj].x();
						if(pars->SolveUnits()!=Normal){
							tmp(i,j)= -(bz*M[cti].y()-by*M[cti].z());
							tmp(i+1, j+1)=-(bx*M[cti].z()-bz*M[cti].x());
							tmp(i+2, j+2)=-(by*M[cti].x()-bx*M[cti].y());
						}else{
							tmp(i,j)= -pars->gamma(cti)*(bz*M[cti].y()-by*M[cti].z());
							tmp(i+1, j+1)=-pars->gamma(cti)*(bx*M[cti].z()-bz*M[cti].x());
							tmp(i+2, j+2)=-pars->gamma(cti)*(by*M[cti].x()-bx*M[cti].y());
						}

					}
				}
				j+=3;
				++ctj;
				++GiterT;
				++ParamIT;
			}
			i+=3;
			++cti;
			++GiterO;
		}
	}
	return tmp;
}


template<class GridEng_t, class ScaleFunc_t>
std::ostream &operator<<(std::ostream &oo, QMDipoleDipole<GridEng_t, ScaleFunc_t> &out)
{
	oo<<"QMDipole-Dipole Interactions"<<std::endl;
	oo<<*(out.Grid());
	oo<<*(out.ScaleFunc());
	return oo;
}

template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, QMDipoleDipole<GridEng_t, NoScaleFunc> &out)
{
	oo<<"QMDipole-Dipole Interactions"<<std::endl;
	oo<<*(out.Grid());
	oo<<"No Scaleing Function";
	return oo;
}


template<class GridEng_t, class ScaleFunc_t>
std::ostream &operator<<(std::ostream &oo, QMDipoleDipole<GridEng_t, ScaleFunc_t> &out);

template<class GridEng_t>
std::ostream &operator<<(std::ostream &oo, QMDipoleDipole<GridEng_t, NoScaleFunc> &out);



END_BL_NAMESPACE


#endif
