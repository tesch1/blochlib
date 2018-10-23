

/* bloch_int_diffu.h ********/


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
 	bloch_int_diffu.h-->Diffusion interactions for the Bloch equations
 */


#ifndef _bloch_interatction_diffusion_h_
#define _bloch_interatction_diffusion_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
//#include "stencils/stencilprep_func.h"

BEGIN_BL_NAMESPACE

/************ the long awatied diffusion...
this class behaves alot lihe the 'DemagField' class
there is the need for a 'precalulation' (we must calculate the Laplacian)

the interacting is quite simple from a functional point of view...
but quite difficult to calculate numerically within this frame work
of numerics....the optimizations up until now have been for
items where the 'neighboring' points did not really need to be known
for this interaction...we need to know neighbors and next nearest neighbors
becuase of this the 'stencil' operations were developed (and stand alone as there own)

BUT becuase this library uses unordered-nonuniform shapes so the neighbors are not
nessesarily known...there is a class that calculates all the neighbors
for any given 'shape' "StencilPrep'

this class embeds the Stencil Prep inside...class class is VERY expensive...but
luckally we only need to calculate it ONCE for a given shape

*/



template<class GridEng_t>
class Diffusion {

	private:

		GridEng_t *gridm;	//ptr to the grid
		Vector<coord<> > Diff_; //the vector of B_d...this only needs to calcualted
							//once upone every 'Bloch' function call
		coord<> D_; //the diffusion constants....
		bool i_on;	// flag to turn the interaction off

	public:

		static const int TotalMag=0;
		static const int TotalVector=0;
		static const int PreCalc=1; //must perform a precalcuation before the main function can be used

		Diffusion():
			gridm(0), i_on(true), D_(0)
		{}

		Diffusion(GridEng_t &gr):
			gridm(&gr), i_on(true), D_(0)
		{
			Diff_.resize(gr.size());
		}

		Diffusion(GridEng_t &gr, double D):
			gridm(&gr), i_on(true), D_(D)
		{
			Diff_.resize(gr.size());
		}

		Diffusion(GridEng_t &gr, double Di, double Dj, double Dk):
			gridm(&gr), i_on(true), D_(Di,Dj,Dk)
		{
			Diff_.resize(gr.size());
		}

		Diffusion(const Diffusion &cp):
			gridm(cp.gridm), Diff_(cp.Bd_), i_on(cp.i_on),
			D_(cp.D_)
		{}

		inline Diffusion &operator=(const Diffusion &cp);

		~Diffusion()
		{
			gridm=NULL;
		}

		inline GridEng_t *Grid(){	return gridm;	}
		inline void SetGrid(GridEng_t &in){	gridm=&in; Diff_.resize(in.size());	}

		double &Dx(){	return D_.x();	}
		double &Dy(){	return D_.y();	}
		double &Dz(){	return D_.z();	}

		coord<> &D(){	return D_;	}

		void Dx(double dx){	D_.x(dx);	}
		void Dy(double dy){	D_.y(dy);	}
		void Dz(double dz){	D_.z(dz);	}

		void D(const coord<> &d){	 D_=d;	}
		void D(double d){	 D_=d;	}


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
};


/******************** DIFFUSION *********************/


//assignments....
template<class GridEng_t>
inline Diffusion<GridEng_t> &
  Diffusion<GridEng_t>::operator=(const Diffusion<GridEng_t> &cp)
{
	if(this==&cp) return *this;
	gridm=(cp.gridm);
	Diff_=(cp.Diff_);
	D_=cp.D_;
	i_on=(cp.i_on);
	return *this;
}

END_BL_NAMESPACE

#endif


