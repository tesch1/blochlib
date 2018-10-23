/* gradimps.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-05-01
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
 	gradimps.h-->a set of specific classes designed to use the
 	classes in 'gradfunc.h' to perform and control implimentation
 	of temperature gradients and B_o inhomogenaity gradients

 	these classes determin weather or not the actuall gradient functions
 	should be implimented at each function call or simply at 'start up'

 	the 'start up' classes are used BEFORE the integration procedure begins

 	and thus are only appiled ONCE...These classes are used by 'ListBlochParams'
 	to set the initial condition

 	the 'Continuous' implimentations are used for variational changes as time
 	marchs on...

 	the Continuous implimentation are used in the 'Diffusion' interaction..
 	(in fact there is NO reason to use the contious version when there is no diffusion
 	as everything maintains is local assigned value for all time)

 */




#ifndef _grad_imps_h_
#define _grad_imps_h_ 1


#include "bloch/gradfunc.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


enum Application_{Startup, Continuous};



template<class GradFunc_t>
class TemperatureGrad :
	public GradFunc_t
{
	private:
		Application_ apply_;

	public:

	//various ctors that link to the GRadFunc
		TemperatureGrad():
			GradFunc_t(),
			apply_(Startup)
		{}

		TemperatureGrad(double mx, double bx):
			GradFunc_t(mx,bx),
			apply_(Startup)
		{}

		TemperatureGrad(double mx, double bx, double off):
			GradFunc_t(mx,bx,off),
			apply_(Startup)
		{}

		TemperatureGrad(double mx, double bx, double off, double sc):
			GradFunc_t(mx,bx,off,sc),
			apply_(Startup)
		{}

		TemperatureGrad(TemperatureGrad &cp):
			GradFunc_t(cp),
			apply_(cp.apply_)
		{}

	//set and get the aplication type
		inline Application_ Application(){		return apply_;	}
		void SetApplication(Application_ in){	apply_=in;	}


	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it
		template<class Params_t>
		void apply(Params_t &in)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				GradFunc_t::operator()(myit.temperature(), myit.Point());
				if(myit.temperature()<0.0){	myit.temperature(1.e-6); }
				myit.calcMo();
				++myit;
			}
			in.calcTotalMo();
		}

	//same func...but in operator style
		template<class Params_t>
		void operator()(Params_t &in)
		{
			apply(in);
		}

	//The 'Shape Expression' functions...will apply the gradient
	// to only those things given by 'ShapeFunc' from the SphapeExpr object

	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it..
		template<class Params_t, class ShapeExpr>
		void apply(Params_t &in, ShapeExpr expr)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				if(expr.ShapeFunc(myit.Point()))
				{
					GradFunc_t::operator()(myit.temperature(), myit.Point());
					if(myit.temperature()<0.0){	myit.temperature(1.e-6); }
				}
				++myit;

			}
		}

	//same func...but in operator style
		template<class Params_t, class ShapeExpr>
		void operator()(Params_t &in, ShapeExpr expr)
		{
			apply(in,expr);
		}


};

template<class GradFunc_t>
std::ostream &operator<<(std::ostream &oo, TemperatureGrad<GradFunc_t> &out)
{
	oo<<"temperature Gradient: ";
	GradFunc_t tmp(out);
	oo<<tmp;
	return oo;
}


//A B_o inhomogenaity...THE VALUES SHOULD BE IN Hz (NOT Telsa)
// alters the OFFSET according to the grid position

template<class GradFunc_t>
class InhomogenaityGrad :
	public GradFunc_t
{
	private:
		Application_ apply_;

	public:

	//various ctors that link to the GRadFunc
		InhomogenaityGrad():
			GradFunc_t(),
			apply_(Startup)
		{}

		InhomogenaityGrad(double mx, double bx):
			GradFunc_t(mx,bx),
			apply_(Startup)
		{}

		InhomogenaityGrad(double mx, double bx, double off):
			GradFunc_t(mx,bx,off),
			apply_(Startup)
		{}

		InhomogenaityGrad(double mx, double bx, double off, double sc):
			GradFunc_t(mx,bx,off,sc),
			apply_(Startup)
		{}

		InhomogenaityGrad(InhomogenaityGrad &cp):
			GradFunc_t(cp),
			apply_(cp.apply_)
		{}

	//set and get the aplication type
		inline Application_ Application(){		return apply_;	}
		void SetApplication(Application_ in){	apply_=in;	}


	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it
		template<class Params_t>
		void apply(Params_t &in)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				GradFunc_t::operator()(myit.spin_offset(), myit.Point());
				++myit;
			}
		}

	//same func...but in operator style
		template<class Params_t>
		void operator()(Params_t &in)
		{
			apply(in);
		}

	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it..
		template<class Params_t, class ShapeExpr>
		void apply(Params_t &in, ShapeExpr expr)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				if(expr.ShapeFunc(myit.Point()))
				{
					GradFunc_t::operator()(myit.spin_offset(), myit.Point());
				}
				++myit;

			}
		}

	//same func...but in operator style
		template<class Params_t, class ShapeExpr>
		void operator()(Params_t &in, ShapeExpr expr)
		{
			apply(in,expr);
		}
};

template<class GradFunc_t>
std::ostream &operator<<(std::ostream &oo, InhomogenaityGrad<GradFunc_t> &out)
{
	oo<<"Inhomogenaity Gradient: ";
	GradFunc_t tmp(out);
	oo<<tmp;
	return oo;
}

//A Mole gradients...THE VALUES SHOULD BE IN moles
// alters the MOLES according to the grid position

template<class GradFunc_t>
class MoleGrad :
	public GradFunc_t
{
	private:
		Application_ apply_;

	public:

	//various ctors that link to the GRadFunc
		MoleGrad():
			GradFunc_t(),
			apply_(Startup)
		{}

		MoleGrad(double mx, double bx):
			GradFunc_t(mx,bx),
			apply_(Startup)
		{}

		MoleGrad(double mx, double bx, double off):
			GradFunc_t(mx,bx,off),
			apply_(Startup)
		{}

		MoleGrad(double mx, double bx, double off, double sc):
			GradFunc_t(mx,bx,off,sc),
			apply_(Startup)
		{}

		MoleGrad(MoleGrad &cp):
			GradFunc_t(cp),
			apply_(cp.apply_)
		{}

	//set and get the aplication type
		inline Application_ Application(){		return apply_;	}
		void SetApplication(Application_ in){	apply_=in;	}


	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it
		template<class Params_t>
		void apply(Params_t &in)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				GradFunc_t::operator()(myit.moles(), myit.Point());
				if(myit.moles()<0.0){	myit.moles(0.0);	}
				++myit;
			}
		}

	//same func...but in operator style
		template<class Params_t>
		void operator()(Params_t &in)
		{
			apply(in);
		}

	//the main function...takes in a 'Parameter Set' (i.e. a ListBlochParams<engine_t>
	// and applies the enitre gradient type to it..
		template<class Params_t, class ShapeExpr>
		void apply(Params_t &in, ShapeExpr expr)
		{
			typename Params_t::iterator myit(in);
			while(myit)
			{
				if(expr.ShapeFunc(myit.Point()))
				{
					GradFunc_t::operator()(myit.moles(), myit.Point());
					if(myit.moles()<0.0){	myit.moles(0.0);	}
				}
				++myit;
				in.calcTotalMo();
			}
		}

	//same func...but in operator style
		template<class Params_t, class ShapeExpr>
		void operator()(Params_t &in, ShapeExpr expr)
		{
			apply(in,expr);
		}
};

template<class GradFunc_t>
std::ostream &operator<<(std::ostream &oo, MoleGrad<GradFunc_t> &out)
{
	oo<<"Mole Gradient: ";
	GradFunc_t tmp(out);
	oo<<tmp;
	return oo;
}


END_BL_NAMESPACE


#endif





