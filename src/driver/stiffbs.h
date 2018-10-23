/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-13-01
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
	NOT WORKING!!!!!

	stiffbs.h-->Bulirsch-Stoer/Richard Extrapolation method diff eq solver
	for STIFF equations...

	this requires a Jacobian function...

*/

#ifndef _stiffbs_h_
#define _stiffbs_h_ 1


#include "container/containers.h"
#include "utils/utils.h"
#include "container/rankType.h"
#include "driver/basicODE.h"

BEGIN_BL_NAMESPACE


/*the Bulirsch-Stoer/Richard Extrapolation method
 *
 *this method should NOT be used for STIFF equations functions
 * STIFF-->diff eqs with 2 wildly different solutions one large one small
 *         They tend to make the integrator scream in pain unless care is taken
 *			to follow the path..i.e. we need a Jacobian to alert us to these
 *		    wild differences
 *
 *this one is based on the idea that the final answer is an analytic
 *solution of the step size...sooo what we do is extropoalate a bunch of
 *step sizes and FIT those points to some polynomial, then evaluate it at
 *a step size of 0
 *
 * The original BS method CANNOT handle stiff equations...this one is meant to
 * fix that problem....
 * Most of it is exactly the same except the need for the Jacobian and the change of
 * the function 'mmid' to a different semi-implicet version..
 *
 * The polynomial fit array also changes ('nseq')
 */

/*******REQUIREMENTS FOR THIS TO WORK!
*
* class func--> this is a class With the very Special Public function called "function"
*  Without that function there can be no usage of this class
*
*  "function"--> must take 3 arguments (Double Time, Container  &y, Container &dydt)...
*   it calculated the acctuall 'diffeq' function
*
*		void function(double time, Container &y, Container &dydt)
*
*  "jacobian"--> must take 3 arguments (Double Time, Container  &y, Matrix &dfdy)...
*   it calculated the acctuall 'diffeq' function
*
*		void jacobian(double time, Container &y, Container &dfdy)
*
* class T--> T can be almost anything so long these functions are defined for that class
*  max, min, astiffbs, +, -, *, /, *= +=, -=
*/
template<class func, class T=double, class Container=Vector<T> >
class stiffbs :
	public basicODE<T,Container>
{
	private:
		static const int kmaxx;
		static const int imaxx;
		static const double SC_max;


		//used in the function mmid (and is called ALOT, so
		//we avoid reinitiallization of these vectors over and over again)
		static Container yn;
		static Container ym;

		func *my_sys;
		//some variables that act as temporary places during the
		//integration from

		//some stuff for the polynomial extrapolation
		Container x;
		_matrix<T, FullMatrix> d;
	public:
		typedef _matrix<typename ObjectTrait<T>::numtype, FullMatrix> mattype;

	private:
		//the jacobian stuff
		mattype jacobi;
		mattype JacTemp; //a temporary

		//the dervative of the function with repect to time
		//calculated vai a simple finite difference
		Container dfdt;

		//items for the LU decomp
		Vector<int> rowPerm; //row permutations for the LUdecomp

		void set_yscal();

		//items for determining the order
		static Vector<double> a; //used to detemin the max order
		static rmatrix alf; //the work matrix
		static Vector<double> err; //local error
		int kopt;
		int kmax;
		int order;
		//we need to calc the max order only once
		bool found_order;

	public:
		typedef T numtype;
		typedef Container container;



	//defalt constructor, used mostly for testing and other classes if nessary
		stiffbs() :basicODE<T, Container>() {};

		stiffbs(func &mighty) :
			basicODE<T, Container>()
		{
			my_sys=&mighty;
		}

		stiffbs(const Container &iny, func &mighty) :
			basicODE<T, Container>()
		{
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		}

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value

		stiffbs(const Container &iny, double int1, double int2,func &mighty) :
			basicODE<T, Container> (int1, int2)
		{
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		};

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instep->initial step size guess

		stiffbs(const Container &iny, double int1, double int2, double instep=0.0, double instepmin=0.0 ) :
			basicODE<T, Container> (int1, int2, instep, instepmin)
		{
			my_sys=NULL;
			y=iny;
			num_func=y.size();
			reset();
		};

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instep->initial step size guess

		stiffbs(const Container &iny, double int1, double int2, double instep, func &mighty) :
			basicODE<T, Container> (int1, int2, instep)
		{
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		};

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instepmin->minimum step size to try
		//instep->initial step size guess

		stiffbs(Container &iny, double int1, double int2, double instep, double instepmin, func &mighty) :
			basicODE<T, Container>(int1, int2, instep, instepmin)
		{
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		};

		~stiffbs(){	my_sys=NULL;	}


	/*** Initial Condition Setters ***/
		void set_y(const Container &in)
		{ 	basicODE<T, Container>::set_y(in);	this->reset(); }

		void setY(const Container &in) { this->set_y(in);}
		void setInitialCondition(const Container &in) { this->set_y(in);}

		void reset();

		//finds the best order to integrate to
		void findorder();

		//this sucker is the meat of the matter
		//it finds the midpoint
		int mmid(Container &yO,double &step, double &nstep, Container &yseq);
		//here is the function that does our polynomial extrapolation
		void polyex(int &xiest, T &xest,Container &estim,Container &extrap);
		//here is a B-S step with internal error and step size
		//fixin'
		void stiffbsstep(Container & yO, Container &yscalO);
		//finally the master checker and solver...
		//the second returns the output array
		void odeint();

};

END_BL_NAMESPACE


#include "driver/stiffbs_meth.h"

#endif


