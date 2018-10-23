
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
	bs.h-->Bulirsch-Stoer/Richard Extrapolation method diff eq solver

*/

#ifndef _bs_h_
#define _bs_h_ 1


#include "container/containers.h"
#include "utils/utils.h"
#include "driver/basicODE.h"


BEGIN_BL_NAMESPACE


/*the Bulirsch-Stoer/Richard Extrapolation method
 *
 *this method should NOT be used for non smooth functions
 *R-K does this much better
 *
 *this one is based on the idea that the final answer is an analytic
 *solution of the step size...sooo what we do is extropoalate a bunch of
 *step sizes and FIT those points to some polynomial, then evaluate it at
 *a step size of 0
 */

/*******REQUIREMENTS FOR THIS TO WORK!
*
* class func--> this is a class With the very Special Public function called "function"
*  Without that function there can be no usage of this class
*  "function"--> must take 3 arguments (Double Time, T  &y, T &dydt)...
*   it calculated the acctuall 'diffeq' function
*
* class T--> T can be almost anything so long these functions are defined for that class
*  max, min, abs, +, -, *, /, *= +=, -=
*/
template<class func, class T=double, class Container=Vector<T> >
class bs :
	public basicODE<T, Container> {
	private:
		//the various constants
		//kmaxx->maximum number of sub divisions per large stepsize H
		//imaxx->kmaxx+1
		static const int kmaxx;
		static const int imaxx;
		static const double SC_max;

		//used in the function mmid (and is called ALOT, so
		//we avoid reinitiallization of these vectors over and over again)
		static Container yn;
		static Container ym;

		func *my_sys;

		//some stuff for the polynomial extrapolation
        coord<T, 15>  x;
		_matrix<T, FullMatrix> d;

	//sets the initial yscal guess values
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
		bs():
			basicODE<T, Container> ()
		{};


		bs(func &mighty) :
			basicODE<T, Container>()
		{
			my_sys=&mighty;
		}

		bs(const Container &iny, func &mighty) :
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
		//instep->initial step size guess

		bs(const Container &iny, double int1, double int2, double instep, func &mighty) :
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

		bs(const Container &iny, double int1, double int2,func &mighty) :
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
		//instepmin->minimum step size to try
		//instep->initial step size guess

		bs(const Container &iny, double int1, double int2, double instep, double instepmin, func &mighty) :
			basicODE<T, Container> (int1, int2, instep, instepmin)
		{
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		};

		virtual ~bs(){	my_sys=NULL;	}

	/*** Initial Condition Setters ***/
		void set_y(const Container &in)
		{ 	basicODE<T, Container>::set_y(in);	this->reset(); }

		void setY(const Container &in) { this->set_y(in);}
		void setInitialCondition(const Container &in) { this->set_y(in); reset();}


		void reset();


	/*** The MASTER function ***/
	//finds the best order to integrate to based on eps
		void findorder();

	//it finds the midpoint
		void mmid(Container &yO,double &step, double &nstep, Container &yseq);
	//here is the function that does our polynomial extrapolation
		void polyex(int &xiest, T &xest,Container &estim,Container &extrap);
	//here is a B-S step with internal error and step size
	//fixin'
		void bsstep(Container & yO, Container &yscalO);


	//finally the master checker and solver...
	//the second returns the output array
		void odeint();

};

END_BL_NAMESPACE

#include "driver/bs_meth.h"

#endif

