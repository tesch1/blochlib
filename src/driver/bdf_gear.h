

BEGIN_BL_NAMESPACE

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
	bdf.h-->a backward differentiation formulas  solver

*/

#ifndef _bdf_gear_h_
#define _bdf_gear_h_ 1


#include "container/containers.h"
#include "utils/utils.h"
#include "driver/basicODE.h"


/* this uses a maximal 5th order solver for solving STIFF DiffEqs
 *
 	It relies on both the Jacobian and a Function from the class 'func'
 	It uses a Gear type backward differentiation formulas (hence the 'bdf')

 	it is designed to be the rough equivilent of the 'ode15s' in
 	MatLab

 */

/*******REQUIREMENTS FOR THIS TO WORK!
*
* class func--> this is a class With the very Special Public function called "function"
*  Without that function there can be no usage of this class
*  "function"--> must take 3 arguments
      void (double Time, Container  &y, Container &dydt)...
   it calculated the acctuall 'diffeq' function
*
*   "jacobian" --> must take 4 arguments
      void jacobian(double time, Container, &y, Container, &dydt, MatType & dfdy)

      where "MatType" is "_matrix<ObjectTrait<T>::numtype, FullMatrix>"

* class T--> T can be almost anything so long these functions are defined for that class
*  max, min, abs, +, -, *, /, *= +=, -=
*/
template<class func, class T=double, class Container=Vector<T> >
class bdf :
	public basicODE<T, Container>
{
	private:
		//the various constants
		// MaxOrder-->the maximal order of the intragation
		// MaxIter-->the maximal sub iteration level
		// G-->holds the coiefs for the integrator
		// inv-->the 'inverse' of G
		// alpha --> the coeifs for the BDF (set to 0 if NO BDF)
		// errConsts--> the error associated with each order params in G and alpha
		// difU --> the transition matrix (weights of each factor
		//          when we are starting with the differentiation formulas
		//          i.e. we need to start at 1, then move to MaxOrder)
		// kI & kJ--> index matrices for 'cumlative product differences'
		//            where each element element difference needs to be 1
		// order--> the order of the integratrion 1...5
		//
		static const int MaxOrder;
		static const int MaxIter;
		static Vector<double> G;
		static Vector<double> invG;
		static Vector<double> alpha;
		static Vector<double> errConsts;
		static rmatrix difU;
		static rmatrix kJ;
		static rmatrix kI;
		int order;
		int k; //our current Order (as we march through the inital points)
		int lastk; //the last Order (as we march through the inital points)
		Range K; //the 1..k vector
		double rate; //the rate of stepsize increase/decrease
		bool have_rate; //wether or not we have a 'rate' os stepping

	//initialized all the above vectors and matrices
		void initConsts(int Order);

		func *my_sys;

		//a holder for the backwards diffs
		Vector<Container> dif;
	public:
		typedef _matrix<ObjectTrait<T>::numtype, FullMatrix> mattype;

	private:

		//the Jacobian matrix
		mattype dfdy;
		mattype JacTemp; //a tmp matrix

		//LU decompo matricies
		Vector<int> rowPerm;
		bool got_jacobi;

	//have we 'initialized' yet?
		bool have_started;


	public:
		typedef T numtype;
		typedef Container container;


	//defalt constructor, used mostly for testing and other classes if nessary
		bdf():
			basicODE<T, Container> (),
			have_rate(false), have_started(false)
		{
			initConsts(MaxOrder);
		};

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instep->initial step size guess

		bdf(const Container &iny, double int1, double int2, double instep, func &mighty) :
			basicODE<T, Container> (int1, int2, instep),
			have_rate(false),have_started(false)
		{
			initConsts(MaxOrder);
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

		bdf(const Container &iny, double int1, double int2, double instep, double instepmin, func &mighty) :
			basicODE<T, Container> (int1, int2, instep, instepmin),
			have_started(false)
		{
			initConsts(MaxOrder);
			my_sys=&mighty;
			y=iny;
			num_func=y.size();
			reset();
		};

		~bdf(){	my_sys=NULL;	}

	/*** Initial Condition Setters ***/
		void set_y(const Container &in)
		{ 	basicODE<T, Container>::set_y(in);	this->reset(); }

		void setY(const Container &in)
		{ this->set_y(in);}
		void setInitialCondition(const Container &in)
		{ this->set_y(in);}


		void reset();


	/*** The MASTER function ***/

	//it finds the midpoint
	//this functions starts the integrators various parameters
		void start();

	//finally the master checker and solver...
	//the second returns the output array
		void odeint();
};

END_BL_NAMESPACE

#include "driver/bdf_gear_meth.h"


#endif


