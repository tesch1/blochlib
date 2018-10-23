


/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	ckrk.h-->a 5th order Runge-Kutta Diff eq solver

*/

#ifndef _ckrk_h_
#define _ckrk_h_ 1

#include "container/containers.h"
#include "utils/utils.h"
#include "driver/basicODE.h"


BEGIN_BL_NAMESPACE

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
template< class func, class T=double, class Container=Vector<T> >
class ckrk :
	public basicODE<T, Container>
{
	private:
		//the various constants
		//grow->the power to grow the step size by
		//shrink->the power to shrink our step size
		//safty->incase a step goes a bit silly, we will help it out
		//err_max->is our accuracy goal
		//max_step->maximum number of steps to take
		//max_step_within->within a given step, this tells us the maximum num times
		//		to keep moving the step size
		//tiny->incase we get stuck around '0' add tiny to keep things goin'
		//eps->our presision goal
		static double grow;
		static double shrink;
		static int max_step_within;



		//these are for the integrator part...they are the various coefs needed
		//for rk part
		static const float a2, a3 , a4,  a5,a6;
		static const float b21;
		static const float b31,b32,b41,b42,b43;
		static const float b51, b52,b53,b54;
		static const float b61,b62,b63;
		static const float b64,b65,c1;
		static const float c3,c4,c6;
		static const float dc1,dc3;
		static const float dc4,dc5,dc6;

		func *my_sys;

	//sets the initial yscal guess values
		void set_yscal();


		public:

			typedef T numtype;
			typedef Container container;

		//defalt constructor, used mostly for testing and other classes if nessary
			ckrk():
				basicODE<T, Container>()
			{}

			ckrk(func &mighty) :
				basicODE<T, Container>()
			{
				my_sys=&mighty;
			}

			ckrk(const Container &iny, func &mighty) :
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
			ckrk(const Container &iny, double int1, double int2,func &mighty) :
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
		//num->number of ODE/functions
			ckrk(Container &iny, double int1, double int2, double instep, double instepmin, func &mighty):
				basicODE<T, Container>(int1, int2, instep, instepmin)
			{
					my_sys=&mighty;
					y=iny;
					num_func=iny.size();
					reset();
			};

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instep->initial step size guess
		//num->number of ODE/functions
			ckrk(const Container &iny, double int1, double int2, double instep, func &mighty):
				basicODE<T, Container>(int1, int2, instep)
			{
					my_sys=&mighty;
					y=iny;
					num_func=iny.size();
					reset();
			};

			~ckrk(){	my_sys=NULL;	}

		/** Initial Condition Setters **/
			void set_y(const Container &in)
			{y=in; num_func=y.size();}
			void setY(const Container &in) { this->set_y(in);}
			void setInitialCondition(const Container &in) { this->set_y(in);}


		//this sucker is the meat of the matter
		//it does the acctuall R-K, the rest of the stuff here is
		//data handling and error checking...so this is important
			void rk();

		//here is the first round of error checking
		//it makes sure that everything is kosher within
		//one step, and changes the step size accordingly if it is not
			void round1();

		//finally the master checker and solver...
		//the second returns the output array
			void odeint();

};

END_BL_NAMESPACE

#include "driver/ckrk_meth.h"

#endif
