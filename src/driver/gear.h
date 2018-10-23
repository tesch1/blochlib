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
	gear.h-->a 5th order Gear Diff eq solver
	  good for STIFF equations....

	  NOTE:: NEEDS WORK...this will not work quite yet....

*/

#ifndef _gear_h_
#define _gear_h_ 1

#include "container/containers.h"
#include "utils/utils.h"
#include "bs.h"


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
template< class func, class T=Vector<double> >
class gear {
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
		static double safty;
		static double err_max;
		static int max_step;
		static int max_step_within;
		static double tiny;
		static double eps;



		//these are for the integrator part...they are the various coefs needed
		//for gr part
		//the 'AA' are the global constants, the a0 types are element constants
		//we need ALL coefs for orders 1-->6
		static const double AA,a0,a1; //the first order gear coifs
		static const double BB, b0,b1, b2; //the second order order gear coifs
		static const double CC, c0,c1, c2, c3; //the third order order gear coifs
		static const double DD, d0,d1, d2, d3,d4; //the forth order order gear coifs
		static const double EE, e0, e1, e2, e3, e4, e5; //the fifth order order gear coifs
		static const double FF, f0, f1, f2, f3, f4, f5, f6; //the sixth order order gear coifs


		//some other variables we need
		//t1-> the start time
		//t2-> the end time
		//stepi->our initial step size
		//cur_step->current step size
		//stepdid->the one we acctualyy used
		//stepnext->estimated next step size
		//stepmin->the minium step size you wnat to do
		//cur_t->the current value of t we happen to be on
		double t1;
		double t2;
		double stepi;
		double cur_step;
		double stepdid;
		double stepnext;
		double stepmin;
		double cur_t;

		//here are some variables that tell you how many steps it took
		//ngood->good ones...ie. ones that just kept goin'
		//nbad->bads ones...ones that had to be adjusted..hopefully fixed

		int ngood;
		int nbad;


		//here are some arrays for saving data as you go if you'd like
		//
		int num_func;

		//the input arrays of equations (the initial conditions)
		//y->our variable list->initially ='y(0)', later 'yn'
		//yerr->error allowed for each of the variables will be |dy|+|step*fxs|
		//yscal->is what we compare our errors (yerr) to: =max_err*cur_step*fxs
		//out->is the tempary storage place for the R-K steps

		Vector<T> y;
		Vector<T> yerr;
		Vector<T> yscal;

		func *my_sys;

		bool hold_step;
		public:
	//defalt constructor, used mostly for testing and other classes if nessary
		gear(): t1(0), t2(0), stepi(0), cur_step(0), stepdid(0),
						stepnext(0), stepmin(0), cur_t(0), ngood(0),nbad(0), num_func(0) {}

		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instepmin->minimum step size to try
		//instep->initial step size guess
		//num->number of ODE/functions
		gear(Vector<T> &iny, double int1, double int2, double instep, double instepmin, func &mighty):
				t1(int1), t2(int2), stepi(instep),cur_step(instep), stepdid(0), stepnext(0),stepmin(instepmin),
				cur_t(int1),ngood(0), nbad(0)
		{
				if(cur_step > abs(int2-int1)/5.0) cur_step=abs(int2-int1)/5.0;
				my_sys=&mighty;
				y=iny;
				num_func=iny.size();
				set_rest_i();
		};
		//this one setup the basic inputs needed, all others take their default value
		//iny->inital y values
		//int1->start value
		//int2->end value
		//instep->initial step size guess
		//num->number of ODE/functions
		gear( Vector<T> &iny, double int1, double int2, double instep, func &mighty):
					t1(int1), t2(int2), stepi(instep), cur_step(instep),stepdid(0), stepnext(0),stepmin(0),
					cur_t(int1),ngood(0), nbad(0)
		{
				if(cur_step > abs(int2-int1)/5.0) cur_step=abs(int2-int1)/5.0;
				my_sys=&mighty;
				y=iny;
				num_func=iny.size();
				set_rest_i();
		};

		/*gear(_matrix<T, FullMatrix> &iny, T &int1, T &int2, T &instep, func &mighty):
			t1(int1), t2(int2), stepi(instep), cur_step(instep), stepdid(0), stepnext(0),stepmin(0),
			cur_t(int1),ngood(0), nbad(0){

				my_sys=&mighty;
				set_y(iny);
				num_func=iny.size();
				set_rest_i();
			}*/
		//some functions to step little bits in our class
		void set_time(double in1, double in2)	{t1=in1; t2=in2;}
		void set_time(double in1, double in2, double in3) { t1=in1; t2=in2; cur_t=t1; stepi=in3; cur_step=stepi;}
		void set_y(Vector<T> &in) 			{y=in; num_func=y.size();}
		//void set_y(matrix &in)							{y=set_y_mat(in);num_func=y.size();}
		//Vector<T> set_y_mat(matrix &iny);
		void set_yscal();
		void set_rest_i();
		void reset() { set_rest_i();}
		void set_stepi(double in)						{stepi=in;}
		void set_stepmin(double in)					{stepmin=in;}
		void set_cur_t(double in)						{cur_t=in;}
		void set_cur_step(double in)					{cur_step=in;}
		void set_stepnext(double in)					{stepnext=in;}
		void set_stepdid(double in)					{stepdid=in;}

		//this sucker is the meat of the matter
		//it does the acctuall Gear itterations,
		void gr(int stepct);

		//here is the first round of error checking
		//it makes sure that everything is kosher within
		//one step, and changes the step size accordingly if it is not

		void round1();
		//finally the master checker and solver...
		//the second returns the output array
		void odeint();
		//various functions to get our bits of stuff
		int get_ngood()									{return ngood;}
		int get_nbad()									{return nbad;}
		double get_cur_t() 							{return cur_t;}
		Vector<T> *get_out()				{return &y;}

		//hold-> tells the driver to limit step size increases to 10x the origina
		void hold() { hold_step=1; }
};

END_BL_NAMESPACE


#include "driver/gear_meth.h"

#endif

