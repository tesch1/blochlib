


/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06.20.02
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
	basicODE.h-->contains the basic data elements for the ODE solvers
*/

#ifndef _basicODE_h_
#define _basicODE_h_ 1


#include "container/containers.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


template<class T=double, class Container=Vector<T> >
class basicODE
{
	public:

	/*** "Error" and Stepper Values ***/
	//safty --> a basic '10% drop in any reduced value
	//safty1->  a 75% safty becuase we do not know our error exactly...we pretend it is smaller
	//safty2->    30%     "
	//red_max->maximum reduction of a step size
	//red_min->minium rediuction of a step size
	//		to keep moving the step size
	//tiny->incase we get stuck around '0' add tiny to keep things goin'
	//scal_max->1/scal_max is the maximum stepsize increase
	//eps->our accuracy goal
	//max_step->the number of steps to be taken
	//relToler->the _relative_ tolerence of errors
	//absToler->the absolute tolerence of errors

		static const double safty;
		static const double safty1;
		static const double safty2;
		static const double red_max;
		static const double red_min;
		static const double tiny;
		static const double scal_max;
		double eps;
		int max_step;
		double relToler;
		double absToler;


	/*** TIME STEP VARIABLES ***/
	//t1-> the start time
	//t2-> the end time
	//stepi->our initial step size
	//cur_step->current step size
	//abs_step->The ABSOLUTE value of the current step size
	//stepdid->the one we acctualyy used
	//stepnext->estimated next step size
	//stepmin->the minium step size you wnat to do
	//cur_t->the current value of t we happen to be on
	//t_direction--> the sign of the time (either -1 (backwards) or +1 (fowards))
	//hold--tells the driver to limit step size increases to 5x the origina

		double t1;
		double t2;
		double stepi;
		double cur_step;
		double abs_step;
		double stepdid;
		double stepnext;
		double stepmin;
		double stepmax;
		double cur_t;
		double tdirection;
		bool hold_step;

	//some stats variables
	//ngood->good ones...ie. ones that just kept goin'
	//nbad->bads ones...ones that had to be adjusted..hopefully fixed
	//func_calls --> number of function calls
	//LU_calls--> number of LU decomp class (if valid)
	//jac_calls--> number of Jacobians calced
		int ngood;
		int nbad;
		int func_calls;
		int LU_calls;
		int jac_calls;

	//here are some arrays for saving data as you go if you'd like
		int num_func;

	//the input arrays of equations (the initial conditions)
	//y->our variable list->initially ='y(t=0)', later 'yn'
	//yerr->error allowed for each of the variables will be |dy|+|step*fxs|
	//yscal->is what we compare our errors (yerr) to: =max_err*cur_step*fxs
	//ysav-> a temporary storing spot for silly thngs
	//ytmp-> a temporary storing spot for silly thngs
		Container y;
		Container yerr;
		Container yscal;
		Container ysav;
		Container ytmp;

		//the FUNCTIONS ptrs typdefs..
		//NOTE:: these are just the type definitions
		//this 'registers' the generation function with a string
		typedef void (*function_t)(double time,Container &y, Container &dydt) ;
		typedef _matrix<typename ObjectTrait<T>::numtype, FullMatrix> mattype;
		typedef void (*jacobian_t)(double time,Container &y, mattype &dfdy) ;



		basicODE():
			eps(1.0e-6),max_step(100000000),relToler(1.0e-3),absToler(1.0e-6),
			t1(0.0), t2(0.0), stepi(0.0), cur_step(0.0), abs_step(0.0), stepdid(0.0),
			stepnext(0.0),stepmin(0.0),stepmax(0.0),cur_t(0.0), tdirection(1.0),hold_step(false),
			ngood(0), nbad(0), func_calls(0),LU_calls(0),jac_calls(0), num_func(0)
		{}

		basicODE(double start, double endt):
			eps(1.0e-6), max_step(100000000),relToler(1.0e-3),absToler(1.0e-6),
			t1(start), t2(endt), stepi(0.0), cur_step(0.0),
			abs_step(0.0), stepdid(0.0),
			stepnext(0.0),stepmin(0.0),stepmax(0.1*(endt-start)),
			cur_t(start),tdirection(double(sign(endt-start))),hold_step(false),
			ngood(0), nbad(0), func_calls(0),LU_calls(0),jac_calls(0), num_func(0)

		{}


		basicODE(double start, double endt,	double instep, double instepmin=0.0, double epsin=1.0e-6, int maxstep=100000000):
			eps(epsin), max_step(maxstep),relToler(1.0e-3),absToler(1.0e-6),
			t1(start), t2(endt), stepi(instep), cur_step(instep),
			abs_step(abs(instep)),stepdid(0.0),
			stepnext(0.0),stepmin(instepmin),stepmax(0.1*(endt-start)),
			cur_t(start), tdirection(double(sign(endt-start))),hold_step(false),
			ngood(0), nbad(0), func_calls(0),LU_calls(0),jac_calls(0), num_func(0)

		{}


		virtual ~basicODE(){}

	/*** Error Functions ***/
		void NoFuncError(std::ostream &out, std::string header="")
		{
			out<<header<<std::endl;
			out<<"ERROR::NO Function Class defined"<<std::endl;
			out<<" cannot integrate anything...."<<std::endl;
			BLEXCEPTION(" cannot integrate anything....")
		}

	/*** Initial Condition Setters ***/
		virtual void set_y(const Container &in);
		virtual void setY(const Container &in) { set_y(in);}
		virtual void setInitialCondition(const Container &in) { set_y(in);}

	/*** Time Setters **/

		//set begining and end times
		void set_time(double in1, double in2)
		{
			t1=in1; t2=in2; cur_t=t1;
			tdirection=double(sign(t2-t1));

		}
		void setTime(double in1, double in2) { set_time(in1, in2);	}

		//set time and initial step size guess
		void set_time(double t1in, double t2in, double stepin){
			t1=t1in;t2=t2in;stepi=stepin;cur_t=t1;cur_step=stepin;
			if(stepin>abs(t1in-t2in)){
				stepi=(t2in-t1in)/10.0;
				cur_step=stepi;
			}

			if(stepin==0.0){
				stepi=(t2in-t1in)/10.0;
				cur_step=stepi;
			}
		}
		void setTime(double t1in, double t2in, double stepin)
		{	set_time(t1in, t2in, stepin);	}

	/*** the 'reset' function ***/
		void reset();

	/** StepSize 'holder ***/
		void hold() {	hold_step=true;	}


	/** the 'solved' data (assuming the integrator was called)	***/
		Container *get_out()	{	return &y;	}
		Container *solvedData()	{	return get_out();	}

	/** Data Dumpers for errors or just info **/
		void error(std::ostream &out, std::string message="");
		void print(std::ostream &out, std::string message="");
		void info(std::ostream &out, std::string message="")
		{	print(out, message);	}

	/** the "integrators" these MUST be defined in the 'sub classes' **/
		virtual void odeint()
		{	NoFuncError(std::cerr, "cannot call odeint from 'basicODE'");	}

		virtual void solve()
		{ odeint();	}

		virtual void solve(double t1, double t2, double step=0.0)
		{ setTime(t1, t2, step); odeint();	}



	//the 'Vector' solver
		virtual Vector<Container> solve(Vector<double> times, double initStep=0.0);

	//the 'Vector' solver with initial conditions y
		virtual Vector<Container> solve(Vector<double> times,const Container &intc, double initStep=0.0);
};



/*** Error And Stepper Constants ***/

template<class T, class Container >
const double  basicODE<T,Container >::safty=0.9;

template<class T, class Container >
const double  basicODE<T,Container >::safty1=0.25;

template<class T, class Container >
const double  basicODE<T,Container >::safty2=0.7;

template<class T, class Container >
const double  basicODE<T,Container >::red_max=1.0e-5;

template<class T, class Container >
const double  basicODE<T,Container >::red_min=0.7;

template<class T, class Container >
const double  basicODE<T,Container >::tiny=1.0e-30;

template<class T, class Container >
const double  basicODE<T,Container >::scal_max=0.1;


/*** Initial Condition Setters ***/
//the 'generic MASK set up'
template<class T, class Container >
void basicODE<T,Container >::set_y(const Container &in)
{	y=in;	num_func=y.size();	}


//the 'generic type' for T=double or complex, etc
template<class T, class Container >
void basicODE<T,Container >::reset()
{
	ysav=y; ysav.fill(0);
	ytmp=ysav;
	yerr=ysav;
	yscal=ysav;
	stepnext=0.;
	stepdid=0.;
	stepmin=0.;
	hold_step=true;
}

/*//the 'generic type' for T=double or complex, etc
template<class T >
void basicODE<T,coord<T> >::reset()
{
	ysav=y; ysav.fill(0);
	ytmp=ysav;
	yerr=ysav;
	yscal=ysav;
	stepnext=0.;
	stepdid=0.;
	stepmin=0.;
	hold_step=true;
}*/
//the reset for T=coord<T1, N>


template<class T, class Container >
void basicODE<T,Container >::error(std::ostream &out, std::string message)
{
	if(message!="")	out<<message<<std::endl;
	out<<"Number of functions: "<<num_func<<std::endl;
	out<<"Last Step Size: "<<stepdid<<std::endl;
	out<<"Next Step Size: "<<stepnext<<std::endl;
	out<<"Minimum Step Size: "<<stepmin<<std::endl;
	out<<"Max Step Size: "<<stepmax<<std::endl;
	out<<"Time direction: "<<tdirection<<std::endl;
	out<<"With Current Step size: "<<cur_step<<std::endl;
	out<<"On Current Time: "<<cur_t<<std::endl;
	print(out);
}

template<class T, class Container >
void basicODE<T,Container >::print(std::ostream &out, std::string message)
{
	if(message!="")	out<<message<<std::endl;
	out<<"The Start Time: "<<t1<<std::endl;
	out<<"The End Time: "<<t2<<std::endl;
	out<<"The max error: "<<eps<<std::endl;
	out<<"The Absolute tolerence: "<<absToler<<std::endl;
	out<<"The Relative tolerence: "<<relToler<<std::endl;
	out<<"Number of 'good' integrator steps: "<<ngood<<std::endl;
	out<<"Number of 'bad' integrator steps: "<<nbad<<std::endl;
	out<<"Number of function calls: "<<func_calls<<std::endl;
	if(LU_calls) out<<"Number of LU decomp calls: "<<LU_calls<<std::endl;
	if(jac_calls) out<<"Number of Jacobians calls: "<<jac_calls<<std::endl;
}

template<class T, class Container >
Vector<Container> basicODE<T,Container >::solve(Vector<double> times, double initStep)
{
	Vector<Container> out(times.size(), 0);
	if(times.size()<1)
	{
		BLEXCEPTION(" the Time Vector must be of length 1 of more...")
	}

	//what to do if there is only one point in the times list
	if(times.size()==1)
	{
		if(times[0]==cur_t){
			std::cerr<<std::endl<<" WARNING: bs::solve(Vector<double>) "<<std::endl;
			std::cerr<<" final time is the same as curent time..doing nothing"<<std::endl;
			out[0]=*(get_out());
			return out;
		}
		setTime(cur_t, times[0], initStep);
		odeint();
		out[0]=*(get_out());
		return out;
	}

	//grab the t=0 point (or where ever we happen to be at cur_t)
	// if that point is desired
	int onpt=0, datapt=0;
	if(times[0]==cur_t)
	{
		out[0]=*(get_out());
		datapt++;
	}

	//set the initial t1 and t2 using the initStep only once
	//so the initgrator determins the stepsize from here on out
	setTime(times[0], times[1], initStep);
	while(onpt<times.size()-1)
	{
		setTime(times[onpt], times[onpt+1]);
		odeint();
		out[datapt++]=*(get_out());
		onpt++;
	}
	return out;
}

template<class T, class Container >
Vector<Container> basicODE<T,Container >::solve(Vector<double> times, const Container &intc, double initStep)
{
	setInitialCondition(intc);
	return solve(times, initStep);
}



END_BL_NAMESPACE


#endif
