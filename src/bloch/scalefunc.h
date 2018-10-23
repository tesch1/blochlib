/* scalefunc.h ********/


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
 	scalefunc.h-->scaling functions for the 'Dipolar Field'
 */

/* --the scale function 'scales' the caclulated Bd by some ammount as determined by the user
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
 * and take in the coord I the coord J, THe DISTANCE bewteen the
 * (or some norm you feel is good) and the current magnitization coord
 *
 * function(coord<> &r, coord<> &rp, doubele dist,  coord<> &Bd[i])
 *
 * Not all the argurments need to be used, this is just
 * to provide some uniformity to the resulting scaling functions
 *
 * */

#ifndef  _scalfunc_h_
#define  _scalfunc_h_

#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


 //a 'dummy' scale function the base class for all other functions
 //this function does nothing......


class NoScaleFunc
{
	private:

	public:
		NoScaleFunc(){}
		NoScaleFunc(const NoScaleFunc &in){}

		//THE MASTER FUNCTION!
		//here is does NOTHING
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd){}
};

std::ostream &operator<<(std::ostream &oo, const NoScaleFunc &out);



//this class is a "hard sphere" simply sets Bd to zero past a certain
//distance 'dist'

class HardSphere : public NoScaleFunc
{

	public:

		double cutoff;
		inline  HardSphere():
			cutoff(1.0e200) //some huge number so most all points are accepted
		{}

		inline HardSphere(const HardSphere &cp):
			NoScaleFunc(), cutoff(cp.cutoff)
		{}
		inline  HardSphere(double cp):
			cutoff(cp)
		{}

		inline HardSphere &operator =(const HardSphere &rhs)
		{
			if(this==&rhs) return *this;
			cutoff=rhs.cutoff;
			return *this;
		}

		void SetCutOff(double in){	cutoff=in;	}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			if(dist>cutoff) Bd=0;
		}

};

std::ostream &operator<<(std::ostream &oo, const HardSphere &out);

//this class is a "Inverse hard sphere" simply sets Bd to zero WITHIN a certain
//cutoff

class InvHardSphere : public NoScaleFunc
{

	public:
		double cutoff;
		inline  InvHardSphere():
			NoScaleFunc(),cutoff(0) //so all points are included
		{}

		inline InvHardSphere(const InvHardSphere &cp):
			NoScaleFunc(), cutoff(cp.cutoff)
		{}

		inline InvHardSphere(double cp):
			NoScaleFunc(), cutoff(cp)
		{}

		void SetCutOff(double in){	cutoff=in;	}

		inline InvHardSphere &operator =(const InvHardSphere &rhs)
		{
			if(this==&rhs) return *this;
			cutoff=rhs.cutoff;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			if(dist<cutoff) Bd=0;
		}

};

std::ostream &operator<<(std::ostream &oo, const InvHardSphere &out);



//this class is a hard shell simply sets Bd to if the
//distance is OUTside the shell

class HardShell : public NoScaleFunc
{

	public:
		double r1;
		double r2;
		inline  HardShell():
			r1(0), r2(1.0e200) //so all points are included
		{}

		inline HardShell(const HardShell &cp):
			NoScaleFunc(), r1(cp.r1), r2(cp.r2)
		{}

		void SetCutOff(double in, double in2){	r1=in; r2=r2;	}

		inline HardShell &operator =(const HardShell &rhs)
		{
			if(this==&rhs) return *this;
			r1=rhs.r1;
			r2=rhs.r2;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			if(dist<r1 || dist>=r2) Bd=0;
		}

};


std::ostream &operator<<(std::ostream &oo, const HardShell &out);




//this class is a Inverse hard shell simply sets Bd to if the
//distance is INSIDE the shell

class InvHardShell : public NoScaleFunc
{

	public:
		double r1;
		double r2;
		inline  InvHardShell():
			r1(0), r2(1.0e200) //so all points are included
		{}

		inline InvHardShell(const InvHardShell &cp):
			NoScaleFunc(), r1(cp.r1), r2(cp.r2)
		{}

		void SetCutOff(double in, double in2){	r1=in; r2=r2;	}

		inline InvHardShell &operator =(const InvHardShell &rhs)
		{
			if(this==&rhs) return *this;
			r1=rhs.r1;
			r2=rhs.r2;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			if(dist>r1 || dist<=r2) Bd=0;
		}

};


std::ostream &operator<<(std::ostream &oo, const InvHardShell &out);

//this class is a inverse roll off scaling
// given a distance 'r' it scales the values by
// 'factor*r^pow'

class Power : public NoScaleFunc
{

	public:
		double power;
		double factor;
		inline  Power():
			power(1), factor(1) //so all points are not scaled
		{}

		inline Power(const Power &cp):
			NoScaleFunc(), power(cp.power), factor(cp.factor)
		{}

		inline Power(double pow, double fact):
			NoScaleFunc(), power(pow), factor(fact)
		{}

		void SetParams(double power_, double factor_){	power=power_; factor=factor_;	}

		inline Power &operator =(const Power &rhs)
		{
			if(this==&rhs) return *this;
			power=rhs.power;
			factor=rhs.factor;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			//cout<<"PRE:"<<factor*std::pow(dist, power)/dist<<" "<<Bd<<endl;
			Bd*=factor*pow(dist, power)/dist;
			//cout<<"POST:" <<factor*std::pow(dist, power)/dist<<" "<<Bd<<endl;

		}

};


std::ostream &operator<<(std::ostream &oo, const Power &out);
//this class scales by SUBTRACTION and a power law
// so the output is Bd-factor*r^pow

class SubtractPower : public NoScaleFunc
{

	public:
		double power;
		double factor;
		inline  SubtractPower():
			power(1), factor(1) //so all points are not scaled
		{}

		inline SubtractPower(const SubtractPower &cp):
			NoScaleFunc(), power(cp.power), factor(cp.factor)
		{}

		inline SubtractPower(const Power &cp):
			NoScaleFunc(), power(cp.power), factor(cp.factor)
		{}

		inline SubtractPower(double pow, double fact):
			NoScaleFunc(), power(pow), factor(fact)
		{}

		void SetParams(double power_, double factor_){	power=power_; factor=factor_;	}

		inline SubtractPower &operator =(const SubtractPower &rhs)
		{
			if(this==&rhs) return *this;
			power=rhs.power;
			factor=rhs.factor;
			return *this;
		}

		inline SubtractPower &operator =(const Power &rhs)
		{
			power=rhs.power;
			factor=rhs.factor;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			Bd-=factor*pow(dist, power);
		}

};


std::ostream &operator<<(std::ostream &oo, const SubtractPower &out);

// A Bolztman distribution scaling function
// it gives a scal of 0 at dist=0, increases to a max of ~0.54 scaling
// at dist=2*factor, then decreases to zero at dist-->infinity
//
// uses a blotzaman veclocity distribution function
// the 'factor' should be LESS THAN 1
// factor^2 Exp[factor*dist]*dist^2

class BoltzmannScale : public NoScaleFunc
{

	public:
		double factor;
		inline  BoltzmannScale():
			factor(1) //so all points are not scaled
		{}

		inline BoltzmannScale(const BoltzmannScale &cp):
			NoScaleFunc(), factor(cp.factor)
		{}

		inline BoltzmannScale(const double &cp):
			NoScaleFunc(), factor(cp)
		{}

		void SetParams(double factor_){ factor=factor_;	}

		inline BoltzmannScale &operator =(const BoltzmannScale &rhs)
		{
			if(this==&rhs) return *this;
			factor=rhs.factor;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
			Bd*=(factor*factor)*exp(factor*dist)*dist*dist;
		}

};


std::ostream &operator<<(std::ostream &oo, const BoltzmannScale &out);


// A hyperbolic tanget scaleing
//  scale factor starts at 0 then increases to 1 (an stays there until infinity)
//  the rise time can be controled by the factor
//  factors LESS then 1 slow the rise, factors BIGGER then 1 seed the rise
//
//  a factor of 1 gives the scaling of 0.5 at ~x=0.55
// scale=tanh(x*factor)

class TanhScale : public NoScaleFunc
{

	public:
		double factor;
		inline  TanhScale():
			factor(1) //so all points are not scaled
		{}

		inline TanhScale(const TanhScale &cp):
			NoScaleFunc(), factor(cp.factor)
		{}

		inline TanhScale(const double &cp):
			NoScaleFunc(), factor(cp)
		{}


		void SetParams(double factor_){ factor=factor_;	}

		inline TanhScale &operator =(const TanhScale &rhs)
		{
			if(this==&rhs) return *this;
			factor=rhs.factor;
			return *this;
		}

		//THE MASTER FUNCTIONS
		inline void function(coord<> &R, coord<> &Rp,double dist, coord<> &Bd)
		{
//			cout<<"PRE: "<<std::tanh(factor*dist)<<endl;
			Bd*=tanh(factor*dist);
		}

};


std::ostream &operator<<(std::ostream &oo, const TanhScale &out);


END_BL_NAMESPACE


#endif



