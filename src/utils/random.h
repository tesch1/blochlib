/* random.h ********/


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
 	random.h-->a collection of random number generators and classes
 */


#ifndef _bloch_random_h_
#define _bloch_random_h_

#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif

#if HAVE_UNISTD_H
#include <unistd.h>
#endif


#include "container/complex.h"

#include <iostream>
#include <time.h>

BEGIN_BL_NAMESPACE


//using namespace std;

/*********************/
//the set of functions below are to be used as single C type calls
// they are later incorperated into the classes...

//a random function (produces good random numbers regardless of seed)
// between 0..1
double Rand();

//gaussian distributed random number (good regardless of seed)
double GaussianRand();

//poisson distributed random number
double GammaLog(double xx); //needed for the PoissonRandom()
double PoissonRand(double mean);


//THis is the real random number generator
class BL_Random{
	private:
		static time_t rand_timer; //for an initial seed
		static const unsigned int N=624;
		static const unsigned int M=397;

		/* constant vector a */
		static const unsigned long MATRIX_A=0x9908b0df;
		/* most significant w-r bits */
		static const unsigned long UPPER_MASK=0x80000000;
		static const unsigned long LOWER_MASK=0x7fffffff;
	/* the array for the state vector  */
		static unsigned long mt[N];
		static int mti;
		static const unsigned long TEMPERING_MASK_B=0x9d2c5680;
		static const unsigned long TEMPERING_MASK_C=0xefc60000;

	public:

	//seed with a timer
		BL_Random();
	//seed CANNOT BE 0
		static void seed(unsigned long seedr);
	//integrer random numbers
		static unsigned long randomI();
	//double random generator
		static double randomD();
};


/**********************/
//The MAIN base class for all the random generators...DOES NOT generate any
// random numbmers...just holds the max and mins
template<class Num_t=double>
class BasicRandom
{
	private:
		Num_t range_;    //the magnitude of the random diviation i.e. range_*Random()
		Num_t offset_;  //the 'center point' for the random number ... offset_+range_*Random()

	public:
		typedef Num_t numtype;

		BasicRandom():
			range_(1), offset_(0)
		{}

		BasicRandom(Num_t low, Num_t high):
			range_(high-low), offset_(low)
		{}

		BasicRandom(const BasicRandom &cp):
			range_(cp.range_), offset_(cp.offset_)
		{}

		inline Num_t &low(){	return offset_;	}
		inline Num_t &range(){	return range_;	}
		inline Num_t high(){	return offset_+range_;	}
		inline Num_t &offset(){	return offset_;	}

		inline void low(Num_t lo){	 offset_=lo; range_=high()-offset_+1;	}
		inline void range(Num_t ra){	 range_=ra;	}
		inline void high(Num_t hi){	 range_=hi-offset_+1;	}
		inline void offset(Num_t off){	 offset_=off;	}

		inline void setLow(Num_t in){	offset_=in;	range_=high()-in+1;	}
		inline void setHigh(Num_t in){	range_=in-offset_+1;	}
		inline void set(Num_t lo, Num_t hi)
		{	range_=hi-lo, offset_=lo;	}

		inline void setOffset(Num_t in){	offset_=in;	}
		inline void setRange(Num_t in){	range_=in;	}


		BasicRandom &operator=(BasicRandom &rhs)
		{
			if(this==&rhs) return *this;
			range_=rhs.range_;
			offset_=rhs.offset_;
			return *this;
		}

};

//special case for 'complex'


//basic UNIFORM random number engine
// uses the 'std::Random() as its main driver

template<class Num_t=double>
class UniformRandom : public BasicRandom<Num_t>
{

	public:

		UniformRandom():
			BasicRandom<Num_t>()
		{}

		UniformRandom(Num_t low, Num_t high):
			BasicRandom<Num_t>(low,high)
		{}

		UniformRandom(const UniformRandom &cp):
			BasicRandom<Num_t>(cp)
		{}

		inline Num_t Random()
		{
			return Num_t(offset()+range()*Rand());
		}

		inline operator Num_t()
		{
			return Random();
		}
};

template<>
class UniformRandom<complex> : public BasicRandom<complex>
{

	public:

		UniformRandom():
			BasicRandom<complex>()
		{}

		UniformRandom(complex low, complex high):
			BasicRandom<complex>(low,high)
		{}

		UniformRandom(const UniformRandom &cp):
			BasicRandom<complex>(cp)
		{}

		inline complex Random()
		{
			return complex(Re(offset())+Re(range())*Rand(), Im(offset())+Im(range())*Rand());
		}

		inline operator complex()
		{
			return Random();
		}
};

template<class Num_t>
std::ostream &operator<<(std::ostream &oo,UniformRandom<Num_t> &out)
{
	oo<<"Uniform Random Distribution: "<<out.low()<<".."<<out.high();
	return oo;
}


//Poisson random number engine
// uses the 'std::PoissonRand() as its main driver

//the template here is not needed and WILL BE ERRORED
// if you use something other then 'double'
template<class Num_t=double>
class PoissonRandom : public BasicRandom<Num_t>
{
	private:
		Num_t mean_;
	public:

		PoissonRandom():
			BasicRandom<Num_t>(),
			mean_(1)
		{}

		PoissonRandom(Num_t mean):
			BasicRandom<Num_t>(),
			mean_(mean)
		{}

		PoissonRandom(Num_t low, Num_t high):
			BasicRandom<Num_t>(low,high),
			mean_(1)
		{}

		PoissonRandom(Num_t low, Num_t high, Num_t mean):
			BasicRandom<Num_t>(low,high),
			mean_(mean)
		{}

		PoissonRandom(const PoissonRandom &cp):
			BasicRandom<Num_t>(cp), mean_(cp.mean_)
		{}

		inline void SetMean(Num_t in){	mean_=in;	}
		inline Num_t &Mean(){	return mean_;	}

		PoissonRandom &operator=(PoissonRandom &rhs)
		{
			if(this==&rhs) return *this;
			BasicRandom<Num_t>::operator=(rhs);
			mean_=rhs.mean_;
			return *this;
		}

 		inline Num_t Random()
		{
			return Num_t(offset()+range()*std::PoissonRand(mean_));
		}

		inline operator Num_t()
		{
			return PoissonRandom::Random();
		}
};

template<class Num_t>
std::ostream &operator<<(std::ostream &oo, PoissonRandom<Num_t> &out)
{
	oo<<"Poisson Random Distribution: "<<out.low()<<".."<<out.high();
	return oo;
}

//Gaussian random number engine
// uses the 'std::PoissonRand() as its main driver

//the template here is not needed and WILL BE ERRORED
// if you use something other then 'double'
template<class Num_t=double>
class GaussianRandom : public BasicRandom<Num_t>
{

	public:

		GaussianRandom():
			BasicRandom<Num_t>()
		{}

		GaussianRandom(Num_t low, Num_t high):
			BasicRandom<Num_t>(low,high)
		{}


		GaussianRandom(const GaussianRandom &cp):
			BasicRandom<Num_t>(cp)
		{}


		GaussianRandom &operator=(GaussianRandom &rhs)
		{
			if(this==&rhs) return *this;
			BasicRandom<Num_t>::operator=(rhs);
			return *this;
		}


 		inline Num_t Random()
		{
			return Num_t(low()+high()*std::GaussianRand());
		}

		inline operator Num_t()
		{
			return GaussianRandom::Random();
		}
};

template<class Num_t>
std::ostream &operator<<(std::ostream &oo, GaussianRandom<Num_t> &out)
{
	oo<<"Gaussian Random Distribution: center="<<out.low()<<"..varience="<<out.high();
	return oo;
}


//this is a prototype Random class such that
//it can be used to produce a Vector of random numbers
//or a matrix of random numbers
//the two MOST important functions here are the
//--> T operator(T, int i) --> for vectors
//--> T operator(T, int i, int j) --> for matrix

//The Random Engines classes are defined above.

template<class Engine_t=UniformRandom<double> >
class Random : public Engine_t{

	public:
		typedef typename Engine_t::numtype numtype;
		Random() :
			Engine_t()
		{}

		Random(numtype mean):
			Engine_t(mean)
		{}

		Random(numtype low, numtype high):
			Engine_t(low,high)
		{}

		Random(numtype low, numtype high, numtype mean):
			Engine_t(low,high, mean)
		{}


		Random(const Random &cp):
			Engine_t(cp)
		{}


		inline operator  numtype()
		{
			return Engine_t::Random();
		}

		inline numtype operator()()
		{
			return Engine_t::Random();
		}


		template<class othernumT>
		numtype operator()(othernumT in, int i)
		{
			return (Engine_t::Random());
		}

		template<class othernumT>
		numtype operator()(othernumT in, int i, int j)
		{
			return Engine_t::Random();
		}
};




END_BL_NAMESPACE





#endif


