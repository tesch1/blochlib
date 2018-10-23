/* utils.h ********/


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
 	utils.h-->a collection of string and other utilities
 */

#ifndef _BL_utils_h_
#define _BL_utils_h_ 1

#include <string>

#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif

#if HAVE_UNISTD_H
	#include <unistd.h>
#endif

/*
#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# else
#  include <time.h>
# endif
*/
#  include <time.h>


#include "container/Vector/Vector.h"
#include "container/complex.h"
#include <iostream>
#include <stdio.h>
#include "utils/constants.h"


BEGIN_BL_NAMESPACE



void query_parameter(int argc,char* argv[],int par,const std::string& Q, char *V);
void query_parameter(int argc,char* argv[],int par,const std::string& Q,std::string& V);
void query_parameter(int argc,char* argv[],int par,const std::string& Q,double& V);
void query_parameter(int argc,char* argv[],int par,const std::string& Q,int& V);

//splits a single std::string into an array given a character to split around
Vector<std::string> parse_param(char * in);
Vector<std::string> parse_param(char * in, char what);
Vector<std::string> parse_param(std::string in);
Vector<std::string> parse_param(std::string in, char what);

void sttoch(std::string in, char *tobe);
std::string collapsVS(const Vector<std::string> &in);
std::string collapsVS(const Vector<std::string> &in, int lim2);
std::string collapsVS(const Vector<std::string> &in, int lim1, int lim2);

//the whitespace remover
std::string removeWhite(std::string);

//gets the stuff inside a set of '(' and ')'
std::string getInside(std::string in);

//parses /outermost set of commas
//i.e. 2, moo(4,5), 3
//will be parsed into::
// 2
// moo(4,5)
// 3
Vector<std::string> parseComma(const std::string &in);

//this strips out comments (lines starting with '#" and any substr starting with '#)
// and strips out any other text inbetween '{' and '}' sooo the input of
//
//  moo=34
//  loo=3 #moknkey
// #commet
// sec{
//	jj=0
// }

//will become

//  moo=34
//  loo=3
//
Vector<std::string> paramStrip(const Vector<std::string> &pset);

std::string itost(int i); //std::string from an int
std::string itost(int i, int len);//std::string from an int of length len
std::string itost_form(const std::string &fmt, int i);	//formatted std::string from an int
std::string dbtost_form(const std::string& fmt, double d); //formated std::string from a double
std::string dbtost(double d,const std::string& fmt="%.4f" ); //formated std::string from a double

//these happen to be already defined on the CodeWarrior ide
//so in case they are not, we define them here in our 'minmax' namespace

template<class T>
inline T sign(const T &in){
	if(in>=0) return T(1);
	return T(-1);
}
template<class T>
inline T max(const T &in,const T &in2){
	return (in)>(in2)?(in):(in2);
}
template<class T>
inline T max(const T &in){
	return in;
}
template<class T>
inline T min(const T &in,const T &in2){
	return (in)>(in2)?(in2):(in);
}
template<class T>
inline T min(const T &in){
	return in;
}


template<class T>
inline void swap_(T &in1, T &in2){
	T tmp=in1;
	in1=in2;
	in2=tmp;
}


template<class T>
inline T cube(const T &in)
{
	return in*in*in;
}


//. wrapers into the 'std' namespace
#ifdef BEGIN_BL_NAMESPACE
template<class T>	inline T sqrt(const T &in){	return std::sqrt(in);	}
template<class T>	inline T exp(const T &in){	return std::exp(in);	}
template<class T>	inline T log(const T &in){	return std::log(in);	}
template<class T>	inline T log10(const T &in){	return std::log10(in);	}
template<class T>	inline T sin(const T &in){	return std::sin(in);	}
template<class T>	inline T cos(const T &in){	return std::cos(in);	}
template<class T>	inline T tan(const T &in){	return std::tan(in);	}
template<class T>	inline T asin(const T &in){	return std::asin(in);	}
template<class T>	inline T acos(const T &in){	return std::acos(in);	}
template<class T>	inline T atan(const T &in){	return std::atan(in);	}
template<class T>	inline T sinh(const T &in){	return std::sinh(in);	}
template<class T>	inline T cosh(const T &in){	return std::cosh(in);	}
template<class T>	inline T tanh(const T &in){	return std::tanh(in);	}
template<class T>	inline T asinh(const T &in){	return std::asinh(in);	}
template<class T>	inline T acosh(const T &in){	return std::acosh(in);	}
template<class T>	inline T atanh(const T &in){	return std::atanh(in);	}
template<class T>	inline T ceil(const T &in){	return std::ceil(in);	}
template<class T>	inline T floor(const T &in){	return std::floor(in);	}
template<class T>	inline T pow(const T &in, T ra){	return std::pow(in, ra);	}
template<class T, class T2>	inline T pow(const T &in, T2 ra){	return std::pow(in, ra);	}
#endif

template<class T>
inline T factorial(const T &in)
{
	if(in<=1){	return 1;	}
	return factorial(in-1)*in;
}

template<class T>
inline T sum_factorial(T in)
{
	T ret=1;
	while(in>1){ ret+=factorial(in); --in;	}
	return ret;
}

inline int atoi(std::string in){	return atoi(in.c_str());	}
inline double atof(std::string in){	return atof(in.c_str());	}


void FFT1D_(complex data[], int size, int isign=1);
void IFFT1D_(complex data[], int size);
void FFTN_(complex data[], int dims[], int ndim, int isign=1);


/* binary i/o functions... */
bool BinaryWriteString(std::string out, std::fstream &oo);
std::string BinaryReadString(std::fstream &oo, int len);

/* typical usage of the timer class

timer stopwatch;

void timereset(){
   stopwatch.reset();
}

//a function that simply prints the time
//the 'times' variable is typically used when one repeats
//a particular operation MANY times to get a good timing sample
//thus the time for that 'one' operation is the 'grand time' divided by the
//number of applications

void print_usage(int times=1){
  std::cout << "Time taken: " << (1e6*stopwatch()/times) << " us\n";
}



*/



class timer {
  std::clock_t time_store;
  std::clock_t gettime() const;

  std::time_t low_time;

public:

  enum timer_type{LowRes, HighRes} ;

  timer_type resolution ;

  timer();
  timer(timer_type myRes);
  void reset();
  double operator()() const ;
};


END_BL_NAMESPACE



#endif

