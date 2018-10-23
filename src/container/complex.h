

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
	complex.h-->complex numbers...easy to understand i  would think....
*/


#ifndef BL_complex_h_
#define BL_complex_h_ 1

//#include <ostream>
//#include <iomanip>
#include <math.h>
#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "utils/blassert.h"
#include "container/rankType.h"

#ifdef ON_WINDOWS
	#include "blochconfigCW.h"
#else
	#include "blochconfig.h"
	#ifdef __MINGW32__
		#undef HAVE_ISNAN
	#endif
#endif

/* functions inside STD */
BEGIN_BL_NAMESPACE
template<class Ctype_T>
class Complex;

typedef Complex<double> complex; //Ctype_T prec
typedef Complex<float> scomplex; //single prec
typedef Complex<int> icomplex; //int prec


END_BL_NAMESPACE

namespace std{
template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::sqrt(const BlochLib::Complex<Ctype_T> & z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::exp(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::log(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::log10(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::pow(const BlochLib::Complex<Ctype_T>& z, const BlochLib::Complex<Ctype_T>& z1);  //z^z1

template<class Ctype_T, class T2>
inline BlochLib::Complex<Ctype_T> std::pow(const BlochLib::Complex<Ctype_T>& z, const T2& z1);  //z^z1

//triggggsss
template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::sin(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::cos(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::tan(const BlochLib::Complex<Ctype_T>& z);

//'a' triggggssss
template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::asin(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::acos(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::atan(const BlochLib::Complex<Ctype_T>& z);

//'h' triggggssss
template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::sinh(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::cosh(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::tanh(const BlochLib::Complex<Ctype_T>& z);

//'a''h' trigggssssss
template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::asinh(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::acosh(const BlochLib::Complex<Ctype_T>& z);

template<class Ctype_T>
inline BlochLib::Complex<Ctype_T> std::atanh(const BlochLib::Complex<Ctype_T>& z);
};

BEGIN_BL_NAMESPACE

//*****************CONJUGATES, RE, IM*******************
// this allows use to overload the Complex<Ctype_T> 'conj' function
// why even bother with a conjugate for REALNAME numbers, you may ask?
// if we have a templated container (a vector, matrix, etc) and wish to take
// the 'conjugate' of it...well suppose the object inside the container is
// an 'int' we still want 'conj' to return something and not give an error uring compilation

//absolute value template function
template<class T>
T abs(T in){
	if(in<T(0)) return -in;
	return in;
}

#ifdef HAVE_FINITE
inline bool hasnan(const double &z){	return !finite((double)z);	}
inline bool hasnan(const float &z){	return !finite((float)z);	}
inline bool hasnan(const int &z){	return !finite((float)z);	}
inline bool hasnan(const short &z){	return !finite((float)z);	}
inline bool hasnan(const char &z){	return !finite((float)z);	}
inline bool hasnan(const bool &z){	return !finite((float)z);	}
#elif HAVE_ISNAN
inline bool hasnan(const double &z){	return isnan((double)z);	}
inline bool hasnan(const float &z){	return isnan((float)z);	}
inline bool hasnan(const int &z){	return isnan((float)z);	}
inline bool hasnan(const short &z){	return isnan((float)z);	}
inline bool hasnan(const char &z){	return isnan((float)z);	}
inline bool hasnan(const bool &z){	return isnan((float)z);	}
#endif




#define compatCom(TYPE) \
	inline TYPE conj(const TYPE &z){ return z;	}	\
	inline TYPE norm(const TYPE &z){ return abs(z);	}		\
	inline TYPE AbsNorm(const TYPE &z){ return abs(z);	}		\
	inline TYPE Re(const TYPE &z){ return z;	}			\
	inline void Re( TYPE &z, TYPE in){ z=in;	}			\
	inline TYPE Im(const TYPE &z){ return TYPE(0);	}			\
	inline TYPE square_norm(const TYPE& z){ return z*z; }		\
	inline void Im(TYPE &z, TYPE in){ return; }			\
	inline TYPE max(const TYPE &z){ return z; }			\
	inline TYPE min(const TYPE &z){ return z; }			\
	inline void mul(TYPE &z, TYPE in1, TYPE in2){ z=in1*in2;	}		\
inline TYPE chop(const TYPE &z, float eps=1.e-12){	return TYPE(std::fabs(float(z)))>eps?z:TYPE(0);	} \
inline TYPE chop(const TYPE &z, double eps=1.e-12){	return TYPE(std::fabs(double(z)))>eps?z:TYPE(0);	} \
	inline TYPE chop(const TYPE &z, int eps=0){	return TYPE(fabs(float(z)))>eps?z:TYPE(0);	} \

	//inline void chop(TYPE &z, Ctype_T eps=1.e-12){	z=abs(z)>eps?z:TYPE(0);	}



compatCom(double);
compatCom(float);
compatCom(int);
compatCom(long);
compatCom(char);

//this is a strange function i know, but it keeps certain 'matrix<string>' functions
// from barfing in Borland C++ Builder
inline std::string Re(const std::string &in){	return in;	}
inline std::string Im(const std::string &in){	return in;	}

//nessesary becuase some code wants
// 'REALNAME', 'IMAGNAME' as the vars,
//some code wnats 'r', 'i' as the var
// and some want 'REALNAME', 'IMAGNAME'
//depending on the BLAS lib uses

#define real REALNAME
#define imag IMAGNAME

//this little guy is a nice formatting trick which
// is based on the gamma lib
/*
I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/
*/
extern const std::string cmx_form;

template<class Ctype_T>
class Complex{
	private:

		void error(int which) const;

	public:

//out bits of data....
		Ctype_T REALNAME, IMAGNAME;

		typedef Ctype_T numtype;

//constructors
		inline Complex<Ctype_T>(){};
		inline Complex<Ctype_T>(Ctype_T r, Ctype_T i=0.0):
		  REALNAME(r), IMAGNAME(i){};

		template<class T1, class T2>
		inline Complex<Ctype_T>(T1 r, T2 i=0.0):
		  REALNAME(r), IMAGNAME(i){};

		inline Complex<Ctype_T>(const Complex<Ctype_T>& z):
		  REALNAME(z.REALNAME), IMAGNAME(z.IMAGNAME){};

		template<class T2>
		inline Complex<Ctype_T>(const Complex<T2>& z):
		  REALNAME(z.REALNAME), IMAGNAME(z.IMAGNAME){};

//assignments

		template<class T2>
		inline Complex<Ctype_T>& operator= (const Complex<T2>& z){
			REALNAME = z.REALNAME;
			IMAGNAME = z.IMAGNAME;
			return *this;
		}

		template<class T2>
		inline Complex<Ctype_T>& operator= (T2 r){
			REALNAME = r;
  			IMAGNAME = 0;
  			return *this;
		}

//comparisons
#if 0
		 inline bool operator== (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline bool operator!= (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline bool operator!= (const Complex<Ctype_T>& z1, const Ctype_T& z2);
#endif
		 template<class T2>
		 inline bool operator== (Complex<T2> z2) const;
		 template<class T2>
		 inline bool operator== (T2 z2) const;
		 template<class T2>
		 inline bool operator!= (Complex<T2> z2) const;
		 template<class T2>
		 inline bool operator!= (T2 z2) const;

		inline bool operator<(const Complex<Ctype_T>& z) const;
		inline bool operator>(const Complex<Ctype_T>& z) const;
		inline bool operator<=(const Complex<Ctype_T>& z) const;
		inline bool operator>=(const Complex<Ctype_T>& z) const;

//addition operators
		 template<class T2>
		 inline void operator+= (const Complex<T2>& z);
		 template<class T2>
		 inline void operator+= (const T2& z);


#if 0
		 inline void add(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void add(Complex<Ctype_T>& z,       Ctype_T    r, const Complex<Ctype_T>& z1);
		 inline void add(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1,       Ctype_T    r);
 		 inline Complex<Ctype_T> operator+ (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline Complex<Ctype_T> operator+ (const Complex<Ctype_T>& z, Ctype_T r);
		 inline void operator+= (Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void operator+= (Complex<Ctype_T>& z, Ctype_T r);
#endif

//subtraction operators
		template<class T2>
		 inline void operator-= (const Complex<T2>& z);
		template<class T2>
		 inline void operator-= (const T2& z);

#if 0
		 inline void sub(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void sub(Complex<Ctype_T>& z,       Ctype_T    r, const Complex<Ctype_T>& z1);
		 inline void sub(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1,       Ctype_T    r);
		 inline Complex<Ctype_T> operator- (const Complex<Ctype_T>& z, Ctype_T r);
		 inline Complex<Ctype_T> operator- (Ctype_T r, const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> operator- (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void operator-= (Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void operator-= (Complex<Ctype_T>& z, Ctype_T r);
#endif

//mulitplication ops
		template<class T2>
		 inline void operator*= (const Complex<T2>& z);
		template<class T2>
		 inline void operator*= (const T2& z);


#if 0
		 inline void mul(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void mul(Complex<Ctype_T>& z,       Ctype_T    r, const Complex<Ctype_T>& z1);
		 inline void mul(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1,       Ctype_T    r);
		 inline Complex<Ctype_T> operator* (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline Complex<Ctype_T> operator* (Ctype_T r, const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> operator* (const Complex<Ctype_T>& z, Ctype_T r);
		 inline void operator*= (Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void operator*= (Complex<Ctype_T>& z, Ctype_T r);
#endif


//divisions
		template<class T2>
		 inline void operator/= (const Complex<T2>& z);
		template<class T2>
		 inline void operator/= (const T2& z);

#if 0
		 inline void div(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline void div(Complex<Ctype_T>& z,       Ctype_T    r, const Complex<Ctype_T>& z1);
		 inline void div(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1,       Ctype_T    r);
		 inline Complex<Ctype_T> operator/ (const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2);
		 inline Complex<Ctype_T> operator/ (Ctype_T r, const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> operator/ (const Complex<Ctype_T>& z, Ctype_T r);
		 inline void operator/= (Complex<Ctype_T>& z, const Complex<Ctype_T>& z1);
		 inline void operator/= (Complex<Ctype_T>& z, Ctype_T r);
#endif



//other Complex<Ctype_T> type functions
		 inline Complex<Ctype_T> conj(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> conj(const Complex<Ctype_T>& a, const Complex<Ctype_T>& b);
		 inline Complex<Ctype_T> conj(const Ctype_T a, const Complex<Ctype_T>& b);
		 inline Complex<Ctype_T> conj(const Complex<Ctype_T>& a, const Ctype_T b);

//multi operation 'quick' function...save a few CPU clicks
		//a+=b*conj(c)
		 inline void muladd_conj(Complex<Ctype_T> &a,const Complex<Ctype_T> &b,const Complex<Ctype_T> &c);
		//a+=conj(b)*c
		 inline void conj_muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> &b,const Complex<Ctype_T> &c);
		//a+=b*c
		 inline void muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> & b,const Complex<Ctype_T> & c);
		//a+=b*c...b is REALNAME
		 inline void muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> &c,Ctype_T b);
       	//a+=b*c...c is REALNAME
      	 inline void muladd(Complex<Ctype_T> &a,Ctype_T c,const Complex<Ctype_T> &b);

		inline Ctype_T Re()const{	return REALNAME;	}
		inline Ctype_T Im()const{ return IMAGNAME;	}
		inline Ctype_T &Re(){	return REALNAME;	}
		inline Ctype_T &Im(){ return IMAGNAME;	}
		template<class T1>
		inline void Re(T1 r){ REALNAME=r; };
		template<class T1>
		inline void Im(T1 r){ IMAGNAME=r; }

		 inline void set_real_part(Complex<Ctype_T>& z, Ctype_T r){ z.REALNAME=r;	}
		 inline void set_imaginary_part(Complex<Ctype_T>& z, Ctype_T r){ z.IMAGNAME=r;	}
#if 0
		 inline Ctype_T real(const Complex<Ctype_T>& z){ return z.REALNAME;	}
		 inline Ctype_T imag(const Complex<Ctype_T>& z){ return z.IMAGNAME;	}
		inline Ctype_T imag()const{ return IMAGNAME;	}
		inline Ctype_T real()const{ return REALNAME;	}

		inline Ctype_T &imag(){ return IMAGNAME;	}
		inline Ctype_T &real(){ return REALNAME;	}

		inline void imag(Ctype_T r){ IMAGNAME=r;	}
		inline void real(Ctype_T r){ REALNAME=r;	}
#endif
		 inline Complex<Ctype_T> chop(const Complex<Ctype_T> &in, Ctype_T eps=1e-12);

		 inline Ctype_T AbsNorm(const Complex<Ctype_T>& z);
		 inline Ctype_T square_norm(const Complex<Ctype_T>& z);
		 inline Ctype_T norm(const Complex<Ctype_T>& z);
		 inline void norm(Complex<Ctype_T>& z, Ctype_T r);  //sets the norm of z to r


		 inline Ctype_T phase(const Complex<Ctype_T>& z);
		 inline void phase(Complex<Ctype_T>& z, Ctype_T r);
	 	inline Complex<Ctype_T> Zexp() const;
#if 0
		 inline Complex<Ctype_T> std::sqrt(const Complex<Ctype_T> & z);
		 inline Complex<Ctype_T> std::exp(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::log(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::log10(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::pow(const Complex<Ctype_T>& z, const Complex<Ctype_T>& z1);  //z^z1

//triggggsss
		 inline Complex<Ctype_T> std::sin(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::cos(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::tan(const Complex<Ctype_T>& z);

//'a' triggggssss
		 Complex<Ctype_T> std::asin(const Complex<Ctype_T>& z);
		 Complex<Ctype_T> std::acos(const Complex<Ctype_T>& z);
		 Complex<Ctype_T> std::atan(const Complex<Ctype_T>& z);

//'h' triggggssss
		 inline Complex<Ctype_T> std::sinh(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::cosh(const Complex<Ctype_T>& z);
		 inline Complex<Ctype_T> std::tanh(const Complex<Ctype_T>& z);

//'a''h' trigggssssss
		 Complex<Ctype_T> std::asinh(const Complex<Ctype_T>& z);
		 Complex<Ctype_T> std::acosh(const Complex<Ctype_T>& z);
		 Complex<Ctype_T> std::atanh(const Complex<Ctype_T>& z);
#endif

//io stuff (ascii)
#if 0
		 std::ostream& operator<< (std::ostream& ostr, const Complex<Ctype_T>& z);		//out
		 std::istream& operator>> (std::istream& istr, Complex<Ctype_T>& z);		//in
#endif
//io stuff (binary)
		void write(const std::string& fn);
		void write(std::ofstream& fp) const;

		void read(const std::string& fn);
		void read(std::ifstream& fp);

};


//defined in Complex.cc
extern const Complex<double> complex0;		//  (0,0)
extern const Complex<double> complex1;		//  (1,0)
extern const Complex<double> complexi;		//  (0,1)
//defined in Complex.cc
extern const Complex<float> scomplex0;		//  (0,0)
extern const Complex<float> scomplex1;		//  (1,0)
extern const Complex<float> scomplexi;		//  (0,1)


/*********** RE and IM ****/
template<class Ctype_T>
 inline Ctype_T Re(const Complex<Ctype_T>& z){ return z.REALNAME;	}
template<class Ctype_T>
 inline Ctype_T Im(const Complex<Ctype_T>& z){ return z.IMAGNAME;	}

template<class Ctype_T, class T1>
inline void Im(Complex<Ctype_T> &z, T1 r){ z.IMAGNAME=r; }
template<class Ctype_T, class T1>
inline void Re(Complex<Ctype_T> &z, T1 r){ z.REALNAME=r; }


//*****************************
//most of our functions are inlined friends...
//so they must be defined outside the class, yet in the same file



//additions******************************************************


template<class Ctype_T>
inline void add(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2){
	z.REALNAME = z1.REALNAME + z2.REALNAME;
	z.IMAGNAME = z1.IMAGNAME + z2.IMAGNAME;
}


template<class Ctype_T>
inline void add(Complex<Ctype_T>& z,Ctype_T r, const Complex<Ctype_T>& z1){
	z.REALNAME = r + z1.REALNAME;
	z.IMAGNAME = z1.IMAGNAME;
}


template<class T1,class T2>
inline void add(Complex<T1>& z, const Complex<T1>& z1, T2   r){
	z.REALNAME = z1.REALNAME + r;
	z.IMAGNAME = z1.IMAGNAME;
}

template<class T2>
//dummy operator does this  '+z'
inline Complex<T2> operator+ (const Complex<T2>& z1){
	return Complex<T2>(z1);

}

template<class T1, class T2>
inline Complex<T1> operator+ (const Complex<T1>& z1, const Complex<T2>& z2){
  return Complex<T1>
      (z1.REALNAME+z2.REALNAME, z1.IMAGNAME+z2.IMAGNAME);
}


template<class T1>
inline Complex<T1> operator+ (const Complex<T1>& z1, float r){
	return Complex<T1>(r+z1.REALNAME, z1.IMAGNAME);
}

template<class T1>
inline Complex<T1> operator+ (float r, const Complex<T1>& z1){
	return Complex<T1>(r+z1.REALNAME, z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator+ (const Complex<T1>& z1, double r){
	return Complex<double>(r+z1.REALNAME, z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator+ (double r, const Complex<T1>& z1){
	return Complex<double>(r+z1.REALNAME, z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator+ (const Complex<T1>& z1, int r){
	return Complex<double>(r+z1.REALNAME, z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator+ (int r, const Complex<T1>& z1){
	return Complex<double>(r+z1.REALNAME, z1.IMAGNAME);
}


template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator+= (const Complex<T2>& z1){
	REALNAME += z1.REALNAME;
	IMAGNAME += z1.IMAGNAME;
}

template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator+= (const T2& z1){
	REALNAME += z1;
}


//subtraction***************************************************

template<class Ctype_T>
inline void sub(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2){
	z.REALNAME = z1.REALNAME - z2.REALNAME;
	z.IMAGNAME = z1.IMAGNAME - z2.IMAGNAME;
}

template<class Ctype_T>
inline void sub(Complex<Ctype_T>& z,Ctype_T r,const Complex<Ctype_T>& z1){
	z.REALNAME = r - z1.REALNAME;
	z.IMAGNAME = z1.IMAGNAME;
}

template<class T1, class T2>
inline void sub(Complex<T1>& z, const Complex<T1>& z1,T2 r){
	z.REALNAME = z1.REALNAME - r;
	z.IMAGNAME = z1.IMAGNAME;
}

template<class Ctype_T>
inline Complex<Ctype_T> operator- (const Complex<Ctype_T>& z1){
	return Complex<Ctype_T>(-z1.REALNAME, -z1.IMAGNAME );
}

template<class T1>
inline Complex<double> operator- (const Complex<T1>& z1, double r){
	return Complex<double>(z1.REALNAME-r, z1.IMAGNAME);
}


template<class T1>
inline Complex<double> operator- (double r, const Complex<T1>& z1){
	return Complex<double>(r-z1.REALNAME, -z1.IMAGNAME);
}

template<class T1>
inline Complex<T1> operator- (const Complex<T1>& z1, float r){
	return Complex<T1>(z1.REALNAME-r, z1.IMAGNAME);
}


template<class T1>
inline Complex<T1> operator- (float r, const Complex<T1>& z1){
	return Complex<T1>(r-z1.REALNAME, -z1.IMAGNAME);
}

template<class T1>
inline Complex<T1> operator- (const Complex<T1>& z1, int r){
	return Complex<T1>(z1.REALNAME-r, z1.IMAGNAME);
}
template<class T1>
inline Complex<T1> operator- (int r, const Complex<T1>& z1){
	return Complex<T1>(r-z1.REALNAME, -z1.IMAGNAME);
}

template<class T1, class T2>
inline Complex<T1> operator- (const Complex<T1>& z1, const Complex<T2>& z2){
	return Complex<T1>(z1.REALNAME-z2.REALNAME, z1.IMAGNAME-z2.IMAGNAME);
}


template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator-= (const Complex<T2>& z2){
  REALNAME -= z2.REALNAME;
  IMAGNAME -= z2.IMAGNAME;
}


template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator-= (const T2 &r)
{ REALNAME -= r; }

//multiplication************************************************

template<class Ctype_T,class T2, class T1>
inline void mul(Complex<Ctype_T>& z, const Complex<T1>& z1, const Complex<T2>& b){
	z.REALNAME = z1.REALNAME*b.REALNAME - z1.IMAGNAME*b.IMAGNAME;
	z.IMAGNAME = z1.REALNAME*b.IMAGNAME + z1.IMAGNAME*b.REALNAME;
}

template<class Ctype_T,class T2, class T1>
inline void mul(Complex<Ctype_T>& z,T2 r,const Complex<T1>& z1){
	z.REALNAME = r * z1.REALNAME;
	z.IMAGNAME = r * z1.IMAGNAME;
}

template<class Ctype_T,class T2, class T1>
inline void mul(Complex<Ctype_T>& z, const Complex<T1>& z1,T2 r){
	z.REALNAME = z1.REALNAME * r;
	z.IMAGNAME = z1.IMAGNAME * r;
}

template<class T1, class T2>
inline Complex<T1> operator* (const Complex<T1>& z1, const Complex<T2>& z2){
	return Complex<T1>
	   (z1.REALNAME*z2.REALNAME - z1.IMAGNAME*z2.IMAGNAME,
	    z1.REALNAME*z2.IMAGNAME + z1.IMAGNAME*z2.REALNAME);
}

template<class T1>
inline Complex<double> operator* (int r, const Complex<T1>& z1){
	return Complex<double>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}


template<class T1>
inline Complex<double> operator* (double r, const Complex<T1>& z1){
	return Complex<double>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}

template<class T1>
inline Complex<T1> operator* (float r, const Complex<T1>& z1){
	return Complex<T1>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator* (const Complex<T1>& z1, int r){
	return Complex<double>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}

template<class T1>
inline Complex<double> operator* (const Complex<T1>& z1, double r){
	return Complex<double>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}

template<class T1>
inline Complex<T1> operator* (const Complex<T1>& z1, float r){
	return Complex<T1>
	  (r*z1.REALNAME, r*z1.IMAGNAME);
}


template<class Ctype_T>
template<class T1>
inline void Complex<Ctype_T>::operator*= (const Complex<T1>& z2){
	T1 r;
	r     = REALNAME*z2.REALNAME - IMAGNAME*z2.IMAGNAME;
	IMAGNAME = REALNAME*z2.IMAGNAME + IMAGNAME*z2.REALNAME;
	REALNAME = r;
}


template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator*= (const T2 &r){
	REALNAME = REALNAME*r;
	IMAGNAME = IMAGNAME*r;
}

//DIVISIONS****************************************************
template<class Ctype_T>
inline void div(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2){
	Ctype_T x = z2.REALNAME*z2.REALNAME + z2.IMAGNAME*z2.IMAGNAME;
	z.REALNAME = (z1.REALNAME*z2.REALNAME + z1.IMAGNAME*z2.IMAGNAME)/x;
	z.IMAGNAME = (z2.REALNAME*z1.IMAGNAME - z2.IMAGNAME*z1.REALNAME)/x;
}

template<class T2, class Ctype_T>
inline void div(Complex<Ctype_T>& z,T2 r, const Complex<Ctype_T>& b){
	Ctype_T x = b.REALNAME*b.REALNAME + b.IMAGNAME*b.IMAGNAME;
	z.REALNAME = (r * b.REALNAME)/x;
	z.IMAGNAME = (-b.IMAGNAME * r)/x;
}


template<class T2, class Ctype_T>
inline void div(Complex<Ctype_T>& z, const Complex<Ctype_T>& z1,T2 r){
	z.REALNAME = z1.REALNAME/r;
	z.IMAGNAME = z1.IMAGNAME/r;
}

template<class T1, class T2>
 inline Complex<T1> operator/
   (const Complex<T1>& z1, const Complex<T2>& z2){
	Complex<T1> z;
	T1 r = z2.REALNAME*z2.REALNAME + z2.IMAGNAME*z2.IMAGNAME;
	z.REALNAME = (z1.REALNAME*z2.REALNAME + z1.IMAGNAME*z2.IMAGNAME)/r;
	z.IMAGNAME = (z2.REALNAME*z1.IMAGNAME - z2.IMAGNAME*z1.REALNAME)/r;
	return z;
}


template<class T1>
inline Complex<T1> operator/ (int r, const Complex<T1>& z1){
	Complex<T1> z;
	T1 d = z1.REALNAME*z1.REALNAME + z1.IMAGNAME*z1.IMAGNAME;
	z.REALNAME = (r*z1.REALNAME)/d;
	z.IMAGNAME = (-z1.IMAGNAME*r)/d;
	return z;
}


template<class T1>
inline Complex<T1> operator/ (const Complex<T1>& z1, int r){
	return Complex<T1>(z1.REALNAME/r, z1.IMAGNAME/r);
}

template<class T1>
inline Complex<T1> operator/ (float r, const Complex<T1>& z1){
	Complex<T1> z;
	T1 d = z1.REALNAME*z1.REALNAME + z1.IMAGNAME*z1.IMAGNAME;
	z.REALNAME = (r*z1.REALNAME)/d;
	z.IMAGNAME = (-z1.IMAGNAME*r)/d;
	return z;
}


template<class T1>
inline Complex<T1> operator/ (const Complex<T1>& z1, float r){
	return Complex<T1>(z1.REALNAME/r, z1.IMAGNAME/r);
}

template<class T1>
inline Complex<double> operator/ (double r, const Complex<T1>& z1){
	Complex<T1> z;
	T1 d = z1.REALNAME*z1.REALNAME + z1.IMAGNAME*z1.IMAGNAME;
	z.REALNAME = (r*z1.REALNAME)/d;
	z.IMAGNAME = (-z1.IMAGNAME*r)/d;
	return z;
}


template<class T1>
inline Complex<double> operator/ (const Complex<T1>& z1, double r){
	return Complex<T1>(z1.REALNAME/r, z1.IMAGNAME/r);
}

template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator/= (const Complex<T2>& z2){
	Ctype_T r = z2.REALNAME*z2.REALNAME+z2.IMAGNAME*z2.IMAGNAME;
	Ctype_T REALNAME = (REALNAME*z2.REALNAME+IMAGNAME*z2.IMAGNAME)/r;
	IMAGNAME = (z2.REALNAME*IMAGNAME-z2.IMAGNAME*REALNAME)/r;
	REALNAME = REALNAME;
}

template<class Ctype_T>
template<class T2>
inline void Complex<Ctype_T>::operator/= (const T2 &r){
	REALNAME = REALNAME/r;
	IMAGNAME = IMAGNAME/r;
}

//*******Comparisons**********
//equality type operators
template<class Ctype_T>
template<class T2>
inline bool Complex<Ctype_T>::operator== (Complex<T2> z1) const{
	return ((REALNAME==z1.REALNAME) && (IMAGNAME==z1.IMAGNAME));
}

template<class Ctype_T>
template<class T2>
inline bool Complex<Ctype_T>::operator== (T2 z1) const {
	return ((REALNAME==z1) && (IMAGNAME==0.0));
}

template<class Ctype_T>
template<class T2>
inline bool Complex<Ctype_T>::operator!= (Complex<T2> z1)const {
	return ((REALNAME!=z1.REALNAME) || (IMAGNAME!=z1.IMAGNAME));
}

template<class Ctype_T>
template<class T2>
inline bool Complex<Ctype_T>::operator!= (T2 z1) const{
	return ((REALNAME!=z1) || (IMAGNAME!=0.0));
}

/*
//equality type operators
template<class Ctype_T>
inline bool operator== (const Complex<Ctype_T>& z, const Complex<Ctype_T>& z1){
	return ((z.REALNAME==z1.REALNAME) && (z.IMAGNAME==z1.IMAGNAME));
}

template<class Ctype_T>
inline bool operator!= (const Complex<Ctype_T>& z, const Complex<Ctype_T>& z1){
	return ((z.REALNAME!=z1.REALNAME) || (z.IMAGNAME!=z1.IMAGNAME));
}


template<class Ctype_T, class T2>
inline bool operator!= (const Complex<Ctype_T>& z, const T2 &z1){
	return ((z.REALNAME!=z1) || (z.IMAGNAME!=z1));
}
*/
template<class Ctype_T>
inline bool Complex<Ctype_T>::operator<(const Complex<Ctype_T>& z) const
{ return (BL_NAMESPACE::norm(*this) < BL_NAMESPACE::norm(z)); }

template<class Ctype_T>
inline bool Complex<Ctype_T>::operator>(const Complex<Ctype_T>& z) const
{ return (BL_NAMESPACE::norm(*this) > BL_NAMESPACE::norm(z)); }

template<class Ctype_T>
inline bool Complex<Ctype_T>::operator<=(const Complex<Ctype_T>& z) const
  { return (BL_NAMESPACE::norm(*this) <= BL_NAMESPACE::norm(z)); }

template<class Ctype_T>
inline bool Complex<Ctype_T>::operator>=(const Complex<Ctype_T>& z) const
  { return (BL_NAMESPACE::norm(*this) >= BL_NAMESPACE::norm(z)); }


//Random optimization functions Complex<Ctype_T> functions
template<class Ctype_T>
inline Complex<Ctype_T> conj(const Complex<Ctype_T>& z1){
	return Complex<Ctype_T>(z1.REALNAME , - z1.IMAGNAME );
}

template<class Ctype_T>
inline Complex<Ctype_T> conj(const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2){
	return Complex<Ctype_T>( z1.REALNAME*z2.REALNAME+z1.IMAGNAME*z2.IMAGNAME, z1.REALNAME*z2.IMAGNAME-z1.IMAGNAME*z2.REALNAME);
}

template<class Ctype_T>
inline Complex<Ctype_T> conj_mul(const Complex<Ctype_T>& z1, const Complex<Ctype_T>& z2){
	return conj(z1,z2);
}

template<class T,class Ctype_T>
inline Complex<Ctype_T> conj(T z1, const Complex<Ctype_T>& z2){
	return Complex<Ctype_T>( z1*z2.REALNAME, -z1*z2.IMAGNAME);
}

template<class T,class Ctype_T>
inline Complex<Ctype_T> conj(const Complex<Ctype_T>& z1, T z2){
	return z2*z1;
}


/*
inline Complex<Ctype_T> conj(const Complex<Ctype_T>& z1, const Ctype_T z2){
	return Complex<Ctype_T>( z1.REALNAME*z2, -z1.IMAGNAME*z2);
}


inline Complex<Ctype_T> conj(const Ctype_T z1, const Complex<Ctype_T> &z2){
	return Complex<Ctype_T>( z2.REALNAME*z1, z2.IMAGNAME*z1);
}
*/
template<class Ctype_T>
inline void muladd_conj(Complex<Ctype_T> &a,const Complex<Ctype_T> &b,const Complex<Ctype_T> &c){
	a.REALNAME+=b.REALNAME*c.REALNAME+b.IMAGNAME*c.IMAGNAME;
	a.IMAGNAME+=b.IMAGNAME*c.REALNAME-b.REALNAME*c.IMAGNAME;
}

template<class Ctype_T>
inline void conj_muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> &b,const Complex<Ctype_T> &c){
	a.REALNAME+=b.REALNAME*c.REALNAME+b.IMAGNAME*c.IMAGNAME;
	a.IMAGNAME+=b.REALNAME*c.IMAGNAME-b.IMAGNAME*c.REALNAME;
}

template<class Ctype_T>
inline void muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> &c,Ctype_T b){
	a.REALNAME+=b*c.REALNAME;
	a.IMAGNAME+=b*c.IMAGNAME;
}

template<class Ctype_T>
inline void muladd(Complex<Ctype_T> &a,Ctype_T b,const Complex<Ctype_T> &c){
	a.REALNAME+=b*c.REALNAME;
	a.IMAGNAME+=b*c.IMAGNAME;
}


template<class Ctype_T>
inline void muladd(Complex<Ctype_T> &a,const Complex<Ctype_T> & b,const Complex<Ctype_T> & c){
	a.REALNAME+=b.REALNAME*c.REALNAME-b.IMAGNAME*c.IMAGNAME;
	a.IMAGNAME+=b.IMAGNAME*c.REALNAME+b.REALNAME*c.IMAGNAME;
}

//**********88these complete template set...since we may
//wish to uses the function with doubles or ints, etc

template<class T>
inline T conj(const T z,const T z1){  return z*z1;	}


template<class T>
inline void muladd_conj(T &a,const T &b,const T &c){
	a+=conj(b)*c;
}

template<class T>
inline void conj_muladd(T &a,const T &b,const T &c){
	a+=b*conj(c);
}

template<class T, class T1, class T2>
inline void muladd(T &a, T1 c,T2 b){
	a+=b*c;
}
//*************


//keeps this function valid for REALNAME numbers too
//template<class T, class T1>
//inline void set_realName_part(T& z,      T1 r) { z=r; }

template<class Ctype_T>
inline Ctype_T AbsNorm(const Complex<Ctype_T>& z)
{ return fabs(z.REALNAME) + fabs(z.IMAGNAME); }

template<class Ctype_T>
inline Ctype_T square_norm(const Complex<Ctype_T>& z)
{ return z.REALNAME*z.REALNAME + z.IMAGNAME*z.IMAGNAME; }

template<class Ctype_T>
inline Ctype_T norm(const Complex<Ctype_T>& z){
	return std::sqrt(z.REALNAME*z.REALNAME+z.IMAGNAME*z.IMAGNAME);
	/*if(z.REALNAME>=z.IMAGNAME){
		return abs(z.REALNAME)*std::sqrt(1.+(z.IMAGNAME/z.REALNAME)*(z.IMAGNAME/z.REALNAME));
	}else{
		return abs(z.IMAGNAME)*std::sqrt(1.+(z.REALNAME/z.IMAGNAME)*(z.REALNAME/z.IMAGNAME));
	}*/
}


template<class Ctype_T>
inline Ctype_T abs(const Complex<Ctype_T> &z)
{ return norm(z);	}

template<class Ctype_T>
inline void norm(Complex<Ctype_T>& z, Ctype_T r){
	r /= std::sqrt(z.REALNAME*z.REALNAME + z.IMAGNAME*z.IMAGNAME );
	z.REALNAME *= r;
	z.IMAGNAME *= r;
}

template<class Ctype_T>
inline Ctype_T phase(const Complex<Ctype_T>& z){
	if(z==0.0) return 0;
	return atan2(z.IMAGNAME, z.REALNAME);
}

template<class Ctype_T, class T1>
inline Complex<Ctype_T> chop(const Complex<Ctype_T> &in, T1 eps)
{	return Complex<Ctype_T>(abs(in.REALNAME)<eps?0.0:in.REALNAME, abs(in.IMAGNAME)<eps?0.0:in.IMAGNAME);	}


template<class Ctype_T>
inline void phase(Complex<Ctype_T>& z, Ctype_T r){
	Ctype_T tmp = norm(z);
	z.REALNAME = std::cos(r) * tmp;
	z.IMAGNAME = std::sin(r) * tmp;
}

template<class Ctype_T>
inline Complex<Ctype_T> Complex<Ctype_T>::Zexp() const
{
	Complex<Ctype_T> z;
	Ctype_T e = std::exp(this->REALNAME);
	z.REALNAME = e * std::cos(IMAGNAME);
	z.IMAGNAME = e * std::sin(IMAGNAME);
	return z;
}

END_BL_NAMESPACE

// FUNCTIONS INSIDE STD NAMESPACE
namespace std{

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::sqrt(const BL_NAMESPACE::Complex<Ctype_T>& z1)
{
	 BL_NAMESPACE::Complex<Ctype_T> z;
	  Ctype_T r = BL_NAMESPACE::norm(z1);
	  if(r==0.0) z = 0;
	  else
	    {
	    Ctype_T p = BL_NAMESPACE::phase(z1)/2;
	    r = std::sqrt(r);
	    z.REALNAME = r*std::cos(p);
	    z.IMAGNAME = r*std::sin(p);
	    }

	return z;

}


template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::exp(const BL_NAMESPACE::Complex<Ctype_T>& z1){
	BL_NAMESPACE::Complex<Ctype_T> z;
	Ctype_T e = std::exp(z1.REALNAME);
	z.REALNAME = e * std::cos(z1.IMAGNAME);
	z.IMAGNAME = e * std::sin(z1.IMAGNAME);
	return z;
}



template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::log(const BL_NAMESPACE::Complex<Ctype_T>& z)
{
	return BL_NAMESPACE::Complex<Ctype_T>(std::log(norm(z)), BL_NAMESPACE::phase(z));
}

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::log10(const BL_NAMESPACE::Complex<Ctype_T>& z)
{
	return BL_NAMESPACE::Complex<Ctype_T>(std::log10(norm(z)), BL_NAMESPACE::phase(z));
}

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::pow(const BL_NAMESPACE::Complex<Ctype_T>& z, const BL_NAMESPACE::Complex<Ctype_T>& z1){
	return BL_NAMESPACE::Complex<Ctype_T>(std::exp(z1 * log(z)));
}

template<class T2, class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::pow(const BL_NAMESPACE::Complex<Ctype_T>& z, const T2& z1)
{
	return BL_NAMESPACE::Complex<Ctype_T>( std::pow(z.Re(), z1), std::pow(z.Im(), z1));
}


//TRIGGGGSSSS**************************************************
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::sin(const BL_NAMESPACE::Complex<Ctype_T>& z)
{
	Ctype_T a = std::exp(z.IMAGNAME);
	Ctype_T b = 1/a;
	return BL_NAMESPACE::Complex<Ctype_T>(std::sin(z.REALNAME) * (a+b)/2, std::cos(z.REALNAME) * (a-b)/2);
}

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::cos(const BL_NAMESPACE::Complex<Ctype_T>& z)
{
	Ctype_T a = std::exp( z.IMAGNAME );
	Ctype_T b = 1/a;
	return BL_NAMESPACE::Complex<Ctype_T>(std::cos(z.REALNAME) * (a+b)/2, -std::sin(z.REALNAME) * (a-b)/2 );
}

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::tan(const BL_NAMESPACE::Complex<Ctype_T>& z )
{
	Ctype_T a = std::exp(2 * z.IMAGNAME );
	Ctype_T b = 1/a;
	Ctype_T d = std::cos(2 * z.REALNAME ) + (a+b)/2;
	return BL_NAMESPACE::Complex<Ctype_T>(std::sin(2*z.REALNAME) / d, (a-b)/(2*d));
}


//'h' TRIGGGGGSSSS*********************************************
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::sinh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return BL_NAMESPACE::Complex<Ctype_T>( std::sinh(z.REALNAME)*std::cos(z.IMAGNAME),std::cosh(z.REALNAME)*std::sin(z.IMAGNAME) ); }

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::cosh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return BL_NAMESPACE::Complex<Ctype_T>( std::cosh(z.REALNAME)*std::cos(z.IMAGNAME),std::sinh(z.REALNAME)*std::sin(z.IMAGNAME) ); }

template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::tanh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{
  Ctype_T tmp = std::cos(2*z.REALNAME) + std::cosh(2*z.IMAGNAME);
  return BL_NAMESPACE::Complex<Ctype_T>( std::sin(2*z.REALNAME)/tmp, std::sinh(2*z.IMAGNAME)/tmp );
}

// 'a' TRIGGGGGSSSS***********************************************
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::asin(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return - BL_NAMESPACE::complexi*std::log(BL_NAMESPACE::complexi*z + sqrt (1-z*z)); }
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::acos(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return BL_NAMESPACE::Complex<Ctype_T>(0,-1)*std::log(z + sqrt (z*z-1)); }
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::atan(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return BL_NAMESPACE::Complex<Ctype_T>(0,-0.5)*std::log((1+BL_NAMESPACE::complexi*z)/(1-BL_NAMESPACE::complexi*z)); }

// 'a''h' TRIGGGGGSSSS***********************************************
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::asinh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return std::log(z + std::sqrt(z*z + 1)); }
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::acosh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return std::log(z + std::sqrt(z*z - 1)); }
template<class Ctype_T>
inline BL_NAMESPACE::Complex<Ctype_T> std::atanh(const BL_NAMESPACE::Complex<Ctype_T>& z)
{ return std::log((1+z)/(1-z))/2; }

};

BEGIN_BL_NAMESPACE

///***************goood oll iioo******************(likely 'o')

template<class Ctype_T>
std::ostream& operator << (std::ostream& otr, const Complex<Ctype_T>& z);

/**** SEVERAL OTHER METHODS ***/

//errors
template<class Ctype_T>
void Complex<Ctype_T>::error(int which) const  {
  switch (which){
    case 0: BLEXCEPTION("...peace be with you...") break;
    case 1: BLEXCEPTION(" evil file troubles...") break;
    case 4: BLEXCEPTION(" cannot read from the file") break;
    case 6: BLEXCEPTION("cannot write to file") break;
    case 7:  BLEXCEPTION("output failure") break;
    default:   BLEXCEPTION("unknown error") break;
  }
}


//reads every two numbers in as num1=REALNAME, num2=Complex<Ctype_T> (NO little charater 'i's)
template<class Ctype_T>
std::istream& operator >> (std::istream& istr, Complex<Ctype_T>& z)
{return istr >> z.REALNAME >> z.IMAGNAME;}

template<class Ctype_T>
void Complex<Ctype_T>::write(const std::string& fn)  //create binary out stream then write
{	
	std::ofstream fp;
	fp.open(fn.c_str(),std::ios::out|std::ios::binary);
	if(!fp) BLEXCEPTION(" evil file troubles...");
	write(fp);
	fp.flush();
	fp.close();
}


template<class Ctype_T>
void Complex<Ctype_T>::write(std::ofstream& fp) const //binary write to out file stream
{
	fp<<&REALNAME;
	fp<<&IMAGNAME;
}


template<class Ctype_T>
void Complex<Ctype_T>::read(const std::string& fn)  //create the stream, then read
{	
	std::ifstream fp;
	fp.open(fn.c_str(),std::ios::in|std::ios::binary);
	if(!fp)	BLEXCEPTION(" evil file troubles...");
	read(fp);
	fp.close();
}


template<class Ctype_T>
void Complex<Ctype_T>::read(std::ifstream& fp)  //readin from binary IN file stream
{	
	if(!fp)  BLEXCEPTION(" evil file troubles...");
	if(fp.eof())	BLEXCEPTION(" evil file troubles...");
	fp>>REALNAME;
	if(fp.eof())	BLEXCEPTION(" evil file troubles...");
	fp>>IMAGNAME;
}

template<class Ctype_T>
std::ostream& operator << (std::ostream& otr, const Complex<Ctype_T>& z){
	char buff1[80], buff2[80];
	sprintf(buff1, cmx_form.c_str(),z.REALNAME);
  	std::string rr=std::string(buff1);
  	sprintf(buff2, cmx_form.c_str(),z.IMAGNAME);
  	std::string ii=std::string(buff2);
	otr<<"("<<rr<<","<<ii<<")";

	return otr;
}

END_BL_NAMESPACE

#endif								// Complex<Ctype_T>.h
