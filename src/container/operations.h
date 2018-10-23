/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-28-01
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

/******
*	operations.h--> contain the expersion template operations classes
*  Here we define the appilcative template 'operations' that are use in both
*  The Vector and Matrix class
******/


#ifndef _AllOps_h_
#define _AllOps_h_ 1

#ifndef ON_WINDOWS
#include "blochconfig.h"
#endif

#include <math.h>

BEGIN_BL_NAMESPACE

/**************************************************/
//Vector Assignment....
/*************************************************/
// a 'quick' meta program (one the compiler performs) to unroll loops completely...
template<int N, int I>
class vecAssign {
public:
    enum { loopFlag = (I < N-1) ? 1 : 0 };

    template<class VecType, class Expr, class Op>
    static inline void assign(VecType& vec, Expr expr, Op u)
    {
        u.apply(vec[I], expr(I));
        vecAssign<N * loopFlag, (I+1) * loopFlag>::assign(vec,expr,u);
    }

    template<class VecType, class T, class Op>
    static inline void assignWithArgs(VecType& vec, Op u,
        T x0, T x1=0, T x2=0, T x3=0,
        T x4=0, T x5=0, T x6=0, T x7=0,
        T x8=0, T x9=0)
    {
        u.apply(vec[I], x0);
        vecAssign<N * loopFlag, (I+1) * loopFlag>::assignWithArgs(vec, u, x1, x2, x3, x4, x5, x6, x7, x8, x9);
    }

};

//the class to 'kill' or stop the above one...we get here we stop the unrolling
template<>
class vecAssign<0,0> {
public:
    template<class VecType, class Expr, class Op>
    static inline void assign(VecType& vec, Expr expr, Op u)
    { }

    template<class VecType, class T, class Op>
    static inline void assignWithArgs(VecType& vec, Op u,
        T x0, T x1=0, T x2=0, T x3=0,
        T x4=0, T x5=0, T x6=0, T x7=0,
        T x8=0, T x9=0)
    {}
};

//****************************************************************************
// Double operators i.e. a+b, a-b, etc applicative templates --
//****************************************************************************

#define MakeDoubleOp(name, op) \
	template<class T1,class T2> \
	class name {	\
		public:	\
			name() { }	\
			static inline OutType(T1,T2) apply(T1 a, T2 b)	\
			{ return a op b; }	\
	};

MakeDoubleOp(ApAdd, +)
MakeDoubleOp(ApSubtract, -)
MakeDoubleOp(ApMultiply, *)
MakeDoubleOp(ApDivide, /)


template<class T1,class T2>
class ApPow {
	public:
		ApPow() { }
		typedef FloatType(OutType(T1,T2)) float_N;
		static inline  float_N apply(T1 a, T2 b)
		{ return pow(float_N(a),float_N(b)); }
};

//****************************************************************************
// BOOLIAN Double operators i.e. a==b, a!=b, etc applicative templates --
//****************************************************************************

#define MakeDoubleOpBool(name, op) \
	template<class T1,class T2> \
	class name {	\
		public:	\
			name() { }	\
			static inline bool apply(T1 a, T2 b)	\
			{ return a op b; }	\
	};

MakeDoubleOpBool(ApEqual, ==)
MakeDoubleOpBool(ApNotEqual, !=)
MakeDoubleOpBool(ApGreater, >)
MakeDoubleOpBool(ApGreaterEqual, >=)
MakeDoubleOpBool(ApLess, <)
MakeDoubleOpBool(ApLessEqual, <=)
MakeDoubleOpBool(ApAnd, &&)
MakeDoubleOpBool(ApOr, ||)





//****************************************************************************
//		the 'operator-2' operations (i.e. +=, *=, etc)
//****************************************************************************

class op2_base { };

#define MakeOp2(name,op) \
	template<class X, class Y> \
	class name : public op2_base { \
	  public:   \
	  static inline void apply(X& x, Y y)  { x op y; }  \
	};

template<class X, class Y>
class ApAssign : public op2_base {
  public:
    static inline void apply(X& x, Y y)
    { x = (X)y; }
};

MakeOp2(ApAdd2, +=)
MakeOp2(ApSubtract2, -=)
MakeOp2(ApDivide2, /=)
MakeOp2(ApMultiply2, *=)



/***************************************************************************/
// Single operators i.e. abs(b), cos(b), etc applicative templates --
//****************************************************************************

#define MakeSingleOp(name, op, convert) \
	template<class T1> \
	class name {	\
		public:	\
			typedef convert numtype; \
			name() { }	\
			static inline convert apply(T1 a)	\
			{ return op(a); }	\
	};

#define MakeSingleOpKeep(name, op) MakeSingleOp(name, op, T1)

 MakeSingleOpKeep(ApAbs, abs)
MakeSingleOpKeep(ApCos, cos)
MakeSingleOpKeep(ApSin, sin)
MakeSingleOpKeep(ApTan, tan)
MakeSingleOpKeep(ApExp, exp)
MakeSingleOpKeep(ApAcos, acos)
MakeSingleOpKeep(ApAsin, asin)
MakeSingleOpKeep(ApAtan, atan)
MakeSingleOpKeep(ApLog, log)
MakeSingleOpKeep(ApLog10, log10)
MakeSingleOpKeep(ApCosh, cosh)
MakeSingleOpKeep(ApSinh, sinh)
MakeSingleOpKeep(ApTanh, tanh)
MakeSingleOpKeep(ApAcosh, acosh)
MakeSingleOpKeep(ApAsinh, asinh)
MakeSingleOpKeep(ApAtanh, atanh)
MakeSingleOpKeep(ApCiel, ceil)
MakeSingleOpKeep(ApFloor, floor)
#ifdef HAVE_FINITE
//defined below
#elif HAVE_ISNAN
MakeSingleOpKeep(ApNan, isnan)
#endif
//MakeSingleOpKeep(ApSqrt, std::sqrt)
MakeSingleOpKeep(ApNeg, -)
MakeSingleOp(ApRe, Re, double)
MakeSingleOp(ApIm, Im, double)
MakeSingleOpKeep(ApConj, conj)

template<class T1>
class ApSqrt {
	public:
		typedef T1 numtype;
		ApSqrt() { }
		static inline T1 apply(T1 a)
		{ return std::sqrt(a); }
};


//special unary operators...not defined in the 'math.h'
template<class T1>
class ApNan {
public:
	typedef T1 numtype;
	ApNan() { }
	static inline T1 apply(T1 a)
	{ return !finite(a); }
};


template<class T1>
class ApSqr {
	public:
		typedef T1 numtype;
		ApSqr() { }
		static inline T1 apply(T1 a)
		{ return a*a; }
};

template<class T1>
class ApCube {
	public:
		typedef T1 numtype;
		ApCube() { }
		static inline T1 apply(T1 a)
		{ return a*a*a; }
};

template<class T1>
class ApForthPow {
	public:
		typedef T1 numtype;
		ApForthPow() { }
		static inline T1 apply(T1 a)
		{ return a*a*a*a; }
};

template<class T1>
class ApFifthPow {
	public:
		typedef T1 numtype;
		ApFifthPow() { }
		static inline T1 apply(T1 a)
		{ return a*a*a*a*a; }
};



//special unary operators...not defined in the 'math.h'
template<class T1>
class ApChop {
	private:
		//static const double eps=1.e-12;
	public:
		typedef T1 numtype;
		ApChop() { }
		static inline T1 apply(T1 a)
		{ return abs(a)>1.e-12?a:0; }
};

/*inline double chop(const double &a, double lim)
{
	return abs(a)>lim?a:0;
}

inline float chop(const float &a, double lim)
{
	return abs(a)>lim?a:0;
}

inline int chop(const int &a, double lim)
{
	return abs(a)>lim?a:0;
}
*/
END_BL_NAMESPACE


#endif

