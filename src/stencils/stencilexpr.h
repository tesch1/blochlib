/* stencilexpr.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-30-01
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

 	stencilexpr.h-> expression templates for the stencil operators

 	allows multiple stencils to be added/subtracted/mult/div, etc
 	without creation of temporaries

 */


#ifndef _stencil_expr_h_
#define _stencil_expr_h_ 1


#include "container/rankType.h"


BEGIN_BL_NAMESPACE


template<class S1, class S2>
class StencilAdd
{
	public:
		StencilAdd() { }

		typedef OutType(typename S1::numtype, typename S2::numtype) numtype;

		template<class container>
		inline numtype apply(S1 &a, S2 &b, container &in)
		{ return a.apply(in) + b.apply(in); }

		template<class container>
		inline numtype operator()(S1 &a, S2 &b, container &in)
		{ return a.apply(in) + b.apply(in); }

};

template<class S1, class S2>
class StencilSub
{
	public:
		StencilSub() { }

		typedef OutType(typename S1::numtype, typename S2::numtype) numtype;

		template<class container>
		inline numtype apply(S1 &a, S2 &b, container &in)
		{ return a.apply(in) - b.apply(in); }

		template<class container>
		inline numtype operator()(S1 &a, S2 &b, container &in)
		{ return a.apply(in) - b.apply(in); }

};


template<class S1, class S2>
class StencilMul
{
	public:
		StencilMul() { }

		typedef OutType(typename S1::numtype, typename S2::numtype) numtype;

		template<class container>
		inline numtype apply(S1 &a, S2 &b, container &in)
		{ return a.apply(in) * b.apply(in); }

		template<class container>
		inline numtype operator()(S1 &a, S2 &b, container &in)
		{ return a.apply(in) * b.apply(in); }

};

template<class S1, class S2>
class StencilDiv
{
	public:
		StencilDiv() { }

		typedef OutType(typename S1::numtype, typename S2::numtype) numtype;

		template<class container>
		inline numtype apply(S1 &a, S2 &b, container &in)
		{ return a.apply(in) / b.apply(in); }

		template<class container>
		inline numtype operator()(S1 &a, S2 &b, container &in)
		{ return a.apply(in) / b.apply(in); }

};

//there lovely operator overloads give the user the applity to
// do "ShapeFunc1 || ShapeFunc2" to create the XYZor<shapefunc1, shapefunc2> set up
// so you can create these large strucutres with normal math notation (a bit easier to read and understand)


template<class expr>
class StencilExpr
{
	private:
		expr iter_;
		StencilExpr() {}

	public:

		typedef typename expr::numtype numtype;

		StencilExpr(const expr& a)  : iter_(a) { }

		template<class container>
		numtype ShapeFunc(container &in)
		{	return iter_.apply(in);	}

		template<class container>
		numtype operator()(container &in)
		{	return iter_.apply(in);	}

};

template<class NumType>
class StencilExprConst
{
	private:
		NumType iter_;
		StencilExprConst() {}

	public:

		typedef NumType numtype;

		StencilExprConst(const NumType& a)  : iter_(a) { }

		template<class container>
		numtype ShapeFunc(container &in)
		{	return iter_;	}

		template<class container>
		numtype operator()(container &in)
		{	return iter_;	}

};


template<class S1, class S2, class Op>
class StencilBinExprOp {

private:
    S1 iter1_;
    S2 iter2_;
    StencilBinExprOp(){}

public:
	typedef OutType(typename S1::numtype, typename S2::numtype) numtype;

    StencilBinExprOp(const S1& a, const S2& b) : iter1_(a), iter2_(b)  { }

	template<class container>
	inline numtype apply(container &in)
	{ return Op::apply(iter1_, iter2_, in);	}

	template<class container>
	inline numtype operator()(container &in)
	{ return Op::apply(iter1_, iter2_, in); }

};

/*

template<class SF1, class SF2, class NumType_t1, class NumType_t2>
StencilExpr<StencilBinExprOp<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2>,
 StencilAdd<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > >
 operator+ (const StencilEng<SF1, NumType_t1> &lhs, const StencilEng<SF2, NumType_t2> &rhs)
{
	typedef StencilBinExprOp<StencilEng<SF1,NumType_t2>, StencilEng<SF2,NumType_t2>,
	        StencilAdd<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > expr;
	return StencilExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2, class NumType_t1, class NumType_t2>
StencilExpr<StencilBinExprOp<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2>,
 StencilSub<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > >
 operator- (const StencilEng<SF1, NumType_t1> &lhs, const StencilEng<SF2, NumType_t2> &rhs)
{
	typedef StencilBinExprOp<StencilEng<SF1,NumType_t2>, StencilEng<SF2,NumType_t2>,
	        StencilSub<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > expr;
	return StencilExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2, class NumType_t1, class NumType_t2>
StencilExpr<StencilBinExprOp<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2>,
 StencilDiv<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > >
 operator/ (const StencilEng<SF1, NumType_t1> &lhs, const StencilEng<SF2, NumType_t2> &rhs)
{
	typedef StencilBinExprOp<StencilEng<SF1,NumType_t2>, StencilEng<SF2,NumType_t2>,
	        StencilDiv<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > expr;
	return StencilExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2, class NumType_t1, class NumType_t2>
StencilExpr<StencilBinExprOp<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2>,
 StencilMul<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > >
 operator* (const StencilEng<SF1, NumType_t1> &lhs, const StencilEng<SF2, NumType_t2> &rhs)
{
	typedef StencilBinExprOp<StencilEng<SF1,NumType_t2>, StencilEng<SF2,NumType_t2>,
	        StencilMul<StencilEng<SF1,NumType_t1>, StencilEng<SF2,NumType_t2> > > expr;
	return StencilExpr<expr >(expr(lhs, rhs));
}

*/

END_BL_NAMESPACE


#endif



