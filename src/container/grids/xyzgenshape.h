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
 	xyzgenshape.h-->the master base class for the 'shape functions' for the class "XYZshape"
 	as of now this does nothing, but it may one day, so i put it here just in case
 	also for those graphical class browsers, this will group all the XYZshapefunc under the
 	same group
 */

#ifndef _GeneralXYZShape_h_
#define _GeneralXYZShape_h_


BEGIN_BL_NAMESPACE


/* hold the base class for the SHAPE parts of all the shapes WITHIN an XYZ grid
   this class acts as a functional class..it does not contain the data for the
   shape, only the function to generate it (the base class, this, does not contain the function
   it simply holds the data and functions global to all XYZtypes...

   for instance the XZYcylinder class would have this as a base class
   and two importent function

   bool operator()(double x, double y, double z)
     and
   bool operator()(coord<> pos)

   these operators take in a single xyz point return true if the point is within the
   functional shape, and false if not....

 */
#include "container/grids/coords.h"
#include "container/Vector/Vector.h"


#ifndef PI2
#define PI2 6.28318530717958647692
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

//the master shape container...Allows the use of expression templates
//specifically to the shapes (and not get confused by other classes)

template<class XYZshapes>
class XYZparts :
	public XYZshapes
{
	public:
		XYZparts():
			XYZshapes()
		{}
		XYZparts(const XYZparts &cp):
			XYZshapes(cp)
		{}

		XYZparts(const coord<> &minin,const coord<> &maxin):
			XYZshapes(minin, maxin, coord<>(0))
		{}

		XYZparts(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
			XYZshapes(xmin,xmax,ymin,ymax,zmin,zmax)
		{}

		XYZparts(double mx, double bx):
			XYZshapes(mx,bx)
		{}

		XYZparts(double mx, double bx, double sc):
			XYZshapes(mx,bx,sc)
		{}


		XYZparts(coord<> mx, coord<> bx, coord<> sc):
			XYZshapes(mx,bx,sc)
		{}
};


//The 'master base' class for all shapes....
 class GeneralXYZShape {
 	private:
 		coord<> min_;  // holds{ x,y,z} minima
		coord<> max_;  // holds {x,y,z} maxima
		coord<> center_;	//holds {x,y,z} CENTER::NOTE X,Y,Z center point....


	public:

		static const int Shape=1;

		GeneralXYZShape():
			min_(0), max_(0), center_(0)
		{}

		GeneralXYZShape(const GeneralXYZShape &copy):
			min_(copy.min_), max_(copy.max_), center_(copy.center_)
		{}

		GeneralXYZShape(const coord<> &minin,const  coord<> &maxin):
			min_(minin), max_(maxin), center_(0)
		{}

		GeneralXYZShape(const coord<> &minin,const coord<> &maxin,const coord<> &cen):
			min_(minin), max_(maxin), center_(cen)
		{}

		GeneralXYZShape(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
			min_(xmin, ymin, zmin), max_(xmax, ymax, zmax), center_(0)
		{}


		GeneralXYZShape &operator=(const GeneralXYZShape &rhs)
		{
			if(this==&rhs) return *this;
			min_=rhs.min_;
			max_=rhs.max_;
			center_=rhs.center_;
			return *this;
		}

		coord<> &min(){	return min_;	}
		coord<> &max(){	return max_;	}
		coord<> &center(){	return center_;	}

		coord<> min() const {	return min_;	}
		coord<> max() const {	return max_;	}
		coord<> center() const {	return center_;	}

};


//----------------MULTI Shape Container Templates-----
// Here we define 'basic pairs' of functions
//  The pairs are designed to take in a single coord<>
//  and perform a "AND" or an "OR"
// i.e. OR--> if the point is in either shape it is accepted
///     AND--> if the point must be in BOTH shapes to be accepted
// this basic container is then used by XYZshape as it driver
// the nice thing about this little set up is that pairs can be grown
//
// like so...
//     XYZor<XYZplaneXZ, XYZand<XYZcyl, XYZrect> >
// this little expression means the point can be either in the plane
// OR in both the cyl and rect
//
//  this makes creating complex shapes relatively easy from the basic blocks

template<class ShapeFunc1, class ShapeFunc2>
class XYZand
{
	private:
		ShapeFunc1 shape1;
		ShapeFunc2 shape2;

	public:
		XYZand():
			shape1(), shape2()
		{}

		XYZand(const ShapeFunc1 &in1, const ShapeFunc2 &in2):
			shape1(in1), shape2(in2)
		{}

	//the shape functions
		bool ShapeFunc(double x, double y, double z)
		{
			return shape1.ShapeFunc(x,y,z) && shape2.ShapeFunc(x,y,z);
		}

		bool ShapeFunc(coord<> &xyz)
		{
			return shape1.ShapeFunc(xyz) && shape2.ShapeFunc(xyz);
		}
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}


	//the altering func
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			shape1(rhs, pos);
			shape2(rhs, pos);
		}

};



template<class ShapeFunc1, class ShapeFunc2>
class XYZor
{
	private:
		ShapeFunc1 shape1;
		ShapeFunc2 shape2;

	public:
		XYZor():
			shape1(), shape2()
		{}

		XYZor(const ShapeFunc1 &in1, const ShapeFunc2 &in2):
			shape1(in1), shape2(in2)
		{}

	//the shape functions
		bool ShapeFunc(double x, double y, double z)
		{
			return shape1.ShapeFunc(x,y,z) || shape2.ShapeFunc(x,y,z);
		}

		bool ShapeFunc(coord<> &xyz)
		{
			return shape1.ShapeFunc(xyz) || shape2.ShapeFunc(xyz);
		}
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}


	//the altering func
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			shape1(rhs, pos);
			shape2(rhs, pos);
		}

};


template<class S1, class S2>
class ShapeOr
{
	public:
		ShapeOr() { }
		static inline bool ShapeFunc(S1 &a, S2 &b, double x, double y, double z)
		{ return a.ShapeFunc(x,y,z) || b.ShapeFunc(x,y,z); }

		static inline bool ShapeFunc(S1 &a, S2 &b, coord<> &pos)
		{ return a.ShapeFunc(pos) || b.ShapeFunc(pos); }

};

template<class S1, class S2>
class ShapeAnd
{
	public:
		ShapeAnd() { }
		static inline bool ShapeFunc(S1 &a, S2 &b, double x, double y, double z)
		{ return a.ShapeFunc(x,y,z) && b.ShapeFunc(x,y,z); }

		static inline bool ShapeFunc(S1 &a, S2 &b, coord<> &pos)
		{ return a.ShapeFunc(pos) && b.ShapeFunc(pos); }

};

//there lovely operator overloads give the user the applity to
// do "ShapeFunc1 || ShapeFunc2" to create the XYZor<shapefunc1, shapefunc2> set up
// so you can create these large strucutres with normal math notation (a bit easier to read and understand)


template<class expr>
class ShapeExpr
{
	private:
		expr iter_;

	public:
		ShapeExpr() {}

		ShapeExpr(const expr& a)  : iter_(a) { }

		bool ShapeFunc(double x, double y, double z)
		{	return iter_.ShapeFunc(x,y,z);	}

		bool ShapeFunc(coord<> &pos)
		{	return iter_.ShapeFunc(pos);	}

		bool operator()(coord<> &pos)
		{	return iter_.ShapeFunc(pos);	}

		bool operator()(double x, double y, double z)
		{	return iter_.ShapeFunc(x,y,z);	}
};



template<class S1, class S2, class Op>
class ShapeBinExprOp {

private:
    S1 iter1_;
    S2 iter2_;
    ShapeBinExprOp(){}

public:

    ShapeBinExprOp(const S1& a, const S2& b) : iter1_(a), iter2_(b)  { }

	bool ShapeFunc(double x, double y, double z)
	{	return Op::ShapeFunc(iter1_, iter2_, x,y,z);	}

	bool ShapeFunc(coord<> &pos)
	{	return Op::ShapeFunc(iter1_, iter2_, pos);	}

	bool operator()(coord<> &pos)
	{	return Op::ShapeFunc(iter1_, iter2_, pos);	}

	bool operator()(double x, double y, double z)
	{	return Op::ShapeFunc(iter1_, iter2_, x,y,z);	}

};



template<class SF1, class SF2>
ShapeExpr<ShapeBinExprOp<XYZparts<SF1>, XYZparts<SF2>, ShapeOr<XYZparts<SF1>, XYZparts<SF2> > > >
	operator||(const XYZparts<SF1> &lhs, const XYZparts<SF2> &rhs)
{
	typedef ShapeBinExprOp<XYZparts<SF1>, XYZparts<SF2>, ShapeOr<XYZparts<SF1>, XYZparts<SF2> > > expr;
	return ShapeExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2>
ShapeExpr<ShapeBinExprOp<XYZparts<SF1>, XYZparts<SF2>, ShapeAnd<XYZparts<SF1>, XYZparts<SF2> > > >
	operator&&(const XYZparts<SF1> &lhs, const XYZparts<SF2> &rhs)
{
	typedef ShapeBinExprOp<XYZparts<SF1>, XYZparts<SF2>, ShapeAnd<XYZparts<SF1>, XYZparts<SF2> > > expr;
	return ShapeExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2>
ShapeExpr<ShapeBinExprOp<ShapeExpr<SF1>, XYZparts<SF2>, ShapeOr<ShapeExpr<SF1>, XYZparts<SF2> > > >
	operator||(const ShapeExpr<SF1> &lhs, const XYZparts<SF2> &rhs)
{
	typedef ShapeBinExprOp<ShapeExpr<SF1>, XYZparts<SF2>, ShapeOr<ShapeExpr<SF1>, XYZparts<SF2> > > expr;
	return ShapeExpr<expr >(expr(lhs, rhs));
}

template<class SF1, class SF2>
ShapeExpr<ShapeBinExprOp<XYZparts<SF1>, ShapeExpr<SF2>, ShapeAnd<XYZparts<SF1>, ShapeExpr<SF2> > > >
	operator&&(const XYZparts<SF1> &lhs, const ShapeExpr<SF2> &rhs)
{
	typedef ShapeBinExprOp<XYZparts<SF1>, ShapeExpr<SF2>, ShapeAnd<XYZparts<SF1>, ShapeExpr<SF2> > > expr;
	return ShapeExpr<expr >(expr(lhs, rhs));
}




//here is a little macro that setup a class that uses the
// expressions above so that we can pass the specialized grid shape around...

//The syntax would be something like
// CreateShape(MyNewSHape) {	return XYZrect || XYZcylinder;	}

#define CreateShape(NAME) \
class NAME ## _cl	\
{	\
	public:	\
		\
		inline bool ShapeFunc(coord<> &pos)	\
		{	return ShapeFunc(pos.x(), pos.y(), pos.z());	}	\
	\
		inline bool operator()(coord<> &pos)	\
		{	return ShapeFunc(pos.x(), pos.y(), pos.z());	}	\
	\
		inline bool operator()(double x, double y, double z)	\
		{	return ShapeFunc(x, y, z);	}	\
	\
		inline bool ShapeFunc(double x, double y, double z);	\
};	\
\
ShapeExpr<NAME ## _cl> NAME;	\
\
\
inline bool NAME ## _cl::ShapeFunc(double x, double y, double z) \


END_BL_NAMESPACE


#endif


