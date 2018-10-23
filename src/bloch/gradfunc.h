/* gradfunc.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-21-01
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
 	gradfunc.h-->simply a set of classes to simplify the
 	application of variational parameters over spacial coordinates

 	things like temperature gradients, Bo inhomogeniety, etc

 	the acctual implimentation of the gradients are performed in other classes


 	the only nessesary piece is the

 		void operator()('thing to alter', 'grid posistion')

 		the 'thing to alter' can be any class you wish
 		that you can change or alter...typically
 		it is a scaler (i.e. a double) but it could be
 		a matrix, coord<>, vector, etc...that choice i leave to you

 		the 'grid position' is a coord<> (x,y,z)...

 	the gradfunc.cc contains any bits that are compilable...

 */


#ifndef _grad_funcs_h_
#define _grad_funcs_h_ 1


#include "container/grids/coords.h"
#include "utils/random.h" 		//for Random()

BEGIN_BL_NAMESPACE



//some globals for the output of things
extern const std::string gf_name;

// this is the 'base' class in essence it is just a way of
// grouping all the 'gradfuncs' together..this one does NOTHING to the input

class NullGradFunc
{

	public:
		NullGradFunc()
		{}

		NullGradFunc(NullGradFunc &cp)
		{}

		// this allows you to set the NullGrad from anything really...
		// as this does nothing
		template<class AnyScaleFunc>
		NullGradFunc &operator=(AnyScaleFunc &rhs)
		{
			return *this;
		}

		template<class Alter_t>
		inline void operator()(Alter_t &in, coord<> pos)
		{
			return;
		}
};

std::ostream &operator<<(std::ostream &oo, NullGradFunc &out);


//a Scaling gradient...sort of 'useless' but SCALES everything by some number
// Alter_t=SCALE*Alter_t

class ScaleGradFunc : public NullGradFunc
{
	private:
		double sc_;

	public:
		ScaleGradFunc():
			NullGradFunc(), sc_(0.)
		{}

		ScaleGradFunc(double sc):
			NullGradFunc(), sc_(sc)
		{}

		ScaleGradFunc(ScaleGradFunc &cp):
			NullGradFunc(), sc_(cp.sc_)
		{}

		inline double Scale(){ return sc_;	}

		template<class num_t>
		void SetScale(num_t in){	sc_=double(in);	}

		ScaleGradFunc &operator=(ScaleGradFunc &rhs);

		template<class Alter_t>
		inline void operator()(Alter_t & in, coord<> &pos)
		{
			in*=sc_;
		}

};

std::ostream &operator<<(std::ostream &oo, ScaleGradFunc &out);



/* Linear gradient types...there are 3 directional ones (X,Y,Z) and
   the 3-D linear one, where all three dirs are included

	the basic linear grad scales simply according to
	its present posistion...the amount of scaling
	is represented by a 'global scaling factor' and
	a line equation y=mx+b

	for X-dir
		Alter*=(m*pos.x()+b)*scale
	for Y-dir
		Alter*=(m*pos.y()+b)*scale
	for Z-dir
		Alter*=(m*pos.z()+b)*scale

	for all-dir
		Alter*=(mx*pos.x()+bx)*(my*pos.y()+by)*(mz*pos.z()+bz)*scale


*/


class BasicLinearGradFunc : public NullGradFunc
{
	private:
		double m_;
		double b_;
		double sc_;

	public:

		BasicLinearGradFunc():
			NullGradFunc(),
			m_(1), b_(0), sc_(1)
		{}

		BasicLinearGradFunc(double mx, double bx):
			NullGradFunc(),
			m_(mx), b_(bx), sc_(1)
		{}

		BasicLinearGradFunc(double mx, double bx, double sc):
			NullGradFunc(),
			m_(mx), b_(bx), sc_(sc)
		{}

		BasicLinearGradFunc(BasicLinearGradFunc &cp):
			NullGradFunc(),
			m_(cp.m_), b_(cp.b_), sc_(cp.sc_)
		{}

		inline double &m(){	return m_;	}
		inline double &b(){	return b_;	}
		inline double &Scale(){	return sc_;	}

		inline double m() const {	return m_;	}
		inline double b() const {	return b_;	}
		inline double Scale() const {	return sc_;	}

		inline void Setm(double mm){	m_=mm;	}
		inline void Setb(double bb){	b_=bb;	}
		inline void SetScale(double sc){	sc_=sc;	}

		BasicLinearGradFunc &operator=(BasicLinearGradFunc &rhs);

};


class LinearXGradFunc : public BasicLinearGradFunc
{

	public:
		static const char id='x';

		LinearXGradFunc():
			BasicLinearGradFunc()
		{}

		LinearXGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearXGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearXGradFunc(LinearXGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs*=Scale()*(pos.x()*m()+b());
		}

};

std::ostream &operator<<(std::ostream &oo,const LinearXGradFunc &out);

class LinearYGradFunc : public BasicLinearGradFunc
{

	public:
		static const char id='y';

		LinearYGradFunc():
			BasicLinearGradFunc()
		{}

		LinearYGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearYGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearYGradFunc(LinearYGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs*=Scale()*(pos.y()*m()+b());
		}

};


std::ostream &operator<<(std::ostream &oo,const LinearYGradFunc &out);


// Linear Z Gradients
class LinearZGradFunc : public BasicLinearGradFunc
{

	public:
		static const char id='z';

		LinearZGradFunc():
			BasicLinearGradFunc()
		{}

		LinearZGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearZGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearZGradFunc(LinearZGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs*=Scale()*(pos.z()*m()+b());
		}

};

std::ostream &operator<<(std::ostream &oo,const LinearZGradFunc &out);



// 3D linear grad types
class Basic3DLinearGradFunc : public NullGradFunc
{
	private:
		coord<> m_;
		coord<> b_;
		coord<> sc_;

	public:

		Basic3DLinearGradFunc():
			NullGradFunc(),
			m_(1), b_(0), sc_(1)
		{}

		Basic3DLinearGradFunc(coord<> mx, coord<> bx):
			NullGradFunc(),
			m_(mx), b_(bx), sc_(1)
		{}

		Basic3DLinearGradFunc(coord<> mx, coord<> bx, coord<> sc):
			NullGradFunc(),
			m_(mx), b_(bx), sc_(sc)
		{}

		Basic3DLinearGradFunc(Basic3DLinearGradFunc &cp):
			NullGradFunc(),
			m_(cp.m_), b_(cp.b_), sc_(cp.sc_)
		{}

		inline coord<> &m(){	return m_;	}
		inline coord<> &b(){	return b_;	}
		inline coord<> &Scale(){	return sc_;	}

		inline coord<> m() const {	return m_;	}
		inline coord<> b() const {	return b_;	}
		inline coord<> Scale() const {	return sc_;	}

		inline void Setm(coord<> mm){	m_=mm;	}
		inline void Setb(coord<> bb){	b_=bb;	}
		inline void SetScale(coord<> sc){	sc_=sc;	}

		Basic3DLinearGradFunc &operator=(Basic3DLinearGradFunc &rhs);

		Basic3DLinearGradFunc &operator=(BasicLinearGradFunc &rhs);

};


class Linear3DGradFunc : public Basic3DLinearGradFunc
{
	public:
		static const char id='a';

		Linear3DGradFunc():
			Basic3DLinearGradFunc()
		{}

		Linear3DGradFunc(coord<> m, coord<> b):
			Basic3DLinearGradFunc(m,b)
		{}

		Linear3DGradFunc(coord<> m, coord<> b, coord<> sc):
			Basic3DLinearGradFunc(m,b,sc)
		{}


		// the mighty driver
		// there is one for a coord<> (where each direction is ltered as it should)
		// and another for scalers

		//the scaler func
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs*=Scale().z()*(pos.z()*m().z()+b().z())*
			     Scale().x()*(pos.x()*m().x()+b().x())*
			     Scale().y()*(pos.y()*m().y()+b().y());
		}

		inline void operator()(coord<> &rhs, coord<> &pos)
		{
			rhs*=Scale()*(pos*m()+b());
		}
};

std::ostream &operator<<(std::ostream &oo,const Linear3DGradFunc &out);

/*****************************************************/

//Linear OFFSET grad funcs...
// these ADD bits to the input (rather then 'scale' them like the previous ones)
// these are probably the most usefull as they allow for the setting
// of some Delta around a constant value
// same 4 basic types, 3 directional ones, and one 3D one

class LinearOffsetXGradFunc : public BasicLinearGradFunc
{

	public:
		static const char id='x';

		LinearOffsetXGradFunc():
			BasicLinearGradFunc()
		{}

		LinearOffsetXGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearOffsetXGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearOffsetXGradFunc(LinearOffsetXGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.x()*m()+b());
		}

};

std::ostream &operator<<(std::ostream &oo,const LinearOffsetXGradFunc &out);

class LinearOffsetYGradFunc : public BasicLinearGradFunc
{

	public:

		static const char id='y';

		LinearOffsetYGradFunc():
			BasicLinearGradFunc()
		{}

		LinearOffsetYGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearOffsetYGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearOffsetYGradFunc(LinearOffsetYGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.y()*m()+b());
		}

};


std::ostream &operator<<(std::ostream &oo,const LinearOffsetYGradFunc &out);


// Linear Z Gradients
class LinearOffsetZGradFunc : public BasicLinearGradFunc
{

	public:
		static const char id='z';

		LinearOffsetZGradFunc():
			BasicLinearGradFunc()
		{}

		LinearOffsetZGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx,bx)
		{}

		LinearOffsetZGradFunc(double mx, double bx, double sc):
			BasicLinearGradFunc(mx,bx,sc)
		{}

		LinearOffsetZGradFunc(LinearOffsetZGradFunc &cp):
			BasicLinearGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.z()*m()+b());
		}

};


//3D linear offset grad func
class LinearOffset3DGradFunc : public Basic3DLinearGradFunc
{
	public:
		static const char id='a';

		LinearOffset3DGradFunc():
			Basic3DLinearGradFunc()
		{}

		LinearOffset3DGradFunc(coord<> m, coord<> b):
			Basic3DLinearGradFunc(m,b)
		{}

		LinearOffset3DGradFunc(coord<> m, coord<> b, coord<> sc):
			Basic3DLinearGradFunc(m,b,sc)
		{}


		// the mighty driver
		// there is one for a coord<> (where each direction is ltered as it should)
		// and another for scalers

		//the scaler func
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale().z()*(pos.z()*m().z()+b().z())*
				Scale().x()*(pos.x()*m().x()+b().x())*
				Scale().y()*(pos.y()*m().y()+b().y());
		}

		inline void operator()(coord<> &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos*m()+b());
		}
};

std::ostream &operator<<(std::ostream &oo,const LinearOffset3DGradFunc &out);


/*********************************************************/

//RANDOM gradients funcs...
// okay, okay, okay...not exactly 'gradient functions', but
// they do perform the same task...altering the input...but not
// based on the posistion..

// the random number generator can be found in 'utils.h'
// it spits out doubles from 0..1
// however to make this generator a bit more functional it would be nice
// to generate some number from (min..max)


//the basic Random Grad Function
// performs A SCALING
//  a*=Random();

template<class RandEngine_t=UniformRandom<double> >
class RandomGradFunc : public RandEngine_t, NullGradFunc
{

	public:
		static const char id='a';

		typedef typename RandEngine_t::numtype numtype;
		RandomGradFunc():
			RandEngine_t()
		{}

		RandomGradFunc(numtype mean):
			RandEngine_t(mean)
		{}

		RandomGradFunc(numtype low, numtype high):
			RandEngine_t(low, high)
		{}

		RandomGradFunc(numtype low, numtype high, numtype mean):
			RandEngine_t(low, high, mean)
		{}

		RandomGradFunc(RandomGradFunc &cp):
			RandEngine_t(cp)
		{}

		RandEngine_t Engine(){	return *this;	}

		template<class Alter_t>
		void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs*=RandEngine_t::Random();
		}

};

template<class RandEngine_t>
std::ostream &operator<<(std::ostream &oo,const RandomGradFunc<RandEngine_t> &out)
{
	oo<<gf_name<<" Random Scale Gradient: a*=(";
	RandEngine_t tmp(out);
	oo<<tmp;
	oo<<")";
	return oo;
}
//the basic Random Grad Function
// performs A Delta OFFSET
//  a+=Random();

template<class RandEngine_t=UniformRandom<double> >
class RandomOffsetGradFunc : public RandEngine_t, NullGradFunc
{

	public:
		static const char id='a';

		typedef typename RandEngine_t::numtype numtype;
		RandomOffsetGradFunc():
			RandEngine_t()
		{}

		RandomOffsetGradFunc(numtype mean):
			RandEngine_t(mean)
		{}

		RandomOffsetGradFunc(numtype low, numtype high):
			RandEngine_t(low, high)
		{}

		RandomOffsetGradFunc(numtype low, numtype high, numtype mean):
			RandEngine_t(low, high, mean)
		{}

		RandomOffsetGradFunc(RandomOffsetGradFunc &cp):
			RandEngine_t(cp)
		{}


		template<class Alter_t>
		void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=RandEngine_t::Random();
		}

};

template<class RandEngine_t>
std::ostream &operator<<(std::ostream &oo,const RandomOffsetGradFunc<RandEngine_t> &out)
{
	oo<<gf_name<<" Random Offset Gradient: a+=(";
	RandEngine_t tmp(out);
	oo<<tmp;
	oo<<")";
	return oo;
}


/********************************************/
// Radial Functions.....
//
// Base class...takes the 'BasicLinear' calss and adds the 'offset' term
// so one can recenter the 0 point for 'R'

class BasicRadialGradFunc : public BasicLinearGradFunc
{
	private:
		double offset_;
	public:

		BasicRadialGradFunc():
			BasicLinearGradFunc(),
			offset_(0)
		{}

		BasicRadialGradFunc(double mx, double bx):
			BasicLinearGradFunc(mx, bx),
			offset_(0)
		{}

		BasicRadialGradFunc(double mx, double bx, double off):
			BasicLinearGradFunc(mx, bx),
			offset_(off)
		{}

		BasicRadialGradFunc(double mx, double bx, double off, double sc):
			BasicLinearGradFunc(mx, bx,sc),
			offset_(off)
		{}

		BasicRadialGradFunc(BasicRadialGradFunc &cp):
			BasicLinearGradFunc(cp),
			offset_(cp.offset_)
		{}

		inline double &offset(){	return offset_;	}

		inline double offset() const {	return offset_;	}

		inline void SetOffset(double sc){	offset_=sc;	}

		BasicRadialGradFunc &operator=(BasicRadialGradFunc &rhs);
		BasicRadialGradFunc &operator=(BasicLinearGradFunc &rhs);
};



// SPHERICAL radial LINEAR offset function
// a+=Scale*m*(sqrt(x^2+y^2+z^2)+offset)+b

class SphericalLinearOffsetGradFunc : public BasicRadialGradFunc
{

	public:
		static const char id='a';

		SphericalLinearOffsetGradFunc():
			BasicRadialGradFunc()
		{}

		SphericalLinearOffsetGradFunc(double mx, double bx):
			BasicRadialGradFunc(mx,bx)
		{}

		SphericalLinearOffsetGradFunc(double mx, double bx, double off):
			BasicRadialGradFunc(mx,bx,off)
		{}

		SphericalLinearOffsetGradFunc(double mx, double bx, double off, double sc):
			BasicRadialGradFunc(mx,bx,off,sc)
		{}

		SphericalLinearOffsetGradFunc(SphericalLinearOffsetGradFunc &cp):
			BasicRadialGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*((norm(pos)+offset())*m()+b());
		}

};

std::ostream &operator<<(std::ostream &oo,const SphericalLinearOffsetGradFunc &out);


// Cylindrical radial LINEAR offset function
// a+=Scale*(m*(sqrt(x^2+y^2)+offset)+b)

class CylindricalLinearOffsetGradFunc : public BasicRadialGradFunc
{

	public:
		static const char id='a';

		CylindricalLinearOffsetGradFunc():
			BasicRadialGradFunc()
		{}

		CylindricalLinearOffsetGradFunc(double mx, double bx):
			BasicRadialGradFunc(mx,bx)
		{}

		CylindricalLinearOffsetGradFunc(double mx, double bx, double off):
			BasicRadialGradFunc(mx,bx,off)
		{}

		CylindricalLinearOffsetGradFunc(double mx, double bx, double off, double sc):
			BasicRadialGradFunc(mx,bx,off,sc)
		{}

		CylindricalLinearOffsetGradFunc(CylindricalLinearOffsetGradFunc &cp):
			BasicRadialGradFunc(cp)
		{}


		// the mighty driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*((std::sqrt(pos.x()*pos.x()+pos.y()*pos.y())+offset())*m()+b());
		}

};

std::ostream &operator<<(std::ostream &oo,const CylindricalLinearOffsetGradFunc &out);


END_BL_NAMESPACE


#endif

