/* gradfunc.cc ********/


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
 	gradfunc.cc-->simply a set of classes to simplify the
 	application of variational parameters over spacial coordinates

 	things like temperature gradients, Bo inhomogeniety, etc

 	the acctual implimentation of the gradients are performed in other classes


 	the only nessesary piece is the

 		void operator()('thing to alter', 'grid posistion')

 */


#ifndef _grad_funcs_cc_
#define _grad_funcs_cc_

#include "bloch/gradfunc.h"

BEGIN_BL_NAMESPACE


const std::string gf_name="Gradient Function: ";


/***********************************************/
//NullGradFunc Bits
std::ostream &operator<<(std::ostream &oo, NullGradFunc &out)
{
	oo<<gf_name<<" Null Gradient Function";
	return oo;
}


/************************************************/
//Scaling Grad function bits

ScaleGradFunc &ScaleGradFunc::operator=(ScaleGradFunc &rhs)
{
	if(&rhs==this) return *this;
	sc_=rhs.sc_;
	return *this;
}

std::ostream &operator<<(std::ostream &oo, ScaleGradFunc &out)
{
	oo<<gf_name<<" Scaling Function: a*="<<out.Scale();
	return oo;
}

/**********************************************/
// Linear grad funcs...

BasicLinearGradFunc &BasicLinearGradFunc::operator=(BasicLinearGradFunc &rhs)
{
	if(&rhs==this) return *this;
	m_=rhs.m_;
	b_=rhs.b_;
	sc_=rhs.sc_;
	return *this;
}

Basic3DLinearGradFunc &Basic3DLinearGradFunc::operator=(Basic3DLinearGradFunc &rhs)
{
	if(&rhs==this) return *this;
	m_=rhs.m_;
	b_=rhs.b_;
	sc_=rhs.sc_;
	return *this;
}

Basic3DLinearGradFunc &Basic3DLinearGradFunc::operator=(BasicLinearGradFunc &rhs)
{
	m_=rhs.m();
	b_=rhs.b();
	sc_=rhs.Scale();
	return *this;
}

std::ostream &operator<<(std::ostream &oo,const LinearXGradFunc &out)
{
	oo<<gf_name<<" Linear X Gradient: a*="<<out.Scale()<<"*("<<out.m()<<"*x+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearYGradFunc &out)
{
	oo<<gf_name<<" Linear Y Gradient: a*="<<out.Scale()<<"*("<<out.m()<<"*y+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearZGradFunc &out)
{
	oo<<gf_name<<" Linear Z Gradient: a*="<<out.Scale()<<"*("<<out.m()<<"*z+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const Linear3DGradFunc &out)
{
	oo<<gf_name<<" Linear 3D Gradient: a*=("<<out.Scale()<<")*(("<<out.m()<<")*(x,y,z)+("<<out.b()<<"))";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearOffsetXGradFunc &out)
{
	oo<<gf_name<<" Linear Offset X Gradient: a+="<<out.Scale()<<"*("<<out.m()<<"*x+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearOffsetYGradFunc &out)
{
	oo<<gf_name<<" Linear Offset Y Gradient: a+="<<out.Scale()<<"*("<<out.m()<<"*y+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearOffsetZGradFunc &out)
{
	oo<<gf_name<<" Linear Offset Z Gradient: a+="<<out.Scale()<<"*("<<out.m()<<"*z+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const LinearOffset3DGradFunc &out)
{
	oo<<gf_name<<" Linear Offset 3D Gradient: a+=("<<out.Scale()<<")*(("<<out.m()<<")*(x,y,z)+("<<out.b()<<"))";
	return oo;
}

/*************************************/
//Radial function bits

BasicRadialGradFunc &BasicRadialGradFunc::operator=(BasicRadialGradFunc &rhs)
{
	if(&rhs==this) return *this;
	BasicLinearGradFunc::operator=(rhs);
	offset_=rhs.offset_;
	return *this;
}

BasicRadialGradFunc &BasicRadialGradFunc::operator=(BasicLinearGradFunc &rhs)
{
	BasicLinearGradFunc::operator=(rhs);
	offset_=0.0;
	return *this;
}

std::ostream &operator<<(std::ostream &oo,const SphericalLinearOffsetGradFunc &out)
{
	oo<<gf_name<<" Spherical Offset Gradient: a+="<<out.Scale()<<"*("<<out.m()<<"*(R+"<<out.offset()<<")+"<<out.b()<<")";
	return oo;
}

std::ostream &operator<<(std::ostream &oo,const CylindricalLinearOffsetGradFunc &out)
{
	oo<<gf_name<<" Cylindrical Offset Gradient: a+="<<out.Scale()<<"*("<<out.m()<<"*(R+"<<out.offset()<<")+"<<out.b()<<")";
	return oo;
}


END_BL_NAMESPACE




#endif


