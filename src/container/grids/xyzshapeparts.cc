 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-11-01
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
 	xyzhapeparts.cc-->a collection shape function classes used to
 	generate the shape grids for the "XYZshape" class (these
 	act as the "ShapeEng_t"
 */

#ifndef _XYZshapeparts_cc_
#define _XYZshapeparts_cc_

#include "container/grids/coords.h"
#include "container/Vector/Vector.h"


//base class file
#include "container/grids/xyzgenshape.h"
#include "container/grids/xyzshapeparts.h"


BEGIN_BL_NAMESPACE

//---------FULL shape-------------------------------------

XYZfull_ &XYZfull_::operator=(const XYZfull_ &rhs)
{
	if(this==&rhs) return *this;
	GeneralXYZShape::operator=(rhs);
	return *this;
}

//-----------End Full-----------------------------------

//-----------------------------------------------------------

//-----------Rectangular solid-----------------------------

XYZrect_ &XYZrect_::operator=(const XYZrect_ &rhs)
{
	if(this==&rhs) return *this;
	GeneralXYZShape::operator=(rhs);
	return *this;
}

bool XYZrect_::ShapeFunc(double x, double y, double z)
{

	if(x>=min().x() && x<=max().x())
	{
		if(y>=min().y() && y<=max().y())
		{
			if(z>=min().z() && z<=max().z())
			{
				return true;
			}
		}
	}
	return false;
}

bool XYZrect_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}

//----------End Rect-----------------------


//-----------Cylinder-----------------------------

XYZcylinder_ &XYZcylinder_::operator=(const XYZcylinder_ &rhs)
{
	if(this==&rhs) return *this;
	GeneralXYZShape::operator=(rhs);
	return *this;
}

bool XYZcylinder_::ShapeFunc(double x, double y, double z)
{
	double tmr=sqrt(x*x+y*y);
	double tmph=0.0;
	if(tmr!=0) tmph=acos(y/tmr);
	if(x<0.0){ tmph=PI2-tmph;	}
	if(tmr>=min().x() && tmr<=max().x())
	{
		if(tmph>=min().y() && tmph<=max().y())
		{
			if(z>=min().z() && z<=max().z())
			{
				//cout<<tmph<<endl;
				return true;
			}
		}
	}
	return false;
}

bool XYZcylinder_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}

//----------End Cylinder-----------------------

//--------------Begin PLANAR Types---------

XYZBasicPlane &XYZBasicPlane::operator=(const XYZBasicPlane &rhs)
{
	if(&rhs==this) return *this;
	m_=rhs.m_;
	b_=rhs.b_;
	sc_=rhs.sc_;
	above_=rhs.above_;
	return *this;
}

//------------XY-test-Plane-----------
XYZplaneXY_ &XYZplaneXY_::operator=(const XYZplaneXY_ &rhs)
{
	if(this==&rhs) return *this;
	XYZBasicPlane::operator=(rhs);
	return *this;
}

bool XYZplaneXY_::ShapeFunc(double x, double y, double z)
{
	if(above()){ if(y>(x*m()+b()))	return true;	}
	else{ if(y<(x*m()+b()))	return true;}
	return false;
}

bool XYZplaneXY_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}


//---------XZ-test-Plane---------------
XYZplaneXZ_ &XYZplaneXZ_::operator=(const XYZplaneXZ_ &rhs)
{
	if(this==&rhs) return *this;
	XYZBasicPlane::operator=(rhs);
	return *this;
}

bool XYZplaneXZ_::ShapeFunc(double x, double y, double z)
{
	if(above()){ if(z>(x*m()+b()))	return true;	}
	else{ if(z<(x*m()+b()))	return true;}
	return false;
}

bool XYZplaneXZ_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}

//---------YZ-test-Plane---------------
XYZplaneYZ_ &XYZplaneYZ_::operator=(const XYZplaneYZ_ &rhs)
{
	if(this==&rhs) return *this;
	XYZBasicPlane::operator=(rhs);
	return *this;
}

bool XYZplaneYZ_::ShapeFunc(double x, double y, double z)
{
	if(above()){ if(z>(y*m()+b()))	return true;	}
	else{ if(z<(y*m()+b()))	return true;}
	return false;
}

bool XYZplaneYZ_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}


//-------------- Basic 3D plane-------------
XYZBasic3Dplane &XYZBasic3Dplane::operator=(const XYZBasic3Dplane &rhs)
{
	if(&rhs==this) return *this;
	GeneralXYZShape::operator=(rhs);
	return *this;
}


XYZBasic3Dplane &XYZBasic3Dplane::operator=(const XYZBasicPlane &rhs)
{
	min()=rhs.m();
	max()=rhs.b();
	center()=rhs.Scale();
	return *this;
}

//--------------3D plane slicer---------------
XYZ3Dplane_ &XYZ3Dplane_::operator=(const XYZ3Dplane_ &rhs)
{
	if(this==&rhs) return *this;
	XYZBasic3Dplane::operator=(rhs);
	return *this;
}

bool XYZ3Dplane_::ShapeFunc(double x, double y, double z)
{
	if(z>(z*m().z()+b().z()) && above()){
		if(y>y*m().y()+b().y() && above()){
			if(x>x*m().x()+b().x() && above()){
				return true;
			}
		}
	}
	return false;
}

bool XYZ3Dplane_::ShapeFunc(coord<> &xyz)
{
	return ShapeFunc(xyz.x(), xyz.y(), xyz.z());
}

END_BL_NAMESPACE


#endif


