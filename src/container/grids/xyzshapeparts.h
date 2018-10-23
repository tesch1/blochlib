
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
 	xyzhapeparts.h-->a collection shape function classes used to
 	generate the shape grids for the "XYZshape" class (these
 	act as the "ShapeEng_t"
 */

#ifndef _XYZshapeparts_h_
#define _XYZshapeparts_h_

#include "container/grids/coords.h"
#include "container/Vector/Vector.h"
#include "container/rankType.h"


//base class file
#include "container/grids/xyzgenshape.h"


BEGIN_BL_NAMESPACE

//file simply includes all those little header parts for the shape frunctions...

/* NOTE::
	This will eventually be a base class for the "XYZShape" class
	the is a special type of grid...one that contains a 'shape'
	inside a general XYZ grid...this class acts only as a function
	generator

	****if you want to add different shapes...you MUST include these two functions
    bool ShapeFunc()(double x, double y, double z)
   		     and
	bool ShapeFunc()(coord<> pos)

	the next functions these hold ARE OPTIONAL
	BUT if you wish to use these function generators for
	'Gradient Functions' or just general 'other' functional applications
	put these three operators in the mix

	void operator()(Alteritem, coord<>); # this alters the 'alteritem' based on its posistion
	bool operator()(coord<>); # like the 'ShapeFunc' function above, but as an operator


	these function take in a single xyz point return true if the point is within the
	functional shape, and false if not....
*/





//------------------------------------------------------------


//---------------The Non-Shape--------------------
class NullXYZ{};

//-----------Full-----------------------
/* This one is a 'null' type shape...it simply allows the use of ALL the grid points it is asked about
  ***it is not advisable to use this one****
  ...simply use the UniformGrid Grid type instead....
  it will save memmory and a bit of initialization speed
*/



class XYZfull_ : public GeneralXYZShape {

	public:
		XYZfull_():
			GeneralXYZShape()
		{}
		XYZfull_(const XYZfull_ &cp):
			GeneralXYZShape()
		{}
		XYZfull_(const coord<> &l1, const coord<> &l2=ZeroType<coord<> >::zero(),const coord<> &l3=ZeroType<coord<> >::zero()): //a dummy
			GeneralXYZShape()
		{}
		inline bool ShapeFunc(double x, double y, double z){ return true;	}
		inline bool ShapeFunc(coord<> &xyz){	return true;	}
		inline bool operator()(coord<> &test){	return true;	}

		template<class Alter_t>
		inline void operator()(Alter_t &in, coord<> &pos){		}

		XYZfull_ &operator=(const XYZfull_ &rhs);
};

typedef XYZparts<XYZfull_> XYZfull;


//-----------End Full-----------------------


//-----------Rectangular solid--------------
class XYZrect_: public GeneralXYZShape {
	private:


	public:
		XYZrect_(){}

		XYZrect_(const XYZrect_ &copy):
			GeneralXYZShape(copy)
		{}

		XYZrect_(const coord<> &minin,const coord<> &maxin ,const coord<> &center=ZeroType<coord<> >::zero()):
			GeneralXYZShape(minin, maxin,center)
		{}

		XYZrect_(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
			GeneralXYZShape(xmin,xmax,ymin,ymax,zmin,zmax)
		{}

		XYZrect_ &operator=(const XYZrect_ &rhs);

		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}

};

typedef XYZparts<XYZrect_> XYZrect;

//-----------End Rectangular solid--------------



//-----------Cylinder-----------------------------


class XYZcylinder_: public GeneralXYZShape {

	public:
		XYZcylinder_(){}

		XYZcylinder_(const XYZcylinder_ &copy):
			GeneralXYZShape(copy)
		{}

		XYZcylinder_(const coord<> &minin,const coord<> &maxin):
			GeneralXYZShape(minin, maxin, coord<>(0))
		{}

		XYZcylinder_(const coord<> &minin,const coord<> &maxin,const coord<> &ce):
			GeneralXYZShape(minin, maxin, coord<>(0))
		{}


		XYZcylinder_(double rmin, double rmax, double phimin, double phimax, double zmin, double zmax):
			GeneralXYZShape(rmin,rmax,phimin,phimax,zmin,zmax)
		{}

		XYZcylinder_ &operator=(const XYZcylinder_ &rhs);

		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}

};

typedef XYZparts<XYZcylinder_> XYZcylinder;


//----------End Cylinder-----------------------

//----------- The Basic 'Plane' slicer class'-----------
// do not use this expeicitly..use the 'directional' vesions..X,Y,Z or 3D

class XYZBasicPlane : public GeneralXYZShape
{
	private:
		double m_;	//slope
		double b_;	//intercept
		double sc_;
		bool above_;	//choose 'above' the plane or 'below' plane

	public:

		XYZBasicPlane():
			GeneralXYZShape(),
			m_(1), b_(0), sc_(1),above_(true)
		{}

		XYZBasicPlane(double mx, double bx):
			GeneralXYZShape(),
			m_(mx), b_(bx), sc_(1),above_(true)
		{}

		XYZBasicPlane(double mx, double bx, double sc):
			GeneralXYZShape(),
			m_(mx), b_(bx), sc_(sc),above_(true)
		{}

		XYZBasicPlane(const XYZBasicPlane &cp):
			GeneralXYZShape(),
			m_(cp.m_), b_(cp.b_), sc_(cp.sc_),above_(cp.above_)
		{}

		inline double &m(){	return m_;	}
		inline double &b(){	return b_;	}
		inline double &Scale(){	return sc_;	}
		inline bool &Above(){	return above_;	}
		inline double &scale(){	return sc_;	}
		inline bool &above(){	return above_;	}

		inline double m() const {	return m_;	}
		inline double b() const {	return b_;	}
		inline double Scale() const {	return sc_;	}
		inline bool Above() const {	return above_;	}
		inline double scale() const {	return sc_;	}
		inline bool above() const {	return above_;	}

		inline void Setm(double mm){	m_=mm;	}
		inline void Setb(double bb){	b_=bb;	}
		inline void SetScale(double sc){	sc_=sc;	}
		inline void SetAbove(){	above_=true;	}
		inline void SetBelow(){	above_=false;	}

		XYZBasicPlane &operator=(const XYZBasicPlane &rhs);

};

//----------- END The Basic 'Plane' slicer class'-----------

//----------- the XY-test-Plane
class XYZplaneXY_ : public XYZBasicPlane
{

	public:

		XYZplaneXY_():
			XYZBasicPlane()
		{}

		XYZplaneXY_(double mx, double bx):
			XYZBasicPlane(mx,bx)
		{}

		XYZplaneXY_(double mx, double bx, double sc):
			XYZBasicPlane(mx,bx,sc)
		{}

		XYZplaneXY_(const XYZplaneXY_ &cp):
			XYZBasicPlane(cp)
		{}

		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}

	// the altering driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.x()*m()+b());
		}

		XYZplaneXY_ &operator=(const XYZplaneXY_ &rhs);

};

typedef XYZparts<XYZplaneXY_> XYZplaneXY;

//----------- the YZ-test-Plane
class XYZplaneYZ_ : public XYZBasicPlane
{

	public:

		XYZplaneYZ_():
			XYZBasicPlane()
		{}

		XYZplaneYZ_(double mx, double bx):
			XYZBasicPlane(mx,bx)
		{}

		XYZplaneYZ_(double mx, double bx, double sc):
			XYZBasicPlane(mx,bx,sc)
		{}

		XYZplaneYZ_(const XYZplaneYZ_ &cp):
			XYZBasicPlane(cp)
		{}

		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}

	// the altering driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.y()*m()+b());
		}

		XYZplaneYZ_ &operator=(const XYZplaneYZ_ &rhs);
};

typedef XYZparts<XYZplaneYZ_> XYZplaneYZ;


//----------- the XZ-test-Plane
class XYZplaneXZ_ : public XYZBasicPlane
{

	public:

		XYZplaneXZ_():
			XYZBasicPlane()
		{}

		XYZplaneXZ_(double mx, double bx):
			XYZBasicPlane(mx,bx)
		{}

		XYZplaneXZ_(double mx, double bx, double sc):
			XYZBasicPlane(mx,bx,sc)
		{}

		XYZplaneXZ_(const XYZplaneXZ_ &cp):
			XYZBasicPlane(cp)
		{}

		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}

	// the altering driver
		template<class Alter_t>
		inline void operator()(Alter_t &rhs, coord<> &pos)
		{
			rhs+=Scale()*(pos.z()*m()+b());
		}

		XYZplaneXZ_ &operator=(const XYZplaneXZ_ &rhs);

};

typedef XYZparts<XYZplaneXZ_> XYZplaneXZ;

//--------------Basic 3D Plane slicer--------------
// 3D linear grad types
class XYZBasic3Dplane : public GeneralXYZShape
{
	private:
		bool above_;

	public:

		XYZBasic3Dplane():
			GeneralXYZShape(),
			above_(true)
		{}

		XYZBasic3Dplane(coord<> mx, coord<> bx):
			GeneralXYZShape(mx,bx),
			above_(true)
		{}

		XYZBasic3Dplane(coord<> mx, coord<> bx, coord<> sc):
			GeneralXYZShape(mx,bx,sc),
			above_(true)
		{}

		XYZBasic3Dplane(const XYZBasic3Dplane &cp):
			GeneralXYZShape(cp),
			above_(cp.above_)
		{}

		inline coord<> &m(){	return min();	}
		inline coord<> &b(){	return max();	}
		inline coord<> &Scale(){	return center();	}
		inline coord<> &scale(){	return center();	}
		inline bool &Above(){	return above_;	}
		inline bool &above(){	return above_;	}

		inline coord<> m() const {	return min();	}
		inline coord<> b() const {	return max();	}
		inline coord<> Scale() const {	return center();	}
		inline coord<> scale() const {	return center();	}
		inline bool Above() const {	return above_;	}
		inline bool above() const {	return above_;	}

		inline void Setm(coord<> mm){	min()=mm;	}
		inline void Setb(coord<> bb){	max()=bb;	}
		inline void SetScale(coord<> sc){	center()=sc;	}
		inline void SetAbove(){	above_=true;	}
		inline void SetBelow(){	above_=false;	}

		XYZBasic3Dplane &operator=(const XYZBasic3Dplane &rhs);

		XYZBasic3Dplane &operator=(const XYZBasicPlane &rhs);

};


// -----------3D Plane-----------------------
// 3D linear grad types
class XYZ3Dplane_ : public XYZBasic3Dplane
{

	public:

		XYZ3Dplane_():
			XYZBasic3Dplane()
		{}

		XYZ3Dplane_(coord<> mx, coord<> bx):
			XYZBasic3Dplane(mx,bx)
		{}

		XYZ3Dplane_(coord<> mx, coord<> bx, coord<> sc):
			XYZBasic3Dplane(mx,bx,sc)
		{}

		XYZ3Dplane_(XYZ3Dplane_ &cp):
			XYZBasic3Dplane(cp)
		{}

		XYZ3Dplane_ &operator=(const XYZ3Dplane_ &rhs);

	//the shape functions
		bool ShapeFunc(double x, double y, double z);
		bool ShapeFunc(coord<> &xyz);
		inline bool operator()(coord<> &rhs){	return ShapeFunc(rhs);	}


	//the altering func
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


typedef XYZparts<XYZ3Dplane_> XYZ3Dplane;

END_BL_NAMESPACE



#endif


