/* xyzshape.h ********/


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
 	xyzshape.h-->A grid in the 'uniform' way (i.e. rectangular x,y,z) that can
 	take a plethera of shape functions....these shape functions CHOP down the
 	large base class (unifom grid) into something that the 'XYZshapeFunc' specifies
 	the iterators over this grid then ONLY work over the validated grid points.

 	this grid type is probably the most usefull of all the grid peices...

 	one can even input shape exrpressions to calculate an abitrary grid
 	(i.e. Plane||Cylinder, or, "rect&&plane" etc

 	the use the most general form/ arbitrary shape use

 		XYZshape<> moo(Grid<UniformGrid>(mins, max, dims)
 		moo.CalcShape("Plane && Rect");

 	where 'Plane' and 'Rect' are the 'XYZshapeparts.h' classes

 	(or one that you have written yourself...)


 */


#ifndef _XYZshape_h_
#define _XYZshape_h_ 1


#include "container/complex.h"
#include "container/range.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"
#include "container/grids/grids.h"
#include "container/grids/xyzgenshape.h"
#include "container/grids/xyzshapeparts.h"


BEGIN_BL_NAMESPACE


/* An XYZ Rectangular Grid with an embeded shape....
  Things to NOTE::
  	1) the "++" operator loops through the contained SHAPE NOT the whole grid
  	2) to loop through the entire grid use "XYZshape::Griditerator"

*/

//-------------------The 'list' of shapes class..up to 10 are allowed-----------------
// HOWEVER...these shapes here are ALWAYS 'ANDED' i.e. a point MUST BE in ALL
//  ten sub shapes...it is better to use the "XYZand" and "XYZor" pair
// functions to create a complex shape becuase you have Much more control
// over the design of your shape....


//pre dec for iterator
template<class ShapeEng_t>
class XYZshapeIter;


template<class Shape_t1=XYZfull, class Shape_t2=NullXYZ, class Shape_t3=NullXYZ,
		 class Shape_t4=NullXYZ, class Shape_t5=NullXYZ, class Shape_t6=NullXYZ,
		 class Shape_t7=NullXYZ,  class Shape_t8=NullXYZ, class Shape_t9=NullXYZ,
		 class Shape_t10=NullXYZ>
class XYZshape{};

template<>
class XYZshape<XYZfull,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public Grid<UniformGrid>
{

	private:
		Vector<coord<> > data_;			//the point inside the uniform grid of the {i,j,k} pos
		typedef double numtype;

		void init();

		template<class ShapeExpr>
		void init(ShapeExpr in)
		{
			data_.resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(in.ShapeFunc(myit.Point()))
				{
					data_(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data_.resizeAndPreserve(ct);
			if(data_.size()<=0) EmptyWarn();
		}

	public:
		typedef XYZshape<XYZfull,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> Shape;

		friend class XYZshapeIter<Shape > ;

		typedef XYZshapeIter<Shape > iterator;

		typedef Grid<UniformGrid>::iterator Griditerator;

		const void EmptyWarn();



		XYZshape():
			Grid<UniformGrid>(),
			data_(0,0)
		{}

		XYZshape(const XYZshape &cp):
			Grid<UniformGrid>(cp),
			data_(cp.data_)
		{}

		XYZshape(const Grid<UniformGrid> &mygrid):
			Grid<UniformGrid>(mygrid),
			data_(0,0)
		{}

		XYZshape(const coord<> &mins, const coord<> &maxs, const coord<int> &dims):
			Grid<UniformGrid>(mins, maxs, dims),
			data_(0,0)
		{}


		template<class ShapeExpr>
		XYZshape(const Grid<UniformGrid> &mygrid, const ShapeExpr &expr):
			Grid<UniformGrid>(mygrid),
			data_(0,0)
		{
			calculate(expr);
		}

		template<class ShapeExpr>
		XYZshape(const coord<> &mins, const coord<> &maxs, const coord<int> &dims,const ShapeExpr &expr):
			Grid<UniformGrid>(mins, maxs, dims),
			data_(0,0)
		{
			calculate(expr);
		}

		XYZshape(const Grid<UniformGrid> &mygrid, const XYZfull &expr):
			Grid<UniformGrid>(mygrid),
			data_(0,0)
		{
			init();
		}

		XYZshape &operator=(const XYZshape &rhs);

	//THE MAIN DRIVER...MUST BE CALLED AFTER DELCATATION!!!
		template<class ShapeExpr>
		inline void CalcShape(const ShapeExpr &in)
		{
			init(in);
		}

	//THE MAIN DRIVER...MUST BE CALLED AFTER DELCATATION!!!
		template<class ShapeExpr>
		inline void calculate(const ShapeExpr &in)
		{
			init(in);
		}

		inline bool empty(){	return data_.size()<=0;	}

		inline double x(int i)	{ return data_(i).x();	}
		inline double y(int i)	{ return data_(i).y();	}
		inline double z(int i)	{ return data_(i).z();	}

		inline coord<> Point(int i, int j, int k)
		{
			return coord<>(data_(i).x(), data_(j).y(), data_(k).z());;
		}

		inline coord<> &operator()(int i)
		{
			return data_(i);
		}

		inline coord<> &Point(int i)
		{
			return data_(i);
		}

		inline coord<> operator()(int i, int j, int k)
		{
			return coord<>(data_(i).x(), data_(j).y(), data_(k).z());
		}

		inline coord<> GridPoint(int i, int j, int k)
		{
			return Grid<UniformGrid>::Point(i,j,k);
		}

		inline coord<> gridPoint(int i, int j, int k)
		{
			return Grid<UniformGrid>::Point(i,j,k);
		}

	//for time dependant shapes
		inline coord<> Point(int i, int j, int k, double t)
		{
			return coord<>(data_(i).x(), data_(j).y(), data_(k).z());;
		}

		inline coord<> &operator()(int i, double t)
		{
			return data_(i);
		}

		inline coord<> &Point(int i, double t)
		{
			return data_(i);
		}

		inline coord<> operator()(int i, int j, int k, double t)
		{
			return coord<>(data_(i).x(), data_(j).y(), data_(k).z());
		}

		inline coord<> GridPoint(int i, int j, int k, double t)
		{
			return Grid<UniformGrid>::Point(i,j,k);
		}

		inline coord<> gridPoint(int i, int j, int k, double t)
		{
			return Grid<UniformGrid>::Point(i,j,k);
		}

		inline Vector<coord<> > &data(){	return data_;	}
		inline coord<>  &data(int i){	return data_(i);	}

		inline int size() const	{	return data_.size();	}
		inline int GridSize() const {	return  Grid<UniformGrid>::size();	}
		inline int gridSize() const {	return  GridSize();	}

		void PrintGrid(std::ostream &oo);

		void PrintShape(std::ostream &oo);

		inline void printGrid(std::ostream &oo)
		{	PrintShape(oo);	}

		inline void printShape(std::ostream &oo)
		{	PrintGrid(oo);	}

};

//typedef to save me some oodles of typing
typedef XYZshape<XYZfull,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> NonShape;



//------------------- Single Shape-------------------
template<class ShapeEng_t1>
class XYZshape<ShapeEng_t1,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public NonShape
{

	private:
		ShapeEng_t1 ins1_;
		void init()
		{
			data().resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(InShape(myit.Point()))
				{
					data()(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data().resizeAndPreserve(ct);
			if(data().size()<=0) EmptyWarn();
		}

	public:
		typedef XYZshape<ShapeEng_t1,NullXYZ,NullXYZ,NullXYZ,NullXYZ,
							NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> Shape;

		friend class XYZshapeIter<Shape >;

		typedef XYZshapeIter<XYZshape<ShapeEng_t1,NullXYZ,NullXYZ,NullXYZ,NullXYZ,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> > iterator;
//the many constructors....

		XYZshape():
			NonShape(),
			ins1_()
		{}

		XYZshape(const XYZshape &cp):
			NonShape(cp),
			ins1_(cp.ins1_)
		{}

		XYZshape(const Grid<UniformGrid> &mygrid):
			NonShape(mygrid),
			ins1_()
		{}

		XYZshape(const Grid<UniformGrid> &mygrid,const ShapeEng_t1 &myshape ):
			NonShape(mygrid),
			ins1_(myshape)
		{
			init();
		}

		~XYZshape()
		{}

		inline void ReCalc()
		{
			init();
		}


		void SetShape(ShapeEng_t1 &in)
		{
			ins1_=in;
			init();
		}

		void SetGrid(Grid<UniformGrid> &in)
		{
			Grid<UniformGrid>::operator=(in);
			init();
		}
		void setShape(ShapeEng_t1 &in)
		{	SetShape(in);	}

		void setGrid(Grid<UniformGrid> &in)
		{	setGrid(in);	}


		inline bool InShape(coord<> &inpoint)
		{
			return ins1_.ShapeFunc(inpoint);
		}


		inline bool InShape(double x, double y, double z)
		{
			return ins1_.ShapeFunc(x,y,z);
		}

		inline bool inShape(coord<> &inpoint)
		{	return InShape(inpoint);	}

		inline bool inShape(double x, double y, double z)
		{	return InShape( x,  y,  z);	}


		XYZshape &operator=(const XYZshape &rhs)
		{
			if(this==&rhs) return *this;
			ins1_=rhs.ins1_;
			NonShape::operator=(rhs);
			return *this;
		}



};

//------------------- Double Shape-------------------
template<class ShapeEng_t1, class ShapeEng_t2>
class XYZshape<ShapeEng_t1,ShapeEng_t2,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public NonShape
{

	private:
		ShapeEng_t1 ins1_;
		ShapeEng_t2 ins2_;
		void init()
		{
			data().resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(InShape(myit.Point()))
				{
					data()(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data().resizeAndPreserve(ct);
			if(data().size()<=0) EmptyWarn();
		}


	public:
		typedef XYZshape<ShapeEng_t1,ShapeEng_t2,NullXYZ,NullXYZ,NullXYZ,
						 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> Shape;

		friend class XYZshapeIter<Shape >;


		typedef XYZshapeIter<XYZshape<ShapeEng_t1,ShapeEng_t2,NullXYZ,NullXYZ,NullXYZ,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> > iterator;
//the many constructors....

		XYZshape():
			NonShape(),
			ins1_(),
			ins2_()
		{}

		XYZshape(XYZshape &cp):
			NonShape(cp),
			ins1_(cp.ins1_),
			ins2_(cp.ins2_)
		{}

		XYZshape(Grid<UniformGrid> &mygrid, ShapeEng_t1 &in1, ShapeEng_t2 &in2):
			NonShape(mygrid),
			ins1_(in1),
			ins2_(in2)
		{
			init();
		}

		~XYZshape()
		{}

		inline void ReCalc()
		{
			init();
		}

		void SetShape(ShapeEng_t1 &in, ShapeEng_t2 &in2)
		{
			ins1_=in;
			ins2_=in;
			init();
		}

		void SetShape1(ShapeEng_t1 &in)
		{
			ins1_=in;
			init();
		}

		void SetShape2(ShapeEng_t2 &in)
		{
			ins2_=in;
			init();
		}

		inline bool InShape(coord<> &inpoint)
		{
			return ins1_.ShapeFunc(inpoint) &&
					ins2_.ShapeFunc(inpoint);
		}

		inline bool InShape(double x, double y, double z)
		{
			return ins1_->ShapeFunc(x,y,z) &&
					ins2_->ShapeFunc(x,y,z);
		}
		inline bool inShape(coord<> &inpoint)
		{	return InShape(inpoint);	}

		inline bool inShape(double x, double y, double z)
		{	return InShape( x,  y,  z);	}


		XYZshape &operator=(const XYZshape &rhs)
		{
			if(this==&rhs) return *this;
			ins1_=rhs.ins1_;
			ins2_=rhs.ins2_;
			NonShape::operator=(rhs);
			return *this;
		}

};

//------------------- Triple Shape-------------------
template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3>
class XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public NonShape
{

	private:
		ShapeEng_t1 ins1_;
		ShapeEng_t2 ins2_;
		ShapeEng_t3 ins3_;

		void init()
		{
			data().resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(InShape(myit.Point()))
				{
					data()(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data().resizeAndPreserve(ct);
			if(data().size()<=0) EmptyWarn();
		}


	public:
		typedef XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,NullXYZ,NullXYZ,
					 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> Shape;

		friend class XYZshapeIter<Shape >;
		typedef XYZshapeIter<XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,NullXYZ,NullXYZ,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> > iterator;
//the many constructors....

		XYZshape():
			NonShape(),
			ins1_(),
			ins2_(),
			ins3_()
		{}

		XYZshape(XYZshape &cp):
			NonShape(cp),
			ins1_(cp.ins1_),
			ins2_(cp.ins2_),
			ins3_(cp.ins3_)
		{}

		XYZshape(Grid<UniformGrid> &mygrid, ShapeEng_t1 &in1, ShapeEng_t2 &in2, ShapeEng_t3 &in3):
			NonShape(mygrid),
			ins1_(in1),
			ins2_(in2),
			ins3_(in3)
		{
			init();
		}

		~XYZshape()
		{}

		inline void ReCalc()
		{
			init();
		}

		void SetShape(ShapeEng_t1 &in, ShapeEng_t2 &in2, ShapeEng_t3 &in3)
		{
			ins1_=in;
			ins2_=in2;
			ins3_=in3;
			init();
		}

		void SetShape1(ShapeEng_t1 &in)
		{
			ins1_=in;
			init();
		}

		void SetShape2(ShapeEng_t2 &in)
		{
			ins2_=in;
			init();
		}

		void SetShape3(ShapeEng_t3 &in)
		{
			ins3_=in;
			init();
		}

		inline bool InShape(coord<> &inpoint)
		{
			return ins1_.ShapeFunc(inpoint) &&
					ins2_.ShapeFunc(inpoint) &&
					ins3_.ShapeFunc(inpoint);
		}

		inline bool InShape(double x, double y, double z)
		{
			return ins1_.ShapeFunc(x,y,z) &&
					ins2_.ShapeFunc(x,y,z) &&
					ins3_.ShapeFunc(x,y,z);
		}

		inline bool inShape(coord<> &inpoint)
		{	return InShape(inpoint);	}

		inline bool inShape(double x, double y, double z)
		{	return InShape( x,  y,  z);	}

		XYZshape &operator=(const XYZshape &rhs)
		{
			if(this==&rhs) return *this;
			ins1_=rhs.ins1_;
			ins2_=rhs.ins2_;
			ins3_=rhs.ins3_;
			NonShape::operator=(in);
			return *this;
		}


};


//------------------- quad Shape-------------------
template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3, class ShapeEng_t4>
class XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public NonShape
{

	private:
		ShapeEng_t1 ins1_;
		ShapeEng_t2 ins2_;
		ShapeEng_t3 ins3_;
		ShapeEng_t4 ins4_;

		void init()
		{
			data().resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(InShape(myit.Point()))
				{
					data()(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data().resizeAndPreserve(ct);
			if(data().size()<=0) EmptyWarn();
		}


	public:

		typedef XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,NullXYZ,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> Shape;

		friend class XYZshapeIter<Shape >;

		typedef XYZshapeIter<Shape > iterator;
//the many constructors....

		XYZshape():
			NonShape(),
			ins1_(),
			ins2_(),
			ins3_(),
			ins4_()
		{}

		XYZshape(XYZshape &cp):
			NonShape(cp),
			ins1_(cp.ins1_),
			ins2_(cp.ins2_),
			ins3_(cp.ins3_),
			ins4_(cp.ins4_)
		{}

		XYZshape(Grid<UniformGrid> &mygrid, ShapeEng_t1 &in1, ShapeEng_t2 &in2, ShapeEng_t3 &in3,ShapeEng_t4 &in4):
			NonShape(mygrid),
			ins1_(in1),
			ins2_(in2),
			ins3_(in3),
			ins4_(in4)
		{
			init();
		}

		~XYZshape()
		{}

		inline void ReCalc()
		{
			init();
		}

		void SetShape(ShapeEng_t1 &in, ShapeEng_t2 &in2, ShapeEng_t3 &in3,ShapeEng_t4 &in4)
		{
			ins1_=in;
			ins2_=in2;
			ins3_=in3;
			ins4_=in4;
			init();
		}

		void SetShape1(ShapeEng_t1 &in)
		{
			ins1_=in;
			init();
		}

		void SetShape2(ShapeEng_t2 &in)
		{
			ins2_=in;
			init();
		}

		void SetShape3(ShapeEng_t3 &in)
		{
			ins3_=in;
			init();
		}

		void SetShape4(ShapeEng_t4 &in)
		{
			ins4_=in;
			init();
		}

		inline bool InShape(coord<> &inpoint)
		{
			return ins1_.ShapeFunc(inpoint) &&
					ins2_.ShapeFunc(inpoint) &&
					ins3_.ShapeFunc(inpoint) &&
					ins4_.ShapeFunc(inpoint);
		}

		inline bool InShape(double x, double y, double z)
		{
			return ins1_.ShapeFunc(x,y,z) &&
					ins2_.ShapeFunc(x,y,z) &&
					ins3_.ShapeFunc(x,y,z) &&
					ins4_.ShapeFunc(x,y,z);
		}

		inline bool inShape(coord<> &inpoint)
		{	return InShape(inpoint);	}

		inline bool inShape(double x, double y, double z)
		{	return InShape( x,  y,  z);	}

		XYZshape &operator=(const XYZshape &rhs)
		{
			if(this==&rhs) return *this;
			ins1_=rhs.ins1_;
			ins2_=rhs.ins2_;
			ins3_=rhs.ins3_;
			ins4_=rhs.ins4_;
			NonShape::operator=(in);
			return *this;
		}

};

//------------------- quad Shape-------------------
template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3, class ShapeEng_t4,class ShapeEng_t5>
class XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,ShapeEng_t5,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ>:
	public NonShape
{
	friend class XYZshapeIter<XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,ShapeEng_t5,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> >;

	private:
		ShapeEng_t1 ins1_;
		ShapeEng_t2 ins2_;
		ShapeEng_t3 ins3_;
		ShapeEng_t4 ins4_;
		ShapeEng_t5 ins5_;

		void init()
		{
			data().resize(Grid<UniformGrid>::size());

			//Grid<UniformGrid>::reset();
			Griditerator myit(*this);
			int ct=0;
			while(myit)
			{
				if(InShape(myit.Point()))
				{
					data()(ct)=(myit.Point());
					ct++;
				}
				++myit;
			}
			data().resizeAndPreserve(ct);
			if(data().size()<=0) EmptyWarn();
		}


	public:

		typedef XYZshapeIter<XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,NullXYZ,
				 	NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> > iterator;
//the many constructors....

		XYZshape():
			NonShape(),
			ins1_(),
			ins2_(),
			ins3_(),
			ins4_(),
			ins5_()
		{}

		XYZshape(XYZshape &cp):
			NonShape(cp),
			ins1_(cp.ins1_),
			ins2_(cp.ins2_),
			ins3_(cp.ins3_),
			ins4_(cp.ins4_),
			ins5_(cp.ins5_)
		{}

		XYZshape(Grid<UniformGrid> &mygrid, ShapeEng_t1 &in1, ShapeEng_t2 &in2,
				ShapeEng_t3 &in3,ShapeEng_t4 &in4, ShapeEng_t5 &in5):
			NonShape(mygrid),
			ins1_(&in1),
			ins2_(&in2),
			ins3_(&in3),
			ins4_(&in4),
			ins5_(&in5)
		{
			init();
		}

		~XYZshape()
		{}

		inline void ReCalc()
		{
			init();
		}

		void SetShape(ShapeEng_t1 &in, ShapeEng_t2 &in2, ShapeEng_t3 &in3,ShapeEng_t4 &in4,ShapeEng_t5 &in5)
		{
			ins1_=in;
			ins2_=in2;
			ins3_=in3;
			ins4_=in4;
			ins5_=in5;
			init();
		}

		void SetShape1(ShapeEng_t1 &in)
		{
			ins1_=in;
			init();
		}

		void SetShape2(ShapeEng_t2 &in)
		{
			ins2_=in;
			init();
		}

		void SetShape3(ShapeEng_t3 &in)
		{
			ins3_=in;
			init();
		}

		void SetShape4(ShapeEng_t4 &in)
		{
			ins4_=in;
			init();
		}

		void SetShape5(ShapeEng_t5 &in)
		{
			ins5_=in;
			init();
		}


		inline bool InShape(coord<> &inpoint)
		{
			return ins1_.ShapeFunc(inpoint) &&
					ins2_.ShapeFunc(inpoint) &&
					ins3_.ShapeFunc(inpoint) &&
					ins4_.ShapeFunc(inpoint) &&
					ins5_.ShapeFunc(inpoint);
		}

		inline bool InShape(double x, double y, double z)
		{
			return ins1_.ShapeFunc(x,y,z) &&
					ins2_.ShapeFunc(x,y,z) &&
					ins3_.ShapeFunc(x,y,z) &&
					ins4_.ShapeFunc(x,y,z) &&
					ins5_.ShapeFunc(x,y,z);
		}

		inline bool inShape(coord<> &inpoint)
		{	return InShape(inpoint);	}

		inline bool inShape(double x, double y, double z)
		{	return InShape( x,  y,  z);	}

		XYZshape &operator=(const XYZshape &rhs)
		{
			if(this==&rhs) return *this;
			ins1_=rhs.ins1_;
			ins2_=rhs.ins2_;
			ins3_=rhs.ins3_;
			ins4_=rhs.ins4_;
			ins5_=rhs.ins5_;
			NonShape::operator=(in);
			return *this;
		}


};



template<class ShapeEng_t=NonShape>
class XYZshapeIter {
	private:
		ShapeEng_t *mys_;
		bool shapeend_;				//bool for the end of the shape loop
		bool gridend_;				//bool for the end of the grid loop
		Range range_;

		int iter_; 	//iteration loop counters for the SHAPE


		coord<> *outPt;				//a static point so we do not have to initiallize a point over and overagain...
	public:

		XYZshapeIter():
			mys_(NULL),
			shapeend_(false),
			iter_(0),
			outPt(NULL)
		{}

		XYZshapeIter(ShapeEng_t &in):
			mys_(&in),
			shapeend_(true),
			range_(0, Range::End),
			iter_(0)
		{
			if(mys_){
				outPt=&(mys_->data()(0));
			}else{
				shapeend_=false;
			}

		}

		XYZshapeIter(ShapeEng_t *in):
			mys_(in),
			shapeend_(true),
			range_(0, Range::End),
			iter_(0),
			outPt(NULL)
		{

			if(in){
				outPt=&(mys_->data()(0));
			}else{
				shapeend_=false;
			}
		}

		XYZshapeIter(ShapeEng_t &in, Range R):
			mys_(&in),
			shapeend_(true),
			range_(R),
			iter_(R.first(0))
		{
			if(mys_){
				outPt=&(mys_->data()(R.first(0)));
			}else{
				shapeend_=false;
			}

		}

		XYZshapeIter(ShapeEng_t *in, Range R):
			mys_(in),
			shapeend_(true),
			range_(R),
			iter_(R.first(0)),
			outPt(NULL)
		{

			if(in){
				outPt=&(mys_->data()(R.first(0)));
			}else{
				shapeend_=false;
			}
		}

		XYZshapeIter(ShapeEng_t &in, int begin, int end=Range::End, int stride=1)  :
			mys_(&in),
			shapeend_(true),
			range_(begin, end, stride),
			iter_(0)
		{
			if(mys_){
				outPt=&(mys_->data()(range_.first(iter_)));
			}else{
				shapeend_=false;
			}

		}

		XYZshapeIter(ShapeEng_t *in, int begin, int end=Range::End, int stride=1):
			mys_(in),
			shapeend_(true),
			range_(begin, end, stride),
			iter_(0),
			outPt(NULL)
		{

			if(in){
				outPt=&(mys_->data()(range_.first(iter_)));
			}else{
				shapeend_=false;
			}
		}

		~XYZshapeIter()
		{
			if(mys_)	mys_=NULL;
			if(outPt)	outPt=NULL;
		}

		void SetShape(ShapeEng_t &in)
		{	mys_=&in;		}

		void operator++()
		{
			if(!shapeend_ || !mys_) mys_->Additerr();
			if(range_(iter_)<range_.last(mys_->size())-1){
				++iter_;
				outPt=&(mys_->data()(range_(iter_)));
				return;
			}else{
				shapeend_=false;
			}
		}

		void operator++(int)
		{
			if(!shapeend_ || !mys_) mys_->Additerr();
			if(range_(iter_)<range_.last(mys_->size())-1){
				++iter_;
				outPt=&(mys_->data()(range_(iter_)));
				return;
			}else{
				shapeend_=false;
			}
		}

		inline operator bool()	{ return shapeend_;	}

		inline void begin(int i)
		{
			RunTimeAssert(i<range_.last(mys_->size())-1);
			iter_=0;
			range_=Range(i, range_.last(mys_->size()), range_.stride());
		}

		inline int begin()const{	return range_.first(0);	}


		inline int end()const{		return range_.last(mys_->size());	}
		inline void end(int i)
		{
			RunTimeAssert(i<mys_->size());
			range_=Range(range_.first(0), i, range_.stride());
		}

		void reset()
		{
			shapeend_=true;
			iter_=0;
			if(mys_){
				outPt=&(mys_->data()(range_.first(iter_)));
			}else{
				shapeend_=false;
			}

		}

		inline coord<> &Point()
		{
			return *outPt;
			//return outPt;
		}

		inline coord<> &Point(double t)
		{
			return *outPt;
			//return outPt;
		}

		inline int curpos(){	return range_(iter_);	}

		inline double x()	{ return outPt->x();	}
		inline double y()	{ return outPt->y();	}
		inline double z()	{ return outPt->z();	}

		inline double x(double t)	{ return outPt->x();	}
		inline double y(double t)	{ return outPt->y();	}
		inline double z(double t)	{ return outPt->z();	}

		inline double dx()	{ return mys_->dx();	}
		inline double dy()	{ return mys_->dy();	}
		inline double dz()	{ return mys_->dz();	}
		inline coord<> dr()	{ return mys_->dr();	}

		inline double CellVolume()	{ return mys_->CellVolume();	}
		inline double cellVolume()	{ return mys_->cellVolume();	}

		XYZshapeIter &operator=(XYZshapeIter &rhs)
		{
			if(*this==rhs) return *this;
			mys_=rhs.mys_;
			range_=rhs.range_;
			shapeend_=rhs.shapeend_;
			iter_=rhs.iter_;
			return *this;
		}
};


std::ostream &operator<<(std::ostream &oo, XYZshape<XYZfull,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out);
template<class ShapeEng_t1>
std::ostream &operator<<(std::ostream &oo, XYZshape<ShapeEng_t1,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out)
{
	out.PrintShape(oo);
	return oo;
}

template<class ShapeEng_t1, class ShapeEng_t2>
std::ostream &operator<<(std::ostream &oo, XYZshape<ShapeEng_t1,ShapeEng_t2,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out)
{
	out.PrintShape(oo);
	return oo;
}

template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3>
std::ostream &operator<<(std::ostream &oo, XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out)
{
	out.PrintShape(oo);
	return oo;
}

template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3, class ShapeEng_t4>
std::ostream &operator<<(std::ostream &oo, XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out)
{
	out.PrintShape(oo);
	return oo;
}

template<class ShapeEng_t1, class ShapeEng_t2, class ShapeEng_t3, class ShapeEng_t4,class ShapeEng_t5>
std::ostream &operator<<(std::ostream &oo, XYZshape<ShapeEng_t1,ShapeEng_t2,ShapeEng_t3,ShapeEng_t4,ShapeEng_t5,NullXYZ,NullXYZ,NullXYZ,NullXYZ,NullXYZ> &out)
{
	out.PrintShape(oo);
	return oo;
}
//#include "container/grids/xyzshape_meth.h"

END_BL_NAMESPACE


#endif

