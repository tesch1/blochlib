

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
 	fullcubegrid.h-->part of the grids classes a full cube grid
 */


#ifndef _fullcubegrid_h_
#define _fullcubegrid_h_


/* Cude Grid Generator Class
	generates a cube of points-->each dimension is the same size
	and each direction has the same maxes and mins... (i.e. a cube)

*/

#include <iostream>
#include <string>
#include "container/grids/gengrid.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


class FullCubeGridIter;

class FullCubeGrid : public GeneralGrid {
	friend class FullCubeGridIter;
	private:
		Vector<coord<> > data_;		//holds the coord data

	public:

		typedef FullCubeGridIter iterator;

		FullCubeGrid(): GeneralGrid() {}

		FullCubeGrid(int dims):
			GeneralGrid(dims,dims,dims), data_(dims,0)
		{
			init();
		}

		FullCubeGrid(const FullCubeGrid &rhs):
			GeneralGrid(rhs), data_(rhs.data_){}

		FullCubeGrid(double min, double max, int dims):
			GeneralGrid(coord<>(min),coord<>(max),coord<int>(dims)), data_(dims,0)
		{
			init();
		}

		~FullCubeGrid(){}

		void init();
		coord<> Point(int i, int j, int k) ;
		coord<> Point(int i, int j, int k, double t) ;

		//dimension extraction
		inline double Xdim()  { return GeneralGrid::Xdim();	}
		inline double Ydim()  { return GeneralGrid::Ydim();	};
		inline double Zdim() { return GeneralGrid::Zdim();	};
		inline double Rdim()  { return GeneralGrid::Rdim();	};
		inline double Thetadim()  { return GeneralGrid::Thetadim();	};
		inline double Phidim() { return GeneralGrid::Phidim();	};
		inline coord<int> dim() { return GeneralGrid::dim();	};

		//Sets the dimensions for each side
		void Xdim(int in)  ;
		void Ydim(int in)  ;
		void Zdim(int in) ;
		void Rdim(int in)  ;
		void Thetadim(int in)  ;
		void Phidim(int in) ;
		void dim(const coord<int> &in) ;

		//min/max extraction
		inline double  Xmax() { return GeneralGrid::Xmax();	}
		inline double  Xmin() { return GeneralGrid::Xmin();	};
		inline double  Ymax() { return GeneralGrid::Ymax();	};
		inline double  Ymin() { return GeneralGrid::Ymin();	};
		inline double  Zmin() { return GeneralGrid::Zmin();	};
		inline double  Zmax() { return GeneralGrid::Zmax();	};
		inline double  Rmax() { return GeneralGrid::Rmax();	};
		inline double  Rmin() { return GeneralGrid::Rmin();	};
		inline double  Thetamax() { return GeneralGrid::Thetamax();	};
		inline double  Thetamin() { return GeneralGrid::Thetamin();	};
		inline double  Phimin() { return GeneralGrid::Phimin();	};
		inline double  Phimax() { return GeneralGrid::Phimax();	};
		inline coord<>  Min() { return GeneralGrid::Min();	};
		inline coord<>  Max() { return GeneralGrid::Max();	};

		//for cartisian based grids
		void Xmax(double in) ;
		void Xmin(double in) ;
		void Ymax(double in) ;
		void Ymin(double in) ;
		void Zmin(double in) ;
		void Zmax(double in) ;
		void Rmax(double in) ;
		void Rmin(double in) ;
		void Thetamax(double in) ;
		void Thetamin(double in) ;
		void Phimin(double in) ;
		void Phimax(double in) ;
		void Min(const coord<> &in) ;
		void Max(const coord<> &in) ;

		inline double x(int i)	{ return data_(i).x();	}
		inline double y(int i)	{ return data_(i).y();	}
		inline double z(int i)	{ return data_(i).z();	}

		inline double x(int i, double t)	{ return data_(i).x();	}
		inline double y(int i, double t)	{ return data_(i).y();	}
		inline double z(int i, double t)	{ return data_(i).z();	}

		coord<> operator()(int i, int j, int k);
		coord<> operator()(int i);
		coord<> operator()(int i, int j, int k,double t);
		coord<> operator()(int i, double t);
		FullCubeGrid &operator=(const FullCubeGrid &rhs);

		void Translate(double xdist, double ydist, double zdist);
		void Translate(coord<> &dist)	{ Translate(dist.x(), dist.y(), dist.z());	}

		void Center(double xC, double yC, double zC);
		void Center(coord<> &cen)	{ Center(cen.x(), cen.y(), cen.z());	}

		void Scale(double xS, double yS, double zC);
		void Scale(coord<> &sca)	{ Scale(sca.x(), sca.y(), sca.z());		}

};

class FullCubeGridIter : public GeneralGridIter{
	private:
		FullCubeGrid *myg_;

	public:

		FullCubeGridIter() :
			myg_(NULL)
		{}

		FullCubeGridIter(FullCubeGrid &in):
			GeneralGridIter(in),
			myg_(&in)
		{}

		FullCubeGridIter(FullCubeGridIter &cp):
			GeneralGridIter(cp),
			myg_(cp.myg_)
		{}

		~FullCubeGridIter()
		{
			myg_=NULL;
		}

		void SetGrid(FullCubeGrid &in)
		{
			GeneralGridIter::SetGrid(in);
			myg_=&in;
			reset();
		}

		inline coord<> &Point()
		{
			Pt_(myg_->data_(iterx).x(), myg_->data_(itery).y(), myg_->data_(iterz).z());
			return Pt_;
		}

		inline double x()	{ return myg_->data_(iterx).x();	}
		inline double y()	{ return myg_->data_(itery).y();	}
		inline double z()	{ return myg_->data_(iterz).z();	}
		inline double r()		{ return  myg_->data_(iterx).r();	}
		inline double theta()	{ return  myg_->data_(itery).theta();	}
		inline double phi()		{ return  myg_->data_(iterz).phi();	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline coord<> &Point(double t)
		{
			Pt_(myg_->data_(iterx).x(), myg_->data_(itery).y(), myg_->data_(iterz).z());
			return Pt_;
		}

		inline double x(double t)	{ return myg_->data_(iterx).x();	}
		inline double y(double t)	{ return myg_->data_(itery).y();	}
		inline double z(double t)	{ return myg_->data_(iterz).z();	}
		inline double r(double t)		{ return  myg_->data_(iterx).r();	}
		inline double theta(double t)	{ return  myg_->data_(itery).theta();	}
		inline double phi(double t)		{ return  myg_->data_(iterz).phi();	}


		FullCubeGridIter &operator=(const FullCubeGridIter &rhs)
		{
			if(this==&rhs) return *this;
			GeneralGridIter::operator=(rhs);
			myg_=rhs.myg_;
			return *this;
		}

		void operator++()
		{
			if(!notended_) myg_->Additerr();
			if(iterx<myg_->Xdim()-1){
				iterx++;
				return;
			}else{
				iterx=0;
			}
			if(itery<myg_->Ydim()-1){
				itery++;
				return;
			}else{
				itery=0;
			}
			if(iterz<myg_->Zdim()-1){
				iterz++;
				return ;
			}else{
				notended_=false;
			}
		}

		void XadvanceGrid()
		{
			xended_=false;
		}


		void YadvanceGrid()
		{
			if(itery<myg_->Ydim()-1){
				itery++;
				return;
			}else{
				yended_=false;
			}
		}

		void ZadvanceGrid()
		{
			if(iterz<myg_->Zdim()-1){
				iterz++;
				return;
			}else{
				zended_=false;
			}
		}

		void RadvanceGrid()
		{
			xended_=false;
		}


		void ThetaadvanceGrid()
		{
			if(itery<myg_->Ydim()-1){
				itery++;
				return;
			}else{
				yended_=false;
			}
		}

		void PhiadvanceGrid()
		{
			if(iterz<myg_->Zdim()-1){
				iterz++;
				return;
			}else{
				zended_=false;
			}
		}

		inline void operator++(int){	operator++();	}
};

END_BL_NAMESPACE


#endif



