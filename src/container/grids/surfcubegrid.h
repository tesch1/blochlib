

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
 	surfcubegrid.h-->part of the grids classes: a surface cube grid
 */

#ifndef _surfcubegrid_h_
#define _surfcubegrid_h_


/* Surface Cube Generator Class
	generates a cube of points-->so it forms a box (not filled)

*/

#include <iostream>
#include <string>
#include "container/grids/gengrid.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


class SurfCubeGridIter;

class SurfCubeGrid : public GeneralGrid {
	friend class SurfCubeGridIter;
	private:
		Vector<double> xdata_;		//holds the x coord data
		Vector<double> ydata_;		//holds the y coord data
		Vector<double> zdata_;		//holds the z coord data


	public:
		typedef SurfCubeGridIter iterator;

		SurfCubeGrid(): GeneralGrid() {}

		SurfCubeGrid(int dims):
			GeneralGrid(dims,dims,dims),
			xdata_(dims,0), ydata_(dims,0), zdata_(dims,0)
		{
			init();
		}

		SurfCubeGrid(const SurfCubeGrid &rhs):
			GeneralGrid(rhs),
			xdata_(rhs.xdata_), ydata_(rhs.ydata_), zdata_(rhs.zdata_)
		{}

		SurfCubeGrid(double min, double max, int dims):
			GeneralGrid(coord<>(min),coord<>(max),coord<int>(dims)),
			xdata_(dims,0), ydata_(dims,0), zdata_(dims,0)
		{
			init();
		}

		SurfCubeGrid(double min, double max, coord<int> dims):
			GeneralGrid(coord<>(min),coord<>(max),dims),
			xdata_(dims[0],0), ydata_(dims[1],0), zdata_(dims[2],0)
		{
			init();
		}

		~SurfCubeGrid(){}

		void init();
		coord<> Point(int i, int j, int k) ;
		inline int size()	{	return (Xdim()*Ydim()+Ydim()*Zdim()+Xdim()*Zdim());	}


		inline double x(int i)	{ return xdata_(i);	}
		inline double y(int i)	{ return ydata_(i);	}
		inline double z(int i)	{ return zdata_(i);	}

		//dimension extraction
		inline int Xdim()  { return GeneralGrid::Xdim();	}
		inline int Ydim()  { return GeneralGrid::Ydim();	};
		inline int Zdim() { return GeneralGrid::Zdim();	};
		inline int Rdim()  { return GeneralGrid::Rdim();	};
		inline int Thetadim()  { return GeneralGrid::Thetadim();	};
		inline int Phidim() { return GeneralGrid::Phidim();	};
		inline coord<int> dim() { return GeneralGrid::dim();	};

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

		//Set the dimesions
		void Xdim(int in)  ;
		void Ydim(int in)  ;
		void Zdim(int in) ;
		void Rdim(int in)  ;
		void Thetadim(int in)  ;
		void Phidim(int in) ;
		void dim(const coord<int> &in) ;

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

		coord<> operator()(int i, int j, int k);
		coord<> operator()(int i);
		SurfCubeGrid &operator=(const SurfCubeGrid &rhs);

		inline void reset(){
			GeneralGrid::reset();
		}

		void Translate(double xdist, double ydist, double zdist);
		void Translate(coord<> &dist)	{ Translate(dist.x(), dist.y(), dist.z());	}

		void Center(double xC, double yC, double zC);
		void Center(coord<> &cen)	{ Center(cen.x(), cen.y(), cen.z());	}

		void Scale(double xS, double yS, double zC);
		void Scale(coord<> &sca)	{ Scale(sca.x(), sca.y(), sca.z());		}

	//for 'teim dependat' grids alothough this is NOT time dep
	// function added for consistancy
		coord<> Point(int i, int j, int k, double t) ;
		inline double x(int i, double t)	{ return xdata_(i);	}
		inline double y(int i, double t)	{ return ydata_(i);	}
		inline double z(int i, double t)	{ return zdata_(i);	}
		coord<> operator()(int i, int j, int k, double t);
		coord<> operator()(int i, double t);


};

class SurfCubeGridIter : public GeneralGridIter{
	private:
		SurfCubeGrid *myg_;
		bool holdx, holdy, holdz;

	public:

		SurfCubeGridIter() :
			myg_(NULL),
			holdx(false), holdy(false), holdz(false)
		{}

		SurfCubeGridIter(SurfCubeGrid &in):
			GeneralGridIter(in),
			myg_(&in),
			holdx(false), holdy(false), holdz(false)

		{}

		SurfCubeGridIter(SurfCubeGridIter &cp):
			GeneralGridIter(cp),
			myg_(cp.myg_),
			holdx(false), holdy(false), holdz(false)

		{}

		~SurfCubeGridIter()
		{
			myg_=NULL;
		}

		void SetGrid(SurfCubeGrid &in)
		{
			GeneralGridIter::SetGrid(in);
			myg_=&in;
			reset();
		}

		void reset()
		{
			GeneralGridIter::reset();
			holdx=false; holdy=false; holdz=false;
		}

		inline coord<> &Point()
		{
			Pt_(myg_->xdata_(iterx), myg_->ydata_(itery), myg_->zdata_(iterz));
			return Pt_;
		}

		inline double x()	{ return myg_->xdata_(iterx);	}
		inline double y()	{ return myg_->ydata_(itery);	}
		inline double z()	{ return myg_->zdata_(iterz);	}
		inline double r()		{ return  myg_->xdata_(iterx);	}
		inline double theta()	{ return  myg_->ydata_(itery);	}
		inline double phi()		{ return  myg_->zdata_(iterz);	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline coord<> &Point(double t)
		{
			Pt_(myg_->xdata_(iterx), myg_->ydata_(itery), myg_->zdata_(iterz));
			return Pt_;
		}

		inline double x(double t)	{ return myg_->xdata_(iterx);	}
		inline double y(double t)	{ return myg_->ydata_(itery);	}
		inline double z(double t)	{ return myg_->zdata_(iterz);	}
		inline double r(double t)		{ return  myg_->xdata_(iterx);	}
		inline double theta(double t)	{ return  myg_->ydata_(itery);	}
		inline double phi(double t)		{ return  myg_->zdata_(iterz);	}


		SurfCubeGridIter &operator=(const SurfCubeGridIter &rhs)
		{
			if(this==&rhs) return *this;
			GeneralGridIter::operator=(rhs);
			myg_=rhs.myg_;
			holdx=rhs.holdx;
			holdy=rhs.holdy;
			holdz=rhs.holdz;
			return *this;
		}

		void operator++()
		{
			if(!notended_) myg_->Additerr();
			if(!holdx){
				if(itery<myg_->Ydim()-1 ){
					itery++;
					return;
				}else{
					itery=0;
				}
				if(iterz<myg_->Zdim()-1 ){
					iterz++;
					return;
				}else{
					iterz=0;
					if(iterx==0) iterx=myg_->Xdim()-1;
					else{ holdx=true; iterx=0; itery=0;	}
				}
			}else if(!holdy){
				if(iterx<myg_->Xdim()-1 ){
					iterx++;
					return;
				}else{
					iterx=0;
				}
				if(iterz<myg_->Zdim()-1 ){
					iterz++;
					return;
				}else{
					iterz=0;
					if(itery==0) itery=myg_->Ydim()-1;
					else{ holdy=true; itery=0; iterx=0; }
				}
			}else if(!holdz){
				if(iterx<myg_->Xdim()-1 ){
					iterx++;
					return;
				}else{
					iterx=0;
				}
				if(itery<myg_->Ydim()-1 ){
					itery++;
					return;
				}else{
					itery=0;
					if(iterz==0) iterz=myg_->Zdim()-1;
					else notended_=false;
				}
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


