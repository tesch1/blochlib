



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
 	fullspheregrid.h-->part of the grids classes: a full sphere grid
 */

#ifndef _fullspheregrid_h_
#define _fullspheregrid_h_


/* Spherical Gird Generator Class
	generates a solid sphere of points-->each dimension is the same size
	and each direction has the same maxes and mins... (i.e. a cube)

*/

#include <iostream>
#include <string>
#include "container/grids/gengrid.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


#define Pi 3.14159265358979323846

//predec for iterator (define below)

class FullSphereGridIter;


class FullSphereGrid : public GeneralGrid {
	friend class FullSphereGridIter;
	private:
		Vector<coord<> > data_;		//holds the coord data
		coord<> offset_;

	public:
		typedef FullSphereGridIter iterator;

		FullSphereGrid(): GeneralGrid() {}

		FullSphereGrid(int dims):
			GeneralGrid(dims,dims,dims), data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		FullSphereGrid(const FullSphereGrid &rhs):
			GeneralGrid(rhs), data_(rhs.data_){}

		FullSphereGrid(double minr, double maxr, int dims):
			GeneralGrid(coord<>(minr,0.,0.),//,spherical),
			coord<>(maxr, Pi, 2.*Pi),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		FullSphereGrid(double minr, double maxr, double minth, double maxth, int dims):
			GeneralGrid(coord<>(minr,minth,0.),//,spherical),
			coord<>(maxr, maxth, 2.*Pi),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}


		FullSphereGrid(double minr, double maxr, double minth, double maxth,double minph, double maxph, int dims):
			GeneralGrid(coord<>(minr,minth,minph),//,spherical),
			coord<>(maxr, maxth, maxph),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		FullSphereGrid(coord<> &mins, coord<> &maxs, coord<int,3> dims):
			GeneralGrid(mins, maxs,coord<int>(dims[0])),
			data_(dims[0],coord<>(0))//,spherical))
		{
			init();
		}


		~FullSphereGrid(){}

		void init();
		coord<> Point(int i, int j, int k) ;

		inline double r()		{ return data_(iterx).r();	}
		inline double theta()	{ return data_(itery).theta();	}
		inline double phi()		{ return data_(iterz).phi();	}


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

		coord<> operator()(int i, int j, int k);
		coord<> operator()(int i);
		FullSphereGrid &operator=(const FullSphereGrid &rhs);

		//void Translate(double xdist, double ydist, double zdist);
		//void Translate(coord<> &dist)	{ Translate(dist.x(), dist.y(), dist.z());	}

		//void Center(double xC, double yC, double zC);
		//void Center(coord<> &cen)	{ Center(cen.x(), cen.y(), cen.z());	}

		//only 'R' is valid (or makes any sence)
		//the two overloads oare for when you want to switch the Engines
		//and do not want to rewrite the entire world
		void Scale(coord<> &Rs){	Scale(Rs[0]);	}
		void Scale(double rs, double ths, double phs){	Scale(rs);	}
		void Scale(double rS);


};


class FullSphereGridIter : public GeneralGridIter{
	private:
		FullSphereGrid *myg_;

	public:

		FullSphereGridIter() :
			myg_(NULL)
		{}

		FullSphereGridIter(FullSphereGrid &in):
			GeneralGridIter(in),
			myg_(&in)
		{}

		FullSphereGridIter(FullSphereGridIter &cp):
			GeneralGridIter(cp),
			myg_(cp.myg_)
		{}

		~FullSphereGridIter()
		{
			myg_=NULL;
		}

		void SetGrid(FullSphereGrid &in)
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


		FullSphereGridIter &operator=(const FullSphereGridIter &rhs)
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
				return;
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



