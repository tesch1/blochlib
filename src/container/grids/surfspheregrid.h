

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
 	surfspheregrid.h-->part of the grids classes: a surface sphere grid
 */


#ifndef _surfspheregrid_h_
#define _surfspheregrid_h_


/* Surface Grid Generator Generator Class
	generates a sphere of points-->not filled

*/

#include <iostream>
#include <string>
#include "container/grids/gengrid.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"

BEGIN_BL_NAMESPACE


#define Pi 3.14159265358979323846

//predec for Iterator
class SurfSphereGridIter;


class SurfSphereGrid : public GeneralGrid {
	friend class SurfSphereGridIter;
	private:
		Vector<coord<> > data_;		//holds the coord data
		coord<> outdat_;

	public:

		typedef SurfSphereGridIter iterator;

		SurfSphereGrid(): GeneralGrid() {}

		SurfSphereGrid(int dims):
			GeneralGrid(dims,dims,dims),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		SurfSphereGrid(const SurfSphereGrid &rhs):
			GeneralGrid(rhs), data_(rhs.data_){}

		SurfSphereGrid(double maxr, int dims):
			GeneralGrid(coord<>(maxr,0.,0.),//,spherical),
			coord<>(maxr, Pi, 2.),//*Pi,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		SurfSphereGrid(double maxth,double maxphi, int dims):
			GeneralGrid(coord<>(1.0,0.,0.),//,spherical),
			coord<>(1.0, maxth, maxphi),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		SurfSphereGrid(double maxr, double minth, double maxth, int dims):
			GeneralGrid(coord<>(maxr,minth,0.),//,spherical),
			coord<>(maxr, maxth, 2.*Pi),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}


		SurfSphereGrid(double maxr, double minth, double maxth,double minph, double maxph, int dims):
			GeneralGrid(coord<>(maxr,minth,minph),//,spherical),
			coord<>(maxr, maxth, maxph),//,spherical),
			coord<int>(dims)),
			data_(dims,coord<>(0))//,spherical))
		{
			init();
		}

		SurfSphereGrid(coord<> &mins, coord<> &maxs, coord<int,3> &dims):
			GeneralGrid(coord<>(maxs(0), mins(1), mins(2)),maxs,dims),
			data_(dims(0),coord<>(0))//,spherical))
		{
			init();
		}

		~SurfSphereGrid(){}

		void init();
		coord<> Point(int i, int j, int k) ;

		inline double x(int i)	{ return data_(i).x();	}
		inline double y(int i)	{ return data_(i).y();	}
		inline double z(int i)	{ return data_(i).z();	}


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


		void operator++();
		coord<> operator()(int i, int j, int k);
		coord<> operator()(int i);

		SurfSphereGrid &operator=(const SurfSphereGrid &rhs);

		//only 'R' is valid (or makes any sence)
		//the two overloads are for when you want to switch the Engines
		//and do not want to rewrite the entire world
		void Scale(coord<> &Rs){	Scale(Rs[0]);	}
		void Scale(double rs, double ths, double phs){	Scale(rs);	}
		void Scale(double rS);

	//for 'teim dependat' grids alothough this is NOT time dep
	// function added for consistancy
		coord<> Point(int i, int j, int k, double t) ;
		inline double x(int i, double t)	{ return data_(i).x();	}
		inline double y(int i, double t)	{ return data_(i).y();	}
		inline double z(int i, double t)	{ return data_(i).z();	}
		coord<> operator()(int i, int j, int k, double t);
		coord<> operator()(int i, double t);


};

class SurfSphereGridIter : public GeneralGridIter{
	private:
		SurfSphereGrid *myg_;

	public:

		SurfSphereGridIter() :
			myg_(NULL)
		{}

		SurfSphereGridIter(SurfSphereGrid &in):
			GeneralGridIter(in),
			myg_(&in)

		{}

		SurfSphereGridIter(SurfSphereGridIter &cp):
			GeneralGridIter(cp),
			myg_(cp.myg_)

		{}

		~SurfSphereGridIter()
		{
			myg_=NULL;
		}

		void SetGrid(SurfSphereGrid &in)
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



		SurfSphereGridIter &operator=(const SurfSphereGridIter &rhs)
		{
			if(this==&rhs) return *this;
			GeneralGridIter::operator=(rhs);
			myg_=rhs.myg_;
			return *this;
		}

		void operator++()
		{
			if(!notended_) myg_->Additerr();

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
