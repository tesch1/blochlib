
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
 	unigrid.h-->part of the grids classes: a uniform rectangular grid
 */

#ifndef _unigrid_h_
#define _unigrid_h_


/* Uniform Gird Generator Class
	generates a 3D rectangular grid of points in a uniform
	(equal step sizes) manner
*/

#include <iostream>
#include <string>
#include "container/grids/gengrid.h"

BEGIN_BL_NAMESPACE


//predeclaration (defined below this class)
class UniforGridIter;


class UniformGrid : public GeneralGrid {
	friend class UniformGridIter;
	private:
		Vector<double> xdata_;		//holds the x coord data
		Vector<double> ydata_;		//holds the y coord data
		Vector<double> zdata_;		//holds the z coord data

		double CellVol_;			//for a uniform grid the cell volume is always the same....
		double dx_, dy_, dz_;		//cell dimensions.....

	public:
		typedef UniformGridIter iterator;

		UniformGrid(): GeneralGrid() {}

		UniformGrid(const UniformGrid &cp): GeneralGrid(cp),
			xdata_(cp.xdata_), ydata_(cp.ydata_), zdata_(cp.zdata_),
			CellVol_(cp.CellVol_),
			dx_(cp.dx_), dy_(cp.dy_), dz_(cp.dz_)

		{}

		UniformGrid(int xst, int yst, int zst):
			GeneralGrid(xst,yst,zst),
			xdata_(xst), ydata_(yst), zdata_(zst)
		{
			init();
		}

		UniformGrid(const coord<> &mins, const coord<> &maxs):
			GeneralGrid(mins,maxs),
			xdata_(1), ydata_(1), zdata_(1)
		{
			init();
		}

		UniformGrid(const coord<> &mins, const coord<> &maxs, const coord<int> &dims):
			GeneralGrid(mins,maxs,dims),
			xdata_(dims.x()), ydata_(dims.y()), zdata_(dims.z())
		{
			init();
		}

		UniformGrid(double xmax ,double ymax, double zmax,int dims):
				GeneralGrid(0,xmax, 0,ymax,0,zmax, dims,dims,dims),
				xdata_(dims), ydata_(dims), zdata_(dims)
		{
			init();
		}

		UniformGrid(double xmax ,double ymax, double zmax,
					int xsize, int ysize, int zsize):
				GeneralGrid(0,xmax, 0,ymax,0,zmax, xsize,ysize,zsize),
				xdata_(xsize), ydata_(ysize), zdata_(zsize)
		{
			init();
		}

		UniformGrid(double xmin,double xmax ,double ymin,
					double ymax,double zmin, double zmax,
					int xsize, int ysize, int zsize):
				GeneralGrid(xmin,xmax, ymin,ymax,zmin,zmax, xsize,ysize,zsize)
		{
			init();
		}


		void init();
		coord<> Point(int i, int j, int k) ;

		inline double x()	{ return xdata_(iterx);	}
		inline double y()	{ return ydata_(itery);	}
		inline double z()	{ return zdata_(iterz);	}

		inline double dx()	{ return dx_;	}
		inline double dy()	{ return dy_;	}
		inline double dz()	{ return dz_;	}
		inline coord<> dr()	{ return coord<>(dx_, dy_, dz_);	}

		inline double d(int dir, int which)	{
			switch(dir){
				case 0: return dx_;
				case 1: return dy_;
				case 2: return dz_;
				default: return 0;
			}
		}

		inline double dx(int i)	{ return dx_;	}
		inline double dy(int i)	{ return dy_;	}
		inline double dz(int i)	{ return dz_;	}

		inline coord<> dr(int i)	{ return coord<>(dx_, dy_, dz_);	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline double dx(int i, double t)	{ return dx_;	}
		inline double dy(int i, double t)	{ return dy_;	}
		inline double dz(int i, double t)	{ return dz_;	}

		inline double x(int i)	{ return xdata_(i);	}
		inline double y(int i)	{ return ydata_(i);	}
		inline double z(int i)	{ return zdata_(i);	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline double x(int i, double t)	{ return xdata_(i);	}
		inline double y(int i, double t)	{ return ydata_(i);	}
		inline double z(int i, double t)	{ return zdata_(i);	}

		//dimensions for each side
		inline int Xdim(){	return  GeneralGrid::Xdim();	}
		inline int Ydim() {	return  GeneralGrid::Ydim();	}
		inline int Zdim() {	return  GeneralGrid::Zdim();	}
		inline coord<int> dim() {	return  GeneralGrid::dim();	}

		inline coord<> Min(){	return GeneralGrid::Min();	}
		inline coord<> Max(){	return GeneralGrid::Max();	}

		inline coord<> min(){	return GeneralGrid::Min();	}
		inline coord<> max(){	return GeneralGrid::Max();	}

		//dimensions for each side
		void Xdim(int in)  ;
		void Ydim(int in)  ;
		void Zdim(int in) ;
		void dim(const coord<int> &in) ;

		//for cartisian based grids
		void Xmax(double in) ;
		void Xmin(double in) ;
		void Ymax(double in) ;
		void Ymin(double in) ;
		void Zmin(double in) ;
		void Zmax(double in) ;
		void Min(const coord<> &in) ;
		void Max(const coord<> &in) ;

		inline double CellVolume(){	return CellVol_;	}
		inline double CellVolume(int i){	return CellVol_;	}
		inline double CellVolume(int i, double t){	return CellVol_;	}

		inline double cellVolume(){	return CellVol_;	}
		inline double cellVolume(int i){	return CellVol_;	}
		inline double cellVolume(int i, double t){	return CellVol_;	}


		int size() const		{	return xdata_.size()*ydata_.size()*zdata_.size();	}


		coord<> operator()(int i, int j, int k);
		coord<> operator()(int i);
		UniformGrid &operator=(const UniformGrid &rhs);

		void Translate(double xdist, double ydist, double zdist);
		void Translate(coord<> &dist)	{ Translate(dist.x(), dist.y(), dist.z());	}

		void Center(double xC, double yC, double zC);
		void Center(coord<> &cen)	{ Center(cen.x(), cen.y(), cen.z());	}

		void Scale(double xS, double yS, double zC);
		void Scale(coord<> &sca)	{ Scale(sca.x(), sca.y(), sca.z());		}


};


class UniformGridIter : public GeneralGridIter{
	private:
		UniformGrid *myg_;

	public:

		UniformGridIter() :
			myg_(NULL)
		{}

		UniformGridIter(UniformGrid &in):
			GeneralGridIter(in),
			myg_(&in)
		{}

		UniformGridIter(UniformGridIter &cp):
			GeneralGridIter(cp),
			myg_(cp.myg_)
		{}

		~UniformGridIter()
		{	myg_=NULL;		}

		void SetGrid(UniformGrid &in)
		{
			GeneralGridIter::SetGrid(in);
			if(myg_){	myg_=&in;	}
			else{
				myg_=new UniformGrid;
				myg_=&in;
			}
			reset();
		}

		inline coord<> &Point()
		{
			Pt_(myg_->xdata_(iterx), myg_->ydata_(itery), myg_->zdata_(iterz));
			return Pt_;
		}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline coord<> &Point(double t)
		{
			Pt_(myg_->xdata_(iterx), myg_->ydata_(itery), myg_->zdata_(iterz));
			return Pt_;
		}

		inline double x()	{ return myg_->xdata_(iterx);	}
		inline double y()	{ return myg_->ydata_(itery);	}
		inline double z()	{ return myg_->zdata_(iterz);	}
		inline double r()	{ return myg_->xdata_(iterx);	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline double x(double t)	{ return myg_->xdata_(iterx);	}
		inline double y(double t)	{ return myg_->ydata_(itery);	}
		inline double z(double t)	{ return myg_->zdata_(iterz);	}
		inline double r(double t)	{ return myg_->xdata_(iterx);	}

		inline double CellVolume(){	return myg_->CellVolume();	}
		inline double cellVolume(){	return myg_->CellVolume();	}
		inline double dx(){	return myg_->dx();	}
		inline double dy(){	return myg_->dy();	}
		inline double dz(){	return myg_->dz();	}
		inline coord<> dr(){	return myg_->dr();	}

	//for grids that depend on time...this one does not
	// but the function is added for consistancy
		inline double CellVolume(double t){	return myg_->CellVolume();	}
		inline double cellVolume(double t){	return myg_->CellVolume();	}
		inline double dx(double t){	return myg_->dx();	}
		inline double dy(double t){	return myg_->dy();	}
		inline double dz(double t){	return myg_->dz();	}
		inline coord<> dr(double t){	return myg_->dr();	}


		inline double theta()	{ return myg_->ydata_(itery);	}
		inline double phi()	{ return myg_->zdata_(iterz);	}

//for grids that depend on time...this one does not
// but the function is added for consistancy
		inline double theta(double t)	{ return myg_->ydata_(itery);	}
		inline double phi(double t)	{ return myg_->zdata_(iterz);	}

		UniformGridIter &operator=(const UniformGridIter &rhs)
		{
			if(this==&rhs) return *this;
			GeneralGridIter::operator=(rhs);
			if(myg_){		myg_=rhs.myg_;	}
			else{
				myg_=new UniformGrid;
				myg_=rhs.myg_;
			}
			return *this;
		}

		inline int curpos(){	return curpos_;	}

		inline void operator++()
		{
			if(!notended_) myg_->Additerr();
			if(iterx<myg_->Xdim()-1){
				iterx++; curpos_++;
				return;
			}else{
				iterx=0;
			}
			if(itery<myg_->Ydim()-1){
				itery++; curpos_++;
				return;
			}else{
				itery=0;
			}
			if(iterz<myg_->Zdim()-1){
				iterz++; curpos_++;
				return ;
			}else{
				notended_=false;
			}
		}

		inline void XadvanceGrid()
		{
			if(iterx<myg_->Xdim()-1){
				iterx++;
				return;
			}else{
				xended_=false;
			}
		}


		inline void YadvanceGrid()
		{
			if(itery<myg_->Ydim()-1){
				itery++;
				return;
			}else{
				yended_=false;
			}
		}

		inline void ZadvanceGrid()
		{
			if(iterz<myg_->Zdim()-1){
				iterz++;
				return;
			}else{
				zended_=false;
			}
		}

		inline void operator++(int){ 	operator++();	}
};


END_BL_NAMESPACE


#endif


