

#ifndef _uniformgrid_cc_
#define _uniformgrid_cc_


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
 	unigrid.cc-->part of the grids classes: a uniform rectangular grid
 */

#include "container/grids/unigrid.h"
#include "container/grids/gengrid.h"


BEGIN_BL_NAMESPACE


void UniformGrid::init()
{
	int i=0;
	dx_=(GeneralGrid::Xmax()-GeneralGrid::Xmin())/GeneralGrid::Xdim();
	dy_=(GeneralGrid::Ymax()-GeneralGrid::Ymin())/GeneralGrid::Ydim();
	dz_=(GeneralGrid::Zmax()-GeneralGrid::Zmin())/GeneralGrid::Zdim();
	double xstart=GeneralGrid::Xmin()+dx_/2.; //center the pt in the middle of the voxel
	double ystart=GeneralGrid::Ymin()+dy_/2.;
	double zstart=GeneralGrid::Zmin()+dz_/2.;
	xdata_.resize(GeneralGrid::Xdim());
	ydata_.resize(GeneralGrid::Ydim());
	zdata_.resize(GeneralGrid::Zdim());



	for(i=0;i<GeneralGrid::Xdim();i++){
		xdata_(i)=(xstart);
		xstart+=dx_;
	}

	for(i=0;i<GeneralGrid::Ydim();i++){
		ydata_(i)=(ystart);
		ystart+=dy_;
	}

	for(i=0;i<GeneralGrid::Zdim();i++){
		zdata_(i)=(zstart);
		zstart+=dz_;
	}
	CellVol_=dx_*dy_*dz_;

	center_((GeneralGrid::Xmax()+GeneralGrid::Xmin())/2., (GeneralGrid::Ymax()+GeneralGrid::Ymin())/2.,(GeneralGrid::Zmax()+GeneralGrid::Zmin())/2.);
}

void UniformGrid::Xdim(int in)			{ GeneralGrid::Xdim(in);	init(); }
void UniformGrid::Ydim(int in)			{ GeneralGrid::Ydim(in);	init();	}
void UniformGrid::Zdim(int in)			{ GeneralGrid::Zdim(in);	init();	}
void UniformGrid::dim(const coord<int> &in) {GeneralGrid::dim(in);	init(); }


//for cartisian based grids
void UniformGrid::Xmax(double in)			{  GeneralGrid::Xmax(in);	init(); }
void UniformGrid::Xmin(double in)			{  GeneralGrid::Xmin(in);	init(); }
void UniformGrid::Ymax(double in)			{  GeneralGrid::Ymax(in);	init(); }
void UniformGrid::Ymin(double in)			{  GeneralGrid::Ymin(in);	init(); }
void UniformGrid::Zmin(double in)		 	{  GeneralGrid::Zmin(in);   init(); }
void UniformGrid::Zmax(double in)			{  GeneralGrid::Zmax(in);	init(); }
void UniformGrid::Min(const coord<> &in)			{  GeneralGrid::Min(in);	init(); }
void UniformGrid::Max(const coord<> &in)			{  GeneralGrid::Max(in);	init(); }


coord<> UniformGrid::Point(int i, int j, int k)
{
	return coord<>(xdata_(i), ydata_(j), zdata_(k));
}


coord<> UniformGrid::operator()(int i, int j, int k)
{
	return coord<>(xdata_(i), ydata_(j),zdata_(k));
}

coord<> UniformGrid::operator()(int i)
{
	return coord<>(xdata_(i), ydata_(i),zdata_(i));
}

UniformGrid &UniformGrid::operator=(const UniformGrid &rhs)
{
	if(this==&rhs)  return *this;
	GeneralGrid::operator=(rhs);
	xdata_ = (rhs.xdata_);
	ydata_ = (rhs.ydata_);
	zdata_ = (rhs.zdata_);
	CellVol_=rhs.CellVol_;
	dx_=rhs.dx_;
	dy_=rhs.dy_;
	dz_=rhs.dz_;

	return *this;
}


void UniformGrid::Translate(double distx, double disty, double distz)
{

	//double xdist=abs(center_.x()-nxc);
	//double ydist=abs(center_.y()-nyc);
	//double zdist=abs(center_.z()-nzc);
	min_.x()+=distx;
	min_.y()+=disty;
	min_.z()+=distz;

	max_.x()+=distx;
	max_.y()+=disty;
	max_.z()+=distz;
	init();

}

void UniformGrid::Center(double xC, double yC, double zC)
{
	double oldxd=abs(center_.x()-min_.x());
	double oldyd=abs(center_.y()-min_.y());
	double oldzd=abs(center_.z()-min_.z());
	min_.x()=xC-oldxd;
	min_.y()=yC-oldyd;
	min_.z()=zC-oldzd;

	max_.x()=xC+oldxd;
	max_.y()=yC+oldyd;
	max_.z()=zC+oldzd;

	init();
}

void UniformGrid::Scale(double xS, double yS, double zS)
{
	double oldxd=abs(center_.x()-min_.x());
	double oldyd=abs(center_.y()-min_.y());
	double oldzd=abs(center_.z()-min_.z());
	if(xS>=0){
		min_.x()=center_.x()-oldxd*xS;
		max_.x()=center_.x()+oldxd*xS;
	}else{
		min_.x()= center_.x()-oldxd/xS;
		max_.x()= center_.x()+oldxd/xS;
	}

	if(yS>=0){
		min_.y()=center_.y()-oldyd*yS;
		max_.y()=center_.y()+oldyd*yS;
	}else{
		min_.y()= center_.y()-oldyd/yS;
		max_.y()= center_.y()+oldyd/yS;
	}

	if(zS>=0){
		min_.z()=center_.z()-oldzd*zS;
		max_.z()=center_.z()+oldzd*zS;
	}else{
		min_.z()= center_.z()-oldzd/zS;
		max_.z()= center_.z()+oldzd/zS;
	}
	init();
}



END_BL_NAMESPACE



#endif


