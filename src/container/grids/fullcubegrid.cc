

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
 	fullcubegrid.cc-->part of the grids classes a full cube grid
 */

#ifndef _fullcube_grid_cc_
#define _fullcude_grid_cc_


#include "fullcubegrid.h"
#include "gengrid.h"


BEGIN_BL_NAMESPACE


void FullCubeGrid::init()
{
	int i=0;
	double dx=(Xmax()-Xmin())/Xdim();
	double dy=(Ymax()-Ymin())/Ydim();
	double dz=(Zmax()-Zmin())/Zdim();
	double xstart=Xmin()+dx/2.; //center the pt in the middle of the voxel
	double ystart=Ymin()+dy/2.;
	double zstart=Zmin()+dz/2.;

	for(i=0;i<Xdim();i++){
		data_(i)(xstart, ystart, zstart);
		xstart+=dx;
		ystart+=dy;
		zstart+=dz;
	}
	center_((Xmax()+Xmin())/2., (Ymax()+Ymin())/2.,(Zmax()+Zmin())/2.);
}



coord<> FullCubeGrid::Point(int i, int j, int k)
{
	return coord<>(data_(i).x(), data_(j).y(),data_(k).z());
}


coord<> FullCubeGrid::Point(int i, int j, int k, double t)
{
	return coord<>(data_(i).x(), data_(j).y(),data_(k).z());
}

void FullCubeGrid::Xdim(int in)			{ GeneralGrid::Xdim(in);	init(); }
void FullCubeGrid::Ydim(int in)			{ GeneralGrid::Ydim(in);	init();	}
void FullCubeGrid::Zdim(int in)			{ GeneralGrid::Zdim(in);	init();	}
void FullCubeGrid::dim(const coord<int> &in) {GeneralGrid::dim(in);	init(); }
void FullCubeGrid::Rdim(int in) 		{ GeneralGrid::Rdim(in);	init(); } ;
void FullCubeGrid::Thetadim(int in) 		{ GeneralGrid::Thetadim(in);	init(); } ;
void FullCubeGrid::Phidim(int in) 		{ GeneralGrid::Phidim(in);	init(); };


//for cartisian based grids
void FullCubeGrid::Xmax(double in)			{  GeneralGrid::Xmax(in);	init(); }
void FullCubeGrid::Xmin(double in)			{  GeneralGrid::Xmin(in);	init(); }
void FullCubeGrid::Ymax(double in)			{  GeneralGrid::Ymax(in);	init(); }
void FullCubeGrid::Ymin(double in)			{  GeneralGrid::Ymin(in);	init(); }
void FullCubeGrid::Zmin(double in)		 	{  GeneralGrid::Zmin(in);   init(); }
void FullCubeGrid::Zmax(double in)			{  GeneralGrid::Zmax(in);	init(); }
void FullCubeGrid::Min(const coord<> &in)			{  GeneralGrid::Min(in);	init(); }
void FullCubeGrid::Max(const coord<> &in)			{  GeneralGrid::Max(in);	init(); }
void FullCubeGrid::Rmax(double in) 			{  GeneralGrid::Rmax(in);	init(); };
void FullCubeGrid::Rmin(double in)			{  GeneralGrid::Rmin(in);	init(); } ;
void FullCubeGrid::Thetamax(double in) 		{  GeneralGrid::Thetamax(in);	init(); };
void FullCubeGrid::Thetamin(double in) 		{  GeneralGrid::Thetamin(in);	init(); };
void FullCubeGrid::Phimin(double in) 			{  GeneralGrid::Phimin(in);	init(); };
void FullCubeGrid::Phimax(double in) 			{  GeneralGrid::Phimax(in);	init(); };

coord<> FullCubeGrid::operator()(int i, int j, int k)
{
	return coord<>(data_(i).x(), data_(j).y(),data_(k).z());
}

coord<> FullCubeGrid::operator()(int i)
{
	return coord<>(data_(i).x(), data_(i).y(),data_(i).z());
}

coord<> FullCubeGrid::operator()(int i, int j, int k, double t)
{
	return coord<>(data_(i).x(), data_(j).y(),data_(k).z());
}

coord<> FullCubeGrid::operator()(int i, double t)
{
	return coord<>(data_(i).x(), data_(i).y(),data_(i).z());
}


FullCubeGrid &FullCubeGrid::operator=(const FullCubeGrid &rhs)
{
	if(this==&rhs)  return *this;
	GeneralGrid::operator=(rhs);
	data_ = (rhs.data_);
	return *this;
}


void FullCubeGrid::Translate(double distx, double disty, double distz)
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

void FullCubeGrid::Center(double xC, double yC, double zC)
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

void FullCubeGrid::Scale(double xS, double yS, double zS)
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


