

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
 	surfspheregrid.cc-->part of the grids classes: a surface sphere grid
 */


#ifndef _surfsphere_grid_cc_
#define _surfsphere_grid_cc_


#include "container/grids/surfspheregrid.h"
#include "container/grids/gengrid.h"


BEGIN_BL_NAMESPACE


void SurfSphereGrid::init()
{
	int i=0;
	double dth=(Ymax()-Ymin())/Ydim();
	double dph=(Zmax()-Zmin())/Zdim();
	double rstart=Xmin();
	double thstart=Ymin();
	double phstart=Zmin();

	for(i=0;i<Ydim();i++){
		data_(i)=coord<>(rstart, thstart, phstart);
		thstart+=dth;
		phstart+=dph;
	}
}

void SurfSphereGrid::Xdim(int in)			{ GeneralGrid::Xdim(in);	init(); }
void SurfSphereGrid::Ydim(int in)			{ GeneralGrid::Ydim(in);	init();	}
void SurfSphereGrid::Zdim(int in)			{ GeneralGrid::Zdim(in);	init();	}
void SurfSphereGrid::Rdim(int in)  			{ GeneralGrid::Rdim(in);	init(); };
void SurfSphereGrid::Thetadim(int in) 		{ GeneralGrid::Thetadim(in);	init(); } ;
void SurfSphereGrid::Phidim(int in) 		{ GeneralGrid::Phidim(in);	init(); };
void SurfSphereGrid::dim(const coord<int> &in) {GeneralGrid::dim(in);	init(); }


//for cartisian based grids
void SurfSphereGrid::Xmax(double in)			{  GeneralGrid::Xmax(in);	init(); }
void SurfSphereGrid::Xmin(double in)			{  GeneralGrid::Xmin(in);	init(); }
void SurfSphereGrid::Ymax(double in)			{  GeneralGrid::Ymax(in);	init(); }
void SurfSphereGrid::Ymin(double in)			{  GeneralGrid::Ymin(in);	init(); }
void SurfSphereGrid::Zmin(double in)		 	{  GeneralGrid::Zmin(in);   init(); }
void SurfSphereGrid::Zmax(double in)			{  GeneralGrid::Zmax(in);	init(); }
void SurfSphereGrid::Min(const coord<> &in)			{  GeneralGrid::Min(in);	init(); }
void SurfSphereGrid::Max(const coord<> &in)			{  GeneralGrid::Max(in);	init(); }
void SurfSphereGrid::Rmin(double in)			{  GeneralGrid::Rmin(in); GeneralGrid::Rmax(in);	init(); } ;
void SurfSphereGrid::Rmax(double in)			{  GeneralGrid::Rmin(in);	GeneralGrid::Rmax(in);init(); } ;
void SurfSphereGrid::Thetamax(double in) 		{  GeneralGrid::Thetamax(in);	init(); };
void SurfSphereGrid::Thetamin(double in) 		{  GeneralGrid::Thetamin(in);	init(); };
void SurfSphereGrid::Phimin(double in) 			{  GeneralGrid::Phimin(in);	init(); };
void SurfSphereGrid::Phimax(double in) 			{  GeneralGrid::Phimax(in);	init(); };



coord<> SurfSphereGrid::Point(int i, int j, int k)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());
}



coord<> SurfSphereGrid::operator()(int i, int j, int k)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());
}

coord<> SurfSphereGrid::operator()(int i)
{
	return coord<>(data_(i).r(), data_(i).theta(),data_(i).phi());
}

coord<> SurfSphereGrid::Point(int i, int j, int k, double t)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());
}



coord<> SurfSphereGrid::operator()(int i, int j, int k, double t)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());
}

coord<> SurfSphereGrid::operator()(int i, double t)
{
	return coord<>(data_(i).r(), data_(i).theta(),data_(i).phi());
}


SurfSphereGrid &SurfSphereGrid::operator=(const SurfSphereGrid &rhs)
{
	if(this==&rhs)  return *this;
	GeneralGrid::operator=(rhs);
	data_ = (rhs.data_);
	return *this;
}

void SurfSphereGrid::Scale(double rS)
{
	if(rS>=0){
		min_.r()*=rS;
		max_.r()*=rS;
	}else{
		min_.r()/=rS;
		max_.r()/=rS;
	}
	init();
}


END_BL_NAMESPACE


#endif



