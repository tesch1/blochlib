


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
 	fullspheregrid.cc-->part of the grids classes: a full sphere grid
 */

#ifndef _fullsphere_grid_cc_
#define _fullsphere_grid_cc_


#include "container/grids/fullspheregrid.h"
#include "container/grids/gengrid.h"


BEGIN_BL_NAMESPACE


void FullSphereGrid::init()
{
	int i=0;
	double drr=(Xmax()-Xmin())/Xdim();
	double dth=(Ymax()-Ymin())/Ydim();
	double dph=(Zmax()-Zmin())/Zdim();
	double rstart=Xmin();
	double thstart=Ymin();
	double phstart=Zmin();

	for(i=0;i<Xdim();i++){
		data_(i)=coord<>(rstart, thstart, phstart);//, spherical);
		rstart+=drr;
		thstart+=dth;
		phstart+=dph;
	}
}




coord<> FullSphereGrid::Point(int i, int j, int k)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());//, spherical);
}


coord<> FullSphereGrid::operator()(int i, int j, int k)
{
	return coord<>(data_(i).r(), data_(j).theta(),data_(k).phi());//, spherical);
}

coord<> FullSphereGrid::operator()(int i)
{
	return coord<>(data_(i).r(), data_(i).theta(),data_(i).phi());//, spherical);
}


FullSphereGrid &FullSphereGrid::operator=(const FullSphereGrid &rhs)
{
	if(this==&rhs)  return *this;
	GeneralGrid::operator=(rhs);
	data_ = (rhs.data_);
	return *this;
}

void FullSphereGrid::Xdim(int in)			{ GeneralGrid::Xdim(in);	init(); }
void FullSphereGrid::Ydim(int in)			{ GeneralGrid::Ydim(in);	init();	}
void FullSphereGrid::Zdim(int in)			{ GeneralGrid::Zdim(in);	init();	}
void FullSphereGrid::dim(const coord<int> &in) {GeneralGrid::dim(in);	init(); }
void FullSphereGrid::Rdim(int in) 		{ GeneralGrid::Rdim(in);	init(); } ;
void FullSphereGrid::Thetadim(int in) 		{ GeneralGrid::Thetadim(in);	init(); } ;
void FullSphereGrid::Phidim(int in) 		{ GeneralGrid::Phidim(in);	init(); };


//for cartisian based grids
void FullSphereGrid::Xmax(double in)			{  GeneralGrid::Xmax(in);	init(); }
void FullSphereGrid::Xmin(double in)			{  GeneralGrid::Xmin(in);	init(); }
void FullSphereGrid::Ymax(double in)			{  GeneralGrid::Ymax(in);	init(); }
void FullSphereGrid::Ymin(double in)			{  GeneralGrid::Ymin(in);	init(); }
void FullSphereGrid::Zmin(double in)		 	{  GeneralGrid::Zmin(in);   init(); }
void FullSphereGrid::Zmax(double in)			{  GeneralGrid::Zmax(in);	init(); }
void FullSphereGrid::Min(const coord<> &in)			{  GeneralGrid::Min(in);	init(); }
void FullSphereGrid::Max(const coord<> &in)			{  GeneralGrid::Max(in);	init(); }
void FullSphereGrid::Rmax(double in) 			{  GeneralGrid::Rmax(in);	init(); };
void FullSphereGrid::Rmin(double in)			{  GeneralGrid::Rmin(in);	init(); } ;
void FullSphereGrid::Thetamax(double in) 		{  GeneralGrid::Thetamax(in);	init(); };
void FullSphereGrid::Thetamin(double in) 		{  GeneralGrid::Thetamin(in);	init(); };
void FullSphereGrid::Phimin(double in) 			{  GeneralGrid::Phimin(in);	init(); };
void FullSphereGrid::Phimax(double in) 			{  GeneralGrid::Phimax(in);	init(); };


/*void FullSphereGrid::Translate(double xdist, double ydist, double zdist)
{
	coord<> tm;
	for(int i=0;i<data_.size();i++){
		tm=data_(i).sph2cart();
		tm.x()+=xdist;
		tm.y()+=ydist;
		tm.z()+=zdist;
		data_(i)=tm.cart2sph();
	}
}*/

/*void FullSphereGrid::Center(double xC, double yC, double zC)
{
	center_(xC, yC, zC);
}
*/
void FullSphereGrid::Scale(double rS)
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



