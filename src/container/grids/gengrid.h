
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
 	gengrid.h-->base class for all grid types
 */

#ifndef _gengrid_h_
#define _gengrid_h_


/* General Gird Generator Class
	the base class for all grid generators
*/

#include <iostream>
#include <string>
#include "container/grids/coords.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


class GeneralGridIter;

class GeneralGrid {
	protected:
		coord<> min_, max_; 	//x,y,z max and mins
		int dimx;
		int dimy;
		int dimz;//number of elements of each directions

		int iterx;	//itterator steps
		int itery;
		int iterz;

		bool notended_;			//have we looped through the entire list?
		bool xended_;			//looped through X
		bool yended_;			//looped through X
		bool zended_;			//looped through X
		coord<> center_;		//center point

	public:

		typedef GeneralGridIter iterator;

		const void Additerr() 
		{
			BLEXCEPTION(" itterated too far...")
		}

		const void Fileerr()
		{
			std::cerr<<std::endl<<"Error: GeneralGrid::write(ofstream)"<<std::endl;
			std::cerr<<" File cannot be written to"<<std::endl;
		}


		GeneralGrid():
			min_(0), max_(0), dimx(1),dimy(1),dimz(1),
			iterx(0), itery(0), iterz(0),notended_(true),
			xended_(true),yended_(true),zended_(true),
			center_(0){}

		GeneralGrid(int xst, int yst, int zst):
			min_(0), max_(1), dimx(xst),dimy(yst),dimz(zst),
			iterx(0), itery(0), iterz(0), notended_(true),
			xended_(true),yended_(true),zended_(true),
			center_(0) {}

		GeneralGrid(const coord<> &mins, const coord<> &maxs):
			min_(mins), max_(maxs), dimx(1),dimy(1),dimz(1),
			iterx(0), itery(0), iterz(0), notended_(true),
			center_(0) {}

		GeneralGrid(const coord<> &mins, const coord<> &maxs, const coord<int> &dims):
			min_(mins), max_(maxs), dimx(dims(0)),
			dimy(dims(1)),dimz(dims(2)),
			iterx(0), itery(0), iterz(0), notended_(true),
			xended_(true),yended_(true),zended_(true),
			center_(0) {}

		GeneralGrid(double xmin,double xmax ,double ymin,
					double ymax,double zmin, double zmax,
					int xsize, int ysize, int zsize):
			min_(xmin,ymin,zmin), max_(xmax, ymax,zmax),
			dimx(xsize),dimy(ysize), dimz(zsize),
			iterx(0), itery(0), iterz(0),  notended_(true),
			xended_(true),yended_(true),zended_(true),
			center_(0) {}

		GeneralGrid(const GeneralGrid &rhs):
			min_(rhs.min_), max_(rhs.max_),
			dimx(rhs.dimx),dimy(rhs.dimy), dimz(rhs.dimz),
			iterx(rhs.iterx), itery(rhs.itery), iterz(rhs.iterz),
			notended_(rhs.notended_),
			xended_(true),yended_(true),zended_(true),
			center_(0) {}

		GeneralGrid &operator=(const GeneralGrid &rhs){
			if(&rhs==this) return *this;
			min_=rhs.min_;
			max_=rhs.max_;
			dimx=rhs.dimx;
			dimy=rhs.dimy;
			dimz=rhs.dimz;
			iterx=rhs.iterx;
			itery=rhs.itery;
			iterz=rhs.iterz;
			notended_=rhs.notended_;
			center_=rhs.center_;
			return *this;
		}

		~GeneralGrid(){}

		inline operator bool()	{ return notended_;	}
		inline bool GridEnd()		{ return notended_;	}
		inline bool XEnd()			{ return xended_;	}
		inline bool YEnd()			{ return yended_;	}
		inline bool ZEnd()			{ return zended_;	}

		inline int Xpos()			{ return iterx;		}
		inline int Ypos()			{ return itery;		}
		inline int Zpos()			{ return iterz;		}
		inline int Pos()			{ return (iterz*dimx*dimy)+(itery*dimy)+iterx;		}


		inline coord<> MaxLimit()	const	{ return max_;	}
		inline coord<> Max()	const		{ return max_;	}
		inline void Max(const coord<> &in)		{ max_=in;	}

		inline coord<> MinLimit()	const	{ return min_;	}
		inline coord<> Min()	const		{ return min_;	}
		inline void Min(const coord<> &in)		{ min_=in;	}

		inline coord<int> Dimension()		{ return coord<int>(dimx, dimy, dimz);	}
		inline coord<int> dim()				{ return coord<int>(dimx, dimy, dimz);	}
		inline void Dimension(const coord<int> &in)	{  dimx=in.x(); dimy=in.y(); dimz=in.z(); }
		inline void dim(const coord<int> &in)		{  dimx=in.x(); dimy=in.y(); dimz=in.z(); }

		inline int size()					{ return dimx*dimy*dimz;	}

//conversions for these are left to data containers (i.e. classes that inherit this guy)
		//for cartisian based grids
		inline double Xmax()			{ return max_.x();	}
		inline double Xmin()			{ return min_.x();	}
		inline double Ymax()			{ return max_.y();	}
		inline double Ymin()			{ return min_.y();	}
		inline double Zmin()		 	{ return min_.z();	}
		inline double Zmax()			{ return max_.z();	}

		//for cartisian based grids
		inline void Xmax(double in)			{ return max_.x(in);	}
		inline void Xmin(double in)			{ return min_.x(in);	}
		inline void Ymax(double in)			{ return max_.y(in);	}
		inline void Ymin(double in)			{ return min_.y(in);	}
		inline void Zmin(double in)		 	{ return min_.z(in);	}
		inline void Zmax(double in)			{ return max_.z(in);	}

		//for spherical based grids
		inline void Rmax(double in)			{ return max_.x(in);	}
		inline void Rmin(double in)			{ return min_.x(in);	}
		inline void Thetamax(double in)			{ return max_.y(in);	}
		inline void Thetamin(double in)			{ return min_.y(in);	}
		inline void Phimin(double in)		 	{ return min_.z(in);	}
		inline void Phimax(double in)			{ return max_.z(in);	}

		//for cylindrical based grids
		inline double Rmax()			{ return max_.x();	}
		inline double Rmin()			{ return min_.x();	}
		inline double Phimax()			{ return max_.z();	}
		inline double Phimin()			{ return min_.z();	}

		//for spherical based grids
		inline double Thetamax()		{ return max_.y();	}
		inline double Thetamin()		{ return min_.y();	}

		//dimensions for each side (carteisian)
		inline int Xdim()			 	{ return dimx;	}
		inline int Ydim()				{ return dimy;	}
		inline int Zdim()				{ return dimz;	}

		//dimensions for each side (carteisian)
		inline void Xdim(int in)			{ dimx=in;	}
		inline void Ydim(int in)			{ dimy=in;	}
		inline void Zdim(int in)			{ dimz=in;	}
		//dimensions for each side (spherical )
		inline void Rdim(int in)			{ dimx=in;	}
		inline void Thetadim(int in)		{ dimy=in;	}
		inline void Phidim(int in)			{ dimz=in;	}

		//dimensions for each side (cylindircal)
		inline int Rdim()			 	{ return dimx;	}
		inline int Phidim()			 	{ return dimy;	}

		//dimensions for each side (spherical)
		inline int Thetadim()			{ return dimy;	}

		inline coord<> Center()			{ return center_;	}

		inline void reset()
		{
			notended_=true;
			iterx=0;itery=0;iterz=0;
			xended_=true;
			yended_=true;
			zended_=true;
		}

		inline void Xreset()
		{
			iterx=0;
			xended_=true;
		}
		inline void Yreset()
		{
			itery=0;
			yended_=true;
		}
		inline void Zreset()
		{
			iterz=0;
			zended_=true;
		}

		inline void GridReset()	{	reset();	}


		bool operator==(GeneralGrid &rhs)
		{
			if(this==&rhs) return true;
			if(min_!=rhs.min_) return false;
			if(max_!=rhs.max_) return false;
			if(dimx!=rhs.dimx) return false;
			if(dimy!=rhs.dimy) return false;
			if(dimz!=rhs.dimz) return false;
			if(center_!=rhs.center_) return false;
			return true;
		}

		bool operator!=(GeneralGrid &rhs)
		{
			return (*this)==rhs? false:true;
		}



		/* binary write::
			file out format
			grid<size><xdim><ydim><zdim><minx><miny><minz><maxz><maxy><maxz>
		*/
		bool write(std::fstream &oo)
		{
			if(!oo || oo.fail()){
				Fileerr(); return false;
			}
			BinaryWriteString("grid", oo);
			int siz=size();	oo.write((char *)&siz, sizeof(int));
			oo.write((char *)&dimx, sizeof(int));
			oo.write((char *)&dimy, sizeof(int));
			oo.write((char *)&dimz, sizeof(int));
			double ll=Xmin();	oo.write((char *)&ll, sizeof(double));
			ll=Ymin();	oo.write((char *)&ll, sizeof(double));
			ll=Zmin();	oo.write((char *)&ll, sizeof(double));
			ll=Xmax();	oo.write((char *)&ll, sizeof(double));
			ll=Ymax();	oo.write((char *)&ll, sizeof(double));
			ll=Zmax();	oo.write((char *)&ll, sizeof(double));
			return true;
		}

		bool read(std::fstream &oo)
		{
			if(!oo || oo.fail()){
				Fileerr(); return false;
			}

			std::string gr=BinaryReadString(oo, 4);
			if(gr!="grid") return false;
			int siz=0;
			oo.read((char *)&siz, sizeof(int));
			oo.read((char *)&dimx, sizeof(int));
			oo.read((char *)&dimy, sizeof(int));
			oo.read((char *)&dimz, sizeof(int));
			oo.read((char *)&min_.x(), sizeof(double));
			oo.read((char *)&min_.y(), sizeof(double));
			oo.read((char *)&min_.z(), sizeof(double));
			oo.read((char *)&max_.x(), sizeof(double));
			oo.read((char *)&max_.y(), sizeof(double));
			oo.read((char *)&max_.z(), sizeof(double));
			return true;

		}

};


class GeneralGridIter {
	private:

	public:
		GeneralGrid *myg_;

		int iterx;	//itterator steps
		int itery;
		int iterz;
		int curpos_;	//the count down the list
		bool notended_;			//have we looped through the entire list?
		bool xended_;			//looped through X
		bool yended_;			//looped through Y
		bool zended_;			//looped through Z
		coord<> Pt_;

		GeneralGridIter() :
			myg_(NULL),
			iterx(0), itery(0), iterz(0),curpos_(0),
			notended_(true), xended_(true), yended_(true), zended_(true),
			Pt_(0,0,0)
		{};

		GeneralGridIter(GeneralGrid &in):
			iterx(0), itery(0), iterz(0),curpos_(0),
			notended_(true), xended_(true), yended_(true), zended_(true),
			Pt_(0,0,0)
		{
			myg_=&in;
		}

		GeneralGridIter(GeneralGridIter &cp):
			iterx(cp.iterx), itery(cp.itery), iterz(cp.iterz),curpos_(cp.curpos_),
			notended_(cp.notended_), xended_(cp.xended_), yended_(cp.yended_), zended_(cp.zended_),
			Pt_(cp.Pt_)
		{
			myg_=cp.myg_;
		}

		~GeneralGridIter()
		{
			myg_=NULL;
		}

		inline void reset()
		{
			notended_=true;
			iterx=0;itery=0;iterz=0;
			curpos_=0;
			xended_=true;
			yended_=true;
			zended_=true;
		}

		inline void Xreset()
		{
			iterx=0;
			xended_=true;
		}
		inline void Yreset()
		{
			itery=0;
			yended_=true;
		}
		inline void Zreset()
		{
			iterz=0;
			zended_=true;
		}

		inline void SetGrid(GeneralGrid &in)
		{
			myg_=&in;
		}


		inline operator bool()	{ return notended_;	}
		inline bool GridEnd()		{ return notended_;	}
		inline bool XEnd()			{ return xended_;	}
		inline bool YEnd()			{ return yended_;	}
		inline bool ZEnd()			{ return zended_;	}

		GeneralGridIter &operator=(const GeneralGridIter &rhs)
		{
			if(this==&rhs) return *this;
			myg_=rhs.myg_;
			iterx=rhs.iterx;
			itery=rhs.itery;
			iterz=rhs.iterz;
			notended_=rhs.notended_;
			xended_=rhs.xended_;
			yended_=rhs.yended_;
			zended_=rhs.zended_;
			curpos_=rhs.curpos_;
			return *this;
		}
};


END_BL_NAMESPACE




#endif


