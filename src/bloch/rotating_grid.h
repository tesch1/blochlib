/* rotating_grid.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 04.28.02
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
 	rotating_grid.h-->uses another gird class as it base,
 	specifiy and axis and a rotation rate, and it will 'rotate'
 	the grid for a given time...

 */


#ifndef _Rotating_grid_shape_h_
#define _Rotating_grid_shape_h_ 1


BEGIN_BL_NAMESPACE


//pre dec for iterator
template<class ShapeEng_t>
class RotatingGridIter;



template<class ShapeEng_t>
class RotatingGrid :
	public ShapeEng_t
{

	private:
		coord<> axis_; //the rotation axis
		double rate_;	//the rotation rate

	//these are to speed up the rotation calculation
	//the entire grid is calculated once for a given time
		Vector<coord<> > rotdata_;	//the rotated grid data
		double curt_;	//the current time

	//this tells us to ignore the 'optimiztion' step
	// and simply calculate the rotation for a given
	//time..the defualt is false
		bool continous_;

	public:

		friend class RotatingGridIter<ShapeEng_t > ;

		typedef RotatingGridIter<ShapeEng_t >  iterator;
		typedef double numtype;

		const void EmptyWarn();

		RotatingGrid():
			ShapeEng_t(),
			axis_(ZeroType<coord<> >::zero()), rate_(ZeroType<coord<> >::zero()),
			rotdata_(ShapeEng_t::data()),
			curt_(-1e30),
			continous_(false)

		{}

		RotatingGrid(const RotatingGrid &cp):
			ShapeEng_t(cp),
			axis_(cp.axis_),
			rate_(cp.rate_),
			rotdata_(cp.rotdata_),
			curt_(-1e30),
			continous_(false)
		{}

		RotatingGrid(const ShapeEng_t &mygrid,
						const coord<> &ina=ZeroType<coord<> >::zero(),
						const double &rte=0.0,
						bool contin=false ):
			ShapeEng_t(mygrid),
			axis_(ina),
			rate_(rte),
			rotdata_(size()),
			curt_(-1e30),
			continous_(contin)
		{
			axis_/=norm(axis_);
		}

		RotatingGrid(const coord<> &mins, const coord<> &maxs, const coord<int> &dims,
						const coord<> &ina=ZeroType<coord<> >::zero(),
						const double &rte=0.0,
						bool contin=false) :
			ShapeEng_t(mins, maxs, dims),
			axis_(ina), rate_(rte),
			rotdata_(size()),
			curt_(-1e30),
			continous_(contin)
		{
			axis_/=norm(axis_);
		}


		void operator=(const RotatingGrid &rhs)
		{
			if(this==&rhs) return;
			ShapeEng_t::operator=(rhs);
			axis_=rhs.axis_;
			rate_=rhs.rate_;
			rotdata_=rhs.rotdata_;
			continous_=rhs.continous_;
			curt_=rhs.curt_;
		}

		void setAxis(const coord<> &ax)
		{
			axis_=ax;
			if(norm(axis_)!=1.0) axis_/=norm(axis_);
		}

		void setRate(const double &ax)
		{		rate_=ax;	}

		void axis(const coord<> &ax)
		{
			axis_=ax;
			if(norm(axis_)!=1.0) axis_/=norm(axis_);
		}

		void rate(const double &ax)
		{		rate_=ax;	}

		coord<> axis() const
		{
			return axis_;
			if(norm(axis_)!=1.0) axis_/=norm(axis_);
		}

		inline double rate() const
		{		return rate_;	}

		inline bool continous() const
		{		return continous_;	}

		inline void continous(bool in)
		{		 continous_=in;	}

		inline double currentTime() const
		{	return curt_;	}

	//this function rotates all the grid points
		void rotateAll(double angle){
			if(curt_!=angle){
				for(int i=0;i<size();++i){
					rotdata_[i]=rotate(angle, ShapeEng_t::Point(i), axis_);
				}
				curt_=angle;
			}
		}

		void rotate(double angle){	rotateAll(angle);	}

	//the rotation function...a 'vector' way of rotating
		coord<> rotate(double angle,const coord<> &in)
		{	return rotate(angle, in, axis_	);		}

		coord<> rotate(double angle,const coord<> &in, const coord<> &axi)
		{
			double c=cos(angle), s=sin(angle), t=(1.0-c);
			return coord<>(
					in.x()*(c+axi.x()*axi.x()*t)+
					in.y()*(axi.x()*axi.y()*t+axi.z()*s)+
					in.z()*(axi.x()*axi.z()*t-axi.y()*s),

					in.x()*(axi.x()*axi.y()*t-axi.z()*s)+
					in.y()*(c+axi.y()*axi.y()*t)+
					in.z()*(axi.y()*axi.z()*t+axi.x()*s),

					in.x()*(axi.x()*axi.z()*t+axi.y()*s)+
					in.y()*(axi.y()*axi.z()*t-axi.x()*s)+
					in.z()*(c+axi.z()*axi.z()*t)
				   );
		}


	//for time dependant shapes
		inline coord<> Point(int i, int j, int k, double t)
		{
			//if(!continous_){
			//	rotateAll(t*rate_);
			//	return coord<>(rotdata_[i].x(),rotdata_[j].y(),rotdata_[k].z());
			//}else{
				return rotate(t*rate_, ShapeEng_t::Point(i,j,k));
			//}
		}

		inline coord<> operator()(int i, double t)
		{	return Point(i,t);	}

		inline coord<> Point(int i, double t)
		{
			//if(!continous_){
			//	rotateAll(t*rate_);
			//	return rotdata_[i];
			//}else{
				return rotate(t*rate_, ShapeEng_t::Point(i));
			//}
		}

		inline coord<> operator()(int i, int j, int k, double t)
		{	return Point(i,j,k,t);		}
};


template<class ShapeEng_t>
class RotatingGridIter :
	public  ShapeEng_t::iterator
{
	private:
		RotatingGrid<ShapeEng_t> *Rot;

	public:

		typedef typename ShapeEng_t::iterator SIter;
		RotatingGridIter():
			SIter(),Rot(NULL)
		{}

		RotatingGridIter(RotatingGrid<ShapeEng_t> &in):
			SIter(in),Rot(&in)
		{}

		RotatingGridIter(RotatingGrid<ShapeEng_t> *in):
			SIter(in),Rot(in)
		{}

		RotatingGridIter(RotatingGrid<ShapeEng_t> &in, Range R):
			SIter(in,R),Rot(&in)
		{}

		RotatingGridIter(RotatingGrid<ShapeEng_t> *in, Range R):
			SIter(in,R),Rot(in)
		{}

		~RotatingGridIter(){	Rot=NULL;	}

		inline coord<> Point(double t)
		{
			//Rot->rotateAll(t*Rot->rate_);
			//return Rot->rotdata_[curpos()];
			return Rot->Point(curpos(), t);
		}

};

END_BL_NAMESPACE


#endif


