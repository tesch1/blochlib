

#ifndef _grids_h_
#define _grids_h_


#include "container/grids/gengrid.h"
#include "container/grids/unigrid.h"
#include "container/grids/fullspheregrid.h"
#include "container/grids/surfspheregrid.h"
#include "container/grids/surfcubegrid.h"
#include "container/grids/fullcubegrid.h"
#include <iostream>
#include <fstream>

BEGIN_BL_NAMESPACE


//A Null grid (empty class container to act as some default widget leter)
class NullGrid {};

template<class Engine_t>
class Grid : public Engine_t {

	public:

		typedef typename Engine_t::iterator iterator;

		Grid(): Engine_t() {}

		Grid(const Grid &cp) :
			Engine_t(cp)
		{}

		Grid(double min, double max, int dims):
			Engine_t(min, max, dims)
		{}

		Grid(const coord<> &min,const coord<> &max,const coord<int> &dims):
			Engine_t(min, max, dims)
		{}

		Grid(const Engine_t &dup):
			Engine_t(dup)
		{}

		Grid(double minr, double maxr, double minth, double maxth, int dims):
			Engine_t(minr, maxr, minth, maxth, dims)
		{}

		Grid(double minr, double maxr, double minth, double maxth,double minph, double maxph, int dims):
			Engine_t(minr, maxr, minth, maxth, minph, maxph, dims)
		{}

		Grid(double maxr, double minth, double maxth, int dims):
			Engine_t(maxr, minth, maxth, dims)
		{}

		Grid(double maxr, double minth, double maxth,double minph, double maxph, int dims):
			Engine_t(maxr, minth, maxth, minph, maxph, dims)
		{}

		Grid(double maxr, double maxth, double maxph, int dimx, int dimy, int dimz):
			Engine_t(maxr, maxth, maxph, dimx, dimy, dimz)
		{}

		Grid(double maxr, double minth, double maxth,double minph, double maxph, int dimx, int dimy, int dimz):
			Engine_t(maxr, minth, maxth, minph, maxph, dimx, dimy, dimz)
		{}

		Grid(double minr, double maxr, double minth, double maxth,double minph, double maxph, int dimx, int dimy, int dimz):
			Engine_t(minr, maxr, minth, maxth, minph, maxph, dimx, dimy, dimz)
		{}

		~Grid(){}
};


template<class Engine_t>
std::ostream &operator<<( std::ostream &oo, Grid<Engine_t> &out){
	typename Grid<Engine_t>::iterator myIt(out);
	while(myIt){
		//coord<> moo=out.Point();
		//moo.ToSphere();
		oo<<myIt.Point()<<endl;
		++myIt;
	}
	return oo;
}


END_BL_NAMESPACE


#endif



