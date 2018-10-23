/* mag_fields.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-15-01
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
 	mag_fields.h-->This acts as the base class and sets of 'normal'
 	Magnetic Field containers that Are meant to be used as the
 	"BCalculator" for the Offset<BCalculator, Offset_T> class

 	Nessesary Functions AND constants for this thing to interface nicely

 	1)
 	//returns the vector of the Bfield at time t and grid index 'Index'
 	// It should in in GUASS
 	coord<> Bfield(double t, int Index);

 	2) IF you are using a STATIC field function base the class off  of 'StaticField'
 	   IF you are usning a DYNAMIC (changes in time) base the class off of 'DynamicField'


 */


#ifndef _Mag_Fields_h_
#define _Mag_Fields_h_ 1

#include <string.h>
#include <iostream>
#include "utils/constants.h"
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"


BEGIN_BL_NAMESPACE



class StaticField{
	private:
		Vector<coord<> > Bfield_;

	public:

		static const int Dynamic=0; //not dynamic
		typedef coord<> Offset_T;
//constructors
		StaticField(){}

		StaticField(int size, double Bz):
			Bfield_(size, coord<>(0.0,0.0,Bz))
		{}

		StaticField(int size, double Bx, double By, double Bz):
			Bfield_(size, coord<>(Bx, By, Bz))
		{}

		StaticField(int size, coord<> B):
			Bfield_(size, B)
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		StaticField(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, double Bz):
			Bfield_(bc.size(), coord<>(0.0,0.0,Bz))
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		StaticField(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, double Bx, double By, double Bz):
			Bfield_(bc.size(), coord<>(Bx, By, Bz))
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		StaticField(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, coord<> B):
			Bfield_(bc.size(), B)
		{}

		StaticField(const Vector<coord<> > &Bs):
			Bfield_(Bs)
		{}

		template<class GridEngine_t, int BPops, class Offset_T>
		StaticField(ListBlochParams<GridEngine_t, BPops, Offset_T> &bc, const Vector<coord<> > &Bs):
			Bfield_(Bs)
		{}


//Element acsess
		inline coord<> &operator()(double t, int indx)
		{	return Bfield_[indx];	}

		inline coord<> operator()(double t, int indx) const
		{	return Bfield_[indx];	}

		inline coord<> &operator()(int indx)
		{	return Bfield_[indx];	}

		inline coord<> operator()(int indx) const
		{	return Bfield_[indx];	}

		inline coord<> &Bfield(int indx)
		{	return Bfield_[indx];	}

		inline coord<> Bfield(int indx) const
		{	return Bfield_[indx];	}

		inline coord<> &Bfield(double t, int indx)
		{	return Bfield_[indx];	}

		inline coord<> Bfield(double t, int indx) const
		{	return Bfield_[indx];	}

//Vector Grabber
		inline Vector<coord<> > &data(){	return Bfield_;	}
		inline Vector<coord<> > data() const {	return Bfield_;	}

//aux functions
		inline int size() const {	return Bfield_.size();	}

};


//because it is all but impossible to tell 'what'
//the dynamic class will contain, there is very little
// i can set here...only the nessesary DUMMY function
// 'Bfield' ans set the correct 'Dynamic' Flag

class DynamicField{

	public:

		static const int Dynamic=1; //not dynamic
		typedef coord<> Offset_T;

		DynamicField(){}


		inline coord<> &operator()(double t, int indx)
		{	return ZeroType<coord<> >::zero();	}

		inline coord<> operator()(double t, int indx) const
		{	return ZeroType<coord<> >::zero();	}


		inline coord<> Bfield(double t, int indx) const
		{	return ZeroType<coord<> >::zero();	}

		inline coord<> &Bfield(double t, int indx)
		{	return ZeroType<coord<> >::zero();	}

};


END_BL_NAMESPACE


#endif



