


/* blochParams.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-30-01
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
 	blochParams.h-->the basic storage block for the SINGLE spin..

 	There are 2 TYPES the 'Particle' BlochParam
 	and the 'Density' BlochParam..

 	the Density BlochParam behaves as if the spin is a density inside
 	a small volume...i.e. the Equilbrium magnitization is a function of
 	temperature, Magnetic field strength, etc

 	the Particle behave like a single spin...the equilibrium magnitization
 	is +1...in order for there to be a 'bulk' type property there must be a
 	list manager (the ListBlochParams)

 	this flag is templated...the class 'BlochParamType' holds the
 	decision enum...the default is a Density
 */



#ifndef _bloch_params_h_
#define _bloch_params_h_ 1



#include "bloch/Isotope.h"
#include "container/grids/coords.h"
#include "container/matrix/matrix.h"
#include "bloch/blochParamsBasic.h"


BEGIN_BL_NAMESPACE


template<int Type=BPoptions::Density | BPoptions::HighField, class Offset_T=double>
class BlochParams;


//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation
// DENSITY!!
template<>
class BlochParams<BPoptions::Density | BPoptions::HighField, double>:
	public BasicBlochParams<double>
{
	friend class BlochParams<BPoptions::Particle |  BPoptions::HighField, double>;
	private:
		double offset_;		// chemical shift 'offset' the offset along x,y,z

	public:

		BlochParams();
		BlochParams(std::string spin);
		BlochParams(const BlochParams &copy);

		void calcMo();
		void calcMoNorm();

		void Mo(double inMo);
		void Bo(double inBo);
		void temperature(double Tempin);
		void moles(double inmoles);

		inline double Mo(){	return	BasicBlochParams<double>::Mo();	}
		inline double Bo(){	return	BasicBlochParams<double>::Bo();	}
		inline double temperature(){	return	BasicBlochParams<double>::temperature();	}
		inline double moles(){	return	BasicBlochParams<double>::moles();	}

		void setInitialCondition(int IC, int HalfFlag=0);

		void operator=(const BlochParams &rhs);
		BlochParams &operator=(const std::string &rhs);

		bool operator==(const BlochParams &rhs);
		bool operator!=(const BlochParams &rhs);

		void print(std::ostream &oo) const;
		inline void print() const { print(std::cout);	}	//text print
		bool write(std::fstream &oo)const;	//binary write
		int read(std::fstream &in);	//binart read
		int binarySize() const;				//number of 'chars' the data fills
};

/************** PARTICEL >>> SINGLE OFFSET *****/

//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation
// PARTICLE!!
template<>
class BlochParams<BPoptions::Particle |  BPoptions::HighField, double>:
	public BlochParams<BPoptions::Density |  BPoptions::HighField, double>
{
	//friend class ListBlochParams<NullGrid>;
	public:

		BlochParams();
		BlochParams(std::string spin);
		BlochParams(const BlochParams &copy);

		void calcMo();
		void calcMoNorm();

		void Mo(double inMo);
		void Bo(double inBo);
		void temperature(double Tempin);
		void moles(double inmoles);

		void operator=(const BlochParams &rhs);
		BlochParams &operator=(const std::string &rhs);

		void setInitialCondition(int IC, int HalfFlag=0);

		bool operator==(const BlochParams &rhs);
		bool operator!=(const BlochParams &rhs);

		void print(std::ostream &oo) const;
		inline void print() const { print(std::cout);	}	//text print

};

/****DENSITY  COORD<> OFFFSET TYPES ******/

//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation
// DENSITY AND COORD<> OFFSET TYPE!!
template<>
class BlochParams<BPoptions::Density |  BPoptions::HighField, coord<> >:
	public BasicBlochParams<coord<> >
{
	friend class BlochParams<BPoptions::Particle |  BPoptions::HighField, coord<> >;
	private:
		//coord<> offset_;		// chemical shift 'offset' the offset along x,y,z

	public:

		BlochParams();
		BlochParams(std::string spin);
		BlochParams(const BlochParams &copy);

		void calcMo();
		void calcMoNorm();

		void setInitialCondition(int IC, int HalfFlag=0);

		void Mo(const coord<> &inMo);
		void Bo(const coord<> &inBo);
		void temperature(double Tempin);

		void moles(double inmoles);

		inline coord<> Mo(){	return	BasicBlochParams<coord<> >::Mo();	}
		inline coord<> Bo(){	return	BasicBlochParams<coord<> >::Bo();	}
		inline double temperature(){	return	BasicBlochParams<coord<> >::temperature();	}
		inline double moles(){	return	BasicBlochParams<coord<> >::moles();	}


		void operator=(const BlochParams &rhs);

		template<int Options_I>
		void operator=(const BlochParams<Options_I, double> &rhs)
		{
			if(this==&rhs) return;
			BasicBlochParams<coord <> >::operator=(rhs);
		}

		BlochParams &operator=(const std::string &rhs);

		bool operator==(const BlochParams &rhs);
		bool operator!=(const BlochParams &rhs);

		void print(std::ostream &oo) const;
		inline void print() const { print(std::cout);	}	//text print
		bool write(std::fstream &oo)const;	//binary write
		int read(std::fstream &in);	//binart read
		int binarySize() const;				//number of 'chars' the data fills
};


//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation
// PARTICLE AND A COORD<> OFFSET TYPE!!
template<>
class BlochParams<BPoptions::Particle | BPoptions::HighField, coord<> >:
	public BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
{
	//friend class ListBlochParams<NullGrid>;
	public:

		BlochParams();
		BlochParams(std::string spin);
		BlochParams(const BlochParams &copy);

		void calcMo();
		void calcMoNorm();

		void Mo(const coord<> &inMo);
		void Bo(const coord<> &inBo);
		void Mo(double inMo);
		void Bo(double inBo);

		void temperature(double Tempin);
		void moles(double inmoles);

		void operator=(const BlochParams &rhs);
		BlochParams &operator=(const std::string &rhs);

		bool operator==(const BlochParams &rhs);
		bool operator!=(const BlochParams &rhs);

		void print(std::ostream &oo) const;
		inline void print() const { print(std::cout);	}	//text print

};


template<int Options_I, class Offset_T>
std::ostream &operator<<(std::ostream &oo,const BlochParams<Options_I, Offset_T> &out)
{
	out.print(oo);
	return oo;
}

template<int Options_I, class Offset_T>
std::fstream &operator<<(std::fstream &oo, const BlochParams<Options_I, Offset_T> &out)
{
	out.write(oo);
	return oo;
}

template<int Options_I, class Offset_T>
std::fstream &operator>>(std::fstream &oo, BlochParams<Options_I, Offset_T> &out)
{
	out.read(oo);
	return oo;
}


typedef BlochParams<BPoptions::Density | BPoptions::HighField, double> BlochParam;


END_BL_NAMESPACE

#endif

