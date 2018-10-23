
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08-06-01
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
	lyapunov.h-->a class to calculate Lyapunov exponents from
	 a chunk of evolved Variational Equation data

	this class is meant to interface with the 'blochsolver' class,
	which contains the following input data format

	 Vector<coord<> > data;


	the 'data' here is set up into 2 disinct chunks
	 suppose there are N differential equations, this means
	 there are N*N Variational equations
	 the data is split up according to

	 difeq 1
	 difeq 2
	 difeq 3
	 ...
	 difeq N
	 varEq 1.1
	 varEq 1.2
	 varEq 1.3
	  ...
	 varEq 1.N
	 varEq 2.1
	 varEq 2.2
	 ....
	 varEq N.N

	for our Bloch Equations each 'difeq' is actually a coordinate of '3' differential equations,
	so this class contains 2 'flavors' the first flavor deals directly with the 'coord<>' strucutre
	the second deals with simple lists of numbers (i.e. if the data was a Vector<double> instead of
	the Vector<coord<> >

	***NOTE::The 'Variational Data' will be treated properly meaning that
	after each calculation of the Exponent at that time step, the data will be
	ORTHONORMALIZED!...to contiue further calculations..

*/

/*  Lyapunov calculations....
		 this simply  calculate the lyapunov exponents
		 from the inegrated variational equations
		 in order for this to function correctly, the 'calcVariational()' function
		 must be specified in the BlochEngine_t parameter set...

		 the class 'blochsolver' will check to make sure this is true

		 'calcstep'--> determins 'how often' to calculate
		 the exponent.  during the simulation process there may be
		 lots of small valid time steps, since the process of lyapunov can be
		 quite expensive, this value can be set to calulate the exponenet every '10' steps
		 if it is set to 0 it will calcualte it at every valid time step

		 'curct'--> the current lyp count (i.e. howmany we have calculated so far)
		 'size'--> the number of DIFFERENTIAL equations (NOT Variational Equations)
		    the class will maintain the correct list counts and posistions

		'ts' --> keeps a running record of 'at what time' was the lyp calulated
		'dts' --> keeps a running record of 'the dts' entcountered in the calculation
		'lyps'--> storage for the exponents themselves

		'data' --> a pointer to the data we wish to analyis.  this is a ptr becuase
		  out diffeq solver maintains the current integrated point

	*/

#ifndef _lyapunov_h_
#define _lyapunov_h_ 1

#include "container/matrix/matrix.h"
#include "container/grids/coords.h"
#include "utils/utils.h"
#include "utils/plotter.h"

BEGIN_BL_NAMESPACE


template<class T >
class Lyapunov{};

template<>
class Lyapunov<coord<> >
{
	private:
		int calcstep_;
		int curct_;
		int pos_;
		int size_;
		Vector<double> ts;
		Vector<double> dts;
		Vector<coord<> > lyps;

		Vector<coord<> > *data_;

		bool continprint_;

	public:

//ctors
		Lyapunov():
			calcstep_(1), curct_(0), pos_(0), size_(0),continprint_(true)
		{}

		Lyapunov(int numLyps):
			calcstep_(1), curct_(0),pos_(0), size_(numLyps),continprint_(true)
		{
			lyps.resize(numLyps,0);
		}

		Lyapunov(int numLyps, int everyNsteps):
			calcstep_(everyNsteps), curct_(0),pos_(0), size_(numLyps),continprint_(true)
		{
			lyps.resize(numLyps,0);
		}

		Lyapunov(int numLyps, Vector<coord<> > &indata, int everyNsteps=1):
			calcstep_(everyNsteps), curct_(0),pos_(0), size_(numLyps),continprint_(true)
		{
			data_=&indata;
			lyps.resize(numLyps,0);
		}

		Lyapunov(int numLyps, Vector<coord<> > *indata, int everyNsteps=1):
			calcstep_(everyNsteps), curct_(0),pos_(0), size_(numLyps),continprint_(true)
		{
			data_=indata;
			lyps.resize(numLyps,0);
		}

		Lyapunov<coord<> > &operator=(Lyapunov<coord<> > &in);

//dtor (null out data_ so the cpu does not try to delete the data contained in it)
		~Lyapunov()
		{
			data_=NULL;
		}

//maintanace functions
		void setData(Vector<coord<> > &indat){		data_=&indat;	}
		void setData(Vector<coord<> > *indat){		data_=indat;	}
		Vector<coord<> > *data(){		return data_;	}
		void setCalcStep(int step) {	calcstep_=(step<=0)?1:step;	}
		void setSize(int insize)
		{
			size_=insize;
			lyps.resize(size_,0);
		}

		//sets the ostream << operator to print the 'current' posistion only
		//(i.e. the 'just' calculated lyapunov)
		void setPrintContinous(){	continprint_=(true);	}


		//sets the ostream << operator to print the entire lyp list posistion
		void setPrintAll(){	continprint_=(false);	}

		void reset();

//zeros out the lyps list
		void Zero(){			lyps.fill(coord<>(0) );		}

		void calcLyapunov(double ont, double dt); //see 'lyapunov_meth.h'

//the exponenets for this particular time step
		Vector<coord<> > Exponents(){	return lyps;	}

//binary writing of the data
// Format::--> <time> <lyp1> <lyp2>...<lypN>
		void write(std::string fname){	write(fname.c_str());	}
		void write(std::fstream &oo);
		void write(const char *fname);

//ASCII write of the data
// Format::--> <time> <lyp1> <lyp2>...<lypN>
		void print(std::string fname){	print(fname.c_str());	}
		void print(std::ostream &oo);
		void print(const char *fname);

};


std::ostream &operator<<(std::ostream &oo, Lyapunov<coord<> > &rhs);

std::fstream &operator<<(std::fstream &oo, Lyapunov<coord<> > &rhs);


END_BL_NAMESPACE


#endif





