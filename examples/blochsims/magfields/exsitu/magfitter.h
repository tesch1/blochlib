
/* magfitter.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 5.16.02
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
	magfitter.h--> a class to aid in setting
	and adjusting the parameters of coils while
	attempting to 'fit' something about the magnetic field

*/

#ifndef _Mag_Fitter_h_
#define _Mag_Fitter_h_ 1

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//MUST include this to use the minimzer
//#include "minuit/minuit.h"

//simple data container
class MagFitterData
{

	friend class MagFitter;

	public:

		double upperBound; //the upper bound for the param
		double lowerBound; //the lower bound for the param
		double currentValue; //the current value

		double error;	//the 'error bound'

		std::string name;	//the parameters name

	//the allowed alterable params in a Biot
		static const int AltLen=10;
		static const std::string Alterable[AltLen];
		enum ValidAlter{
			R=0,
			Z=1,
			turns=2,
			theta1=3,
			theta2=4,
			amps=5,
			loops=6,
			Xcenter=7,
			Ycenter=8,
			Zcenter=9
		};


		std::string pSection; //the section name in a 'parameters' chunk

	private:
		int alterInt;  //the enum value
		std::string alter; //the name of the thing to alter in Biot
		bool use; //use the param or not?

	public:
		MagFitterData():
			upperBound(0), lowerBound(0), currentValue(0),	error(0),
			name("NA"), pSection(""), alterInt(-1),alter(""), use(false)
		{}

		MagFitterData(Parameters &pset);
		void read(Parameters &pset);

		MagFitterData(const Vector<std::string> &pset);
		void read(const Vector<std::string> &pset);


	//need to check to make sure it is usable
		void setAlter(std::string in);

	//this would be a valid MINUIT string
	//to be used with the function "MNPARS(char *, int err)"
	//the MagFitter supplies the <var #>
		std::string minuitString();
		bool setCoilParams(Parameters &pset);

	//print it out
		void print();
		void print(std::ostream &out);
		void simplePrint(std::ostream &oo);


};

std::ostream &operator<<(std::ostream &oo, MagFitterData &out);


class MagFitter
{
	private:
	//The Coil Parameters chunk
		Parameters CoilPars;

	//the Fitting parameter chunk
		Parameters FitPars;
		std::string subbase; //the base for the multiple params entry

		int numFits; //the number of different parameters
		int numUseable; //number of useabl parameters

	//the actual parameters values
		Vector<MagFitterData> fitparams;

	public:


		MagFitter();
		MagFitter(const Vector<std::string> &Fpset,const  Vector<std::string> &Cpset);
		MagFitter(Parameters &Fpset, Parameters &Cpset);

		void read(const Vector<std::string> &Fpset,const  Vector<std::string> &Cpset);
		void read(Parameters &Fpset, Parameters &Cpset);

		inline int size() const {	return fitparams.size();	}

	//this would be a valid MINUIT string
	//to be used with the function "MNPARS(char *, int err)"
		std::string minuitString(int i);

	//sets the Coil Parameters from the MINUIT
	// variable array
		void setCoilParams(double *vars);

	//is the parameters at index i used
		bool used(int i);
		bool used(std::string name);

		inline Parameters &coilParams() {	return CoilPars;	}

		void print(std::ostream &oo, int header=1);
};

std::ostream &operator<<(std::ostream &oo, MagFitter &out);



#endif

