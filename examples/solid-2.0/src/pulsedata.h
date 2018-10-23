
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-25-02
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

#ifndef __Pulse_Data_h__
#define __Pulse_Data_h__ 1

#include "blochlib.h"
using namespace BlochLib;
using namespace std;


class PulseData{
		private:
			void PulseErr(std:: string,const char *, int line);
			void DelayErr(std::string, const char *, int line);
		public:
			double rotorangle;

			Vector<double> angle; //pulses amps for each spin type
			Vector<double> phase;
			Vector<double> offset;

			Vector<double> amp;

			double t1; //start time
			double t2; //end time

			Vector<double> delay; //of >0 will provide no pulse
			Vector<double> dtDelay;
			Vector<bool> incDelay;

			Vector<std::string> on;

			std::string usesys;
			std::string usepow;
			std::string cycler;

		//this is a ptr to a Parser object
		// so that we can do global or local variable translation
		// another class should set this
		// if the expression in the iunput string are going
		// to be expression
			Parser *myParse;

		//teh empty constructor
			PulseData();

		//create a pulse data from a pulse input
			PulseData(double ina, double inp, double ino,double inamp, std::string inon);

		//create a Pulse data from a delay
			PulseData(double del);

		//copy constructor
			PulseData(const PulseData &rhs);

		//destructor
			~PulseData();

		//the assignment operator
			void operator=(const PulseData &rhs);

		//calculates the time for a pulse or a dely input
			void calcTime(double start=0.0);

		//standard pulse on a 'spin'
			hmatrixs Pulse(SolidSys &A);

		//simple printer
			void print(std::ostream &oo,bool printhead=false) const;
			static void printHead(std::ostream &oo) ;

		//reads a pulse line...
			bool read(std::string inSec);
};

//simple Text based ouput of a pulse
std::ostream &operator<<(std::ostream &oo,const PulseData &out);
std::ostream &operator<<(std::ostream &oo, Vector<PulseData> &out);


//Reads the Pulse Data vector from an input Parameter set
// given the pulse spin 'the section name' on
Vector<PulseData> readPulses(Parameters &pset, Parser *myP=NULL);

//Reads the Pulse Data vector from an input Parameter set
// given the pulse spin 'the section name' on
Vector<PulseData> readPulses(const Vector<std::string> &pset, Parser *myP=NULL);

#endif

