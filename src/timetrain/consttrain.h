

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-11-01
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
	consttrain.h--> header for the constant uniform time train Engine
*/


#ifndef _consttrain_h_
#define _consttrain_h_ 1

#include "container/Vector/Vector.h"
#include "timetrain/gentrain.h"


BEGIN_BL_NAMESPACE



/* CONSTANT uniform dt array defs
		beginT--> the start time
		endT-->the end time
		nsteps-->number of divisions between the beginT and endT
		dstep-->number of divisions WITHIN a (endT-beingT)/nsteps

*/
template<int nsteps>
class ConstUniformTimeEngine: public GenTimeEngine{
	private:
		double dt_;
		double beginT_;
		double endT_;
		int dstep_;

	public:
		ConstUniformTimeEngine():
			GenTimeEngine(nsteps), beginT_(0), endT_(0), dstep_(1)
		{}

		ConstUniformTimeEngine(double endT):
			GenTimeEngine(nsteps), beginT_(0), endT_(endT), dstep_(1)
		{}

		ConstUniformTimeEngine(double beginT,double endT):
			GenTimeEngine(nsteps), beginT_(beginT), endT_(endT), dstep_(1)
		{}

		ConstUniformTimeEngine(double beginT,double endT, int dstep):
			GenTimeEngine(nsteps), beginT_(beginT), endT_(endT), dstep_(dstep)
		{}

		~ConstUniformTimeEngine(){};

/************************************
* The Master Function...the rest of everything here
* is 'fluff' To make the TimeTrain work, all you really need
* is this function....takes 3 vectors
* tO_--> The times starting at beginT---endT
* dtO_--> the time steps between each step in tO
* dstepO_--> the 'sub steps' for each dt step
*****************************/
		int TimeFunc(Vector<double> &tO_, Vector<double> &dtO_, Vector<int> &dstepO_);

/***********Functions for 'after' initialization modification
* theses are called from the Master Class "TimeTrain"
* and are thus 'overridden' by those function...typically
* after the modification the time lists must be recalulated
*************************************/

		ConstUniformTimeEngine &operator =(const ConstUniformTimeEngine &rhs);


//sets the steps WITHIN a dt CANNOT be set to '0'  will be set to '1' instead
// AVOIDS the overflow errors in 'stepsize'
		void setSubStep(int newdt);

//set the begin Time after initialization
		void setBeginTime(double st);

//set the end Time after initialization
		void setEndTime(double st);


//reading function reads a BINARY FILE in the following format
//ConstUniformTimeEngine
//<beingT><endT><dstep>
//endConstUniformTimeEngine

		int read(std::fstream &in);

//reading function reads a ASCII FILE in the following format

//ConstUniformTimeEngine
//<beingT><endT><dstep>
//endConstUniformTimeEngine

		int readASCII(std::ifstream &in);
		int read(std::ifstream &in){	return readASCII(in);	}

//the Binary Writing function
//writes the following format to the file...

//ConstUniformTimeEngine
//<beingT><endT><dstep>
//endConstUniformTimeEngine

		int write(std::fstream &in);

//this is the ASCII writer function
//write the following format

//ConstUniformTimeEngine
//<beingT><endT><nstep><dstep>
//endConstUniformTimeEngine

//but in ASCII not binary

		int writeASCII(std::ofstream &out);
		int write(std::ofstream &out){	return writeASCII(out);	}
		void print(std::ostream &out);

//binary size=(22+25)*sizeof(char)+2*sizeof(double)+2*sizeof(int)
int binarySize();


};

END_BL_NAMESPACE

#include "timetrain/consttrain_meth.h"

#endif



