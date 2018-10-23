

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
	unitrain.h--> uniform Time Train Engine
*/

/* NON CONSTANT uniform dt array defs
		beginT--> the start time
		endT-->the end time
		nsteps-->number of divisions between the beginT and endT
		dstep-->number of divisions WITHIN a (endT-beingT)/nsteps

		Here you can ADD extra steps to the array (the above one, you cannot...but it is faster
		then this one becuase the compile can expand loops on the above one
		And you can drop steps (off the end)

*/

#ifndef _unitrain_h_
#define _unitrain_h_ 1

#include "timetrain/gentrain.h"
#include "timetrain/consttrain.h"
#include <string>
#include <iostream>


BEGIN_BL_NAMESPACE


class UniformTimeEngine: public GenTimeEngine{

	private:
		double dt_;
		double beginT_;
		double endT_;
		int dstep_;


	public:


		UniformTimeEngine():
			GenTimeEngine(1), beginT_(0), endT_(0), dstep_(1)
		{}

		UniformTimeEngine(const UniformTimeEngine &cp):
			GenTimeEngine(cp), beginT_(cp.beginT_), endT_(cp.endT_), dstep_(cp.dstep_)
		{}


		UniformTimeEngine(double endT):
			GenTimeEngine(1), beginT_(0), endT_(endT), dstep_(1)
		{}

		UniformTimeEngine(double beginT,double endT):
			GenTimeEngine(1), beginT_(beginT), endT_(endT), dstep_(1)
		{}

		UniformTimeEngine(double beginT,double endT, int nstep):
			GenTimeEngine(nstep), beginT_(beginT), endT_(endT), dstep_(1)
		{}

		UniformTimeEngine(double beginT,double endT,int nstep, int dstep):
			GenTimeEngine(nstep), beginT_(beginT), endT_(endT), dstep_(dstep)
		{}

		~UniformTimeEngine(){};

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
		UniformTimeEngine &operator =(const UniformTimeEngine &rhs);

		template<int nts>
		UniformTimeEngine &operator =(const ConstUniformTimeEngine<nts> &rhs);

//add a step to the array (becuase it is a unifor array, it only increases the array by 'dt')
		void addStep();
		void push_back(){ addStep();	}

//drop a last step to the array
		void dropStep();
		void pop(){ dropStep();	}

//set the steps size after initialization
		void setStep(int st);

//set the begin Time after initialization
		void setBeginTime(double st);

//set the end Time after initialization
		void setEndTime(double st);


//sets the SSubStep size
//cannot be negative nor 0...either 1 or bigger
		void setSubStep(int newSub);



//reading function reads a BINARY FILE in the following format
//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

		int read(std::fstream &in);

//reading function reads a ASCII FILE in the following format

//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

		int readASCII(std::ifstream &in);
		int read(std::ifstream &in){	return readASCII(in);	}

//the Binary Writing function
//writes the following format to the file...

//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

		int write(std::fstream &in);

//this is the ASCII writer function
//write the following format
//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine
//but in ASCII not binary

		int writeASCII(std::ofstream &out);
		int write(std::ofstream &out){	return writeASCII(out);	}
		void print(std::ostream &out);

//binary size=(16+19)*sizeof(char)+2*sizeof(double)+2*sizeof(int)
		int binarySize();

};


END_BL_NAMESPACE




#endif
