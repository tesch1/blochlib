
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
	timetrain.h--> The Master Time Train class
	also contains the Master Time Train iterator
*/

/*a class that holds a simple range list
 t starts at 0, and any 'added' time acts as break in time
 and can signify 'anything'  it allows for easy looping over
 certain time steps (espeically if the time steps are not uniform)

 it also holds a "dt" step, where you can specify the number of steps
 to go from "t1-->t2"
*/


#ifndef _timetrain_h_
#define _timetrain_h_ 1

#include <iostream>
#include "timetrain/consttrain.h"
#include "timetrain/unitrain.h"
#include "timetrain/gentrain.h"
#include "timetrain/filetrain.h"
#include <string>

BEGIN_BL_NAMESPACE


template<class Engine_t>
class TimeTrainIter;

/* MAIN Time Train Class */
template<class Engine_t>
class TimeTrain : public Engine_t{
	friend class TimeTrainIter<Engine_t>;
	private:
		Vector<double> t_;
		Vector<double> dt_;
		double beginT_;
		double endT_;
		int nele_;
		Vector<int> dstep_;

		const void InitErr();

		void init();


	public:
		typedef TimeTrainIter<Engine_t> iterator;

		TimeTrain():
			Engine_t(),
			t_(0,0), dt_(0,0), beginT_(0), endT_(0), nele_(1), dstep_(0,1)
		{
			init();
		}

		TimeTrain(const TimeTrain &copy):
			Engine_t(copy)
		{
			init();
		}

		TimeTrain(double it1, double it2):
			Engine_t(it1, it2)
		{
			init();
		}

		TimeTrain(double it1, double it2, int it3):
			Engine_t(it1, it2, it3)
		{
			init();
		}

		TimeTrain(double it1, double it2, int it3, int it4 ):
			Engine_t(it1, it2, it3, it4)
		{
			init();
		}

		TimeTrain(Engine_t in):
			Engine_t(in)
		{
			init();
		}

		void calculate(){	init();	}

/*********Set The Engine To A New one of the smae type
* REQUIRES the 'operator=' be defined in the Engine
*
*/
		template<class NewEng>
		void setEngine(const NewEng &in);

/******************************************
* These acses the vector elements on a
* 'on' call basis.....
* These should be ENGINE DEPENDANT
*/

//the delta time at position 'i'
		double dt(int i) const;

//the Very Begining time (typically 0)
		double beginTime() const;

//the Very End Time (i.e. the last element in the time list)
		double endTime() const ;

//time at pos 'i'
		double t(int i);
//current steps WITHIN a dt at pos 'i' (this will always be the same for this engine)
		int substep(int i) const;

//returns dt/dstep (this will always be the same for this engine)
		double stepsize(int i) const;

//the begining time for the step at pos 'i'
		double beginStepTime(int i) const;

//the end time for the at pos 'i+1'
		double endStepTime(int i) const;


/*******************************************
*  The functions MUST be defined in the
* Inherited class "Engine_t"
* Lucky for us...these are no expressed until
* they are used by the user...so if they do not
* exsist in the base class, and are not used
* there will be no consequences
*/

//reads in a file...(typically used for the FileTimeTrain...but can be used
// for any other type you wish to write a read function for
//'rereads' any exssting files already inputted into the Engine
		void read(){	init();	}
		void read(const char *in){	read(std::string(in)); 	}
		void read(std::string in)	{ Engine_t::read(in); init();	}

//Ascii reader...the ABOVE functions do read ASCII, but it is assumed that
//the data in those files is ONLY for that parameter
//this function allows a file with MORE THEN ONE parameter set to be read
//each Engine has a different file format (or can)
		int readASCII(std::ifstream &in)
		{
			int t=Engine_t::readASCII(in);
			if(t){ init(); return 1;	}
			std::cerr<<std::endl<<"Error: TimeTrain.read(std::ifstream)"<<std::endl;
			std::cerr<<" cannot read file...."<<std::endl;
			return 0;
		}
		int read(std::ifstream &in){	return readASCII(in);	}


//this is the BINARY reading function...it is past to the Engine for binary reading
//returns 1 if 'happy' and '0' if fail
		int read(std::fstream &in)
		{
			int t=Engine_t::read(in);
			if(t){ init(); return 1;	}
			std::cerr<<std::endl<<"Error: TimeTrain.read(std::fstream)"<<std::endl;
			std::cerr<<" cannot read file...."<<std::endl;
			return 0;
		}


//this is the BINARY writing function...it passes it to the Engin  whose formats
//can varry...the the various tEngine header files for BINARY each file format
		int write(std::fstream &in)
		{
			int t=Engine_t::write(in);
			if(t){ init(); return 1;	}
			std::cerr<<std::endl<<"Error: TimeTrain.write(std::fstream)"<<std::endl;
			std::cerr<<" cannot write file...."<<std::endl;
			return 0;
		}

//this is the ASCII writing function...it passes it to the Engin  whose formats
//can varry...the the various tEngine header files for ASCII each file format
		int writeASCII(std::ofstream &in)
		{
			int t=Engine_t::writeASCII(in);
			if(t){ init(); return 1;	}
			std::cerr<<std::endl<<"Error: TimeTrain.writeASCII(std::fstream)"<<std::endl;
			std::cerr<<" cannot write file...."<<std::endl;
			return 0;
		}
		int write(std::ofstream &out){	return writeASCII(out);	}

//set the steps size after initialization
		void setStep(int st);

//set the begin Time after initialization
		void setBeginTime(double st);

//set the end Time after initialization
		void setEndTime(double st);


//sets the steps WITHIN a dt CANNOT be set to '0'  will be set to '1' instead
// AVOIDS the overflow errors in 'stepsize'
		void setSubStep(int newSub);

//The size of the train
		inline int size()const{	return nele_;	}

//add a step to the array (becuase it is a unifor array, it only increases the array by 'dt')
		void addStep();
		void push_back(){ addStep();	}

//add a step to the array (becuase it is a unifor array, it only increases the array by 'dt')
		void addStep(double newEnd);
		void push_back(double newEnd){ addStep(newEnd);	}

//drop a last step to the array
		void dropStep();
		void pop(){ dropStep();	}

/***************End Aux Funcs ***************/
		template<class Engine_t2>
		TimeTrain &operator=(const TimeTrain<Engine_t2> &rhs);

		TimeTrain &operator=(const TimeTrain &copy);

		template<class Engine_t2>
		TimeTrain &operator=(const Engine_t2 &copy);

		~TimeTrain(){}

//these are the ASCII write functions (output formats that the
//Engines can LATER READ
		void print(std::ostream &oo){	Engine_t::print(oo);	}

//prints the FIle Formate Valid for the "FileTimeEngine" Engine
// This ASCII AND is assumed to be the ONLY thing in the file....
		void printFileTimeEngine(std::ostream &oo);

//the display function...prints out the list in its entirity in a nice formated
//display...this is usually used for debugging code...i.e. prints to COUT
//This format CANNOT be read be the Engines
		void display(){	display(cout);	}
		void display(std::ostream &oo);

};

template<class Engine_t>
class TimeTrainIter : public GenTimeIter{
	private:
		TimeTrain<Engine_t> *mye_;

	public:

		TimeTrainIter():
			GenTimeIter(),
			mye_(NULL)
		{}

		//copy constructor
		TimeTrainIter(TimeTrainIter &in):
			GenTimeIter(in),
			mye_(in.mye_)
		{}

		TimeTrainIter(TimeTrain<Engine_t> &in):
			GenTimeIter(in),
			mye_(&in)
		{}

		TimeTrainIter &operator=(const TimeTrainIter &rhs);

		~TimeTrainIter()
		{
			mye_=NULL;
		}

//current time
		double t() const ;
//current dt
		double dt() const;

//current steps WITHIN a dt
		int substep() const;

//returns dt/dstep;
		double stepsize() const;

//the begining time for the CURRENT step
		double beginStepTime() const;

//the end time for the CURRENT step
		double endStepTime() const;

	//increment operator
		void operator++();
		void operator++(int) {	operator++();	}
	//decrament operator;
		void operator--();
		void operator--(int){	operator--();	}
};

//Display operator...
template<class Engine_t>
std::ostream &operator << (std::ostream &oo,TimeTrain<Engine_t> &out);

//ASCII WRITE
template<class Engine_t>
std::ofstream &operator << (std::ofstream &OutFile,TimeTrain<Engine_t> &ObjToWrite);


//BINARY WRITE
template<class Engine_t>
std::fstream &operator << (std::fstream &in, TimeTrain<Engine_t> &out);


//BINARY READ
template<class Engine_t>
std::fstream &operator >> (std::fstream &in, TimeTrain<Engine_t> &out);

//ASCII READ
template<class Engine_t>
std::ifstream &operator >> (std::ifstream &in, TimeTrain<Engine_t> &out);


END_BL_NAMESPACE


#include "timetrain/timetrain_meth.h"

#endif

