

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
	timetrain_meth.h--> The Master Time Train class methods
	also contains the Master Time Train iterator methods
*/


/* methods for Time Train */

#ifndef _TimeTrain_Meth_h_
#define _TimeTrain_Meth_h_ 1

#include "timetrain/timetrain.h"
#include "timetrain/gentrain.h"
#include "timetrain/consttrain.h"
#include "timetrain/unitrain.h"
#include "container/Vector/Vector.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


template<class Engine_t>
const void TimeTrain<Engine_t>::InitErr()
{
	std::cerr<<std::endl<<"Error:: TimeTrain"<<std::endl;
	BLEXCEPTION(" Initialization Error...i must exit the program....")
}


template<class Engine_t>
void TimeTrain<Engine_t>::init()
{
	if(Engine_t::TimeFunc(t_, dt_, dstep_)!=1) InitErr();

	nele_=t_.size();
	beginT_=t_(0);
	endT_=t_(nele_-1);
}

template<class Engine_t>
double TimeTrain<Engine_t>::dt(int i) const
{
	return dt_(i);
}

template<class Engine_t>
int TimeTrain<Engine_t>::substep(int i) const
{
	return dstep_(i);
}

template<class Engine_t>
double TimeTrain<Engine_t>::beginTime() const
{
	return beginT_;
}

template<class Engine_t>
double TimeTrain<Engine_t>::endTime() const
{
	return endT_;
}

template<class Engine_t>
double TimeTrain<Engine_t>::beginStepTime(int i) const
{
	return t_(i);
}

template<class Engine_t>
double TimeTrain<Engine_t>::endStepTime(int i) const
{
	return t_(i+1);
}

template<class Engine_t>
void TimeTrain<Engine_t>::setBeginTime(double st)
{
	Engine_t::setBeginTime(st);
	init();
}

template<class Engine_t>
void TimeTrain<Engine_t>::setEndTime(double st)
{
	Engine_t::setEndTime(st);
	init();
}

template<class Engine_t>
void TimeTrain<Engine_t>::setSubStep(int st)
{
	Engine_t::setSubStep(st);
	init();
}

//add a step to the array (becuase it is a unifor array, it only increases the array by 'dt')
template<class Engine_t>
void TimeTrain<Engine_t>::addStep()
{
	Engine_t::addStep();
	init();
}

//add a step to the array (becuase it is a unifor array, it only increases the array by 'dt')
template<class Engine_t>
void TimeTrain<Engine_t>::addStep(double newEnd)
{
	Engine_t::addStep(newEnd);
	init();
}

//drop a last step to the array
template<class Engine_t>
void TimeTrain<Engine_t>::dropStep()
{
	Engine_t::dropStep();
	init();
}

template<class Engine_t>
template<class Engine_t2>
TimeTrain<Engine_t> &TimeTrain<Engine_t>::operator=(const TimeTrain<Engine_t2> &rhs)
{
	if(this==&rhs) return *this;
	beginT_=rhs.beginT_;
	nele_=rhs.nele_;
	endT_=rhs.endT_;
	this->init();
	return *this;
}

template<class Engine_t>
TimeTrain<Engine_t> &TimeTrain<Engine_t>::operator=(const TimeTrain<Engine_t> &rhs)
{
	if(this==&rhs) return *this;
	Engine_t::operator=(rhs);
	beginT_=rhs.beginT_;
	nele_=rhs.nele_;
	endT_=rhs.endT_;
	this->init();
	return *this;
}

template<class Engine_t>
template<class NewEng>
void TimeTrain<Engine_t>::setEngine(const NewEng &in)
{
	Engine_t::operator=(in);
	init();
}

template<class Engine_t>
template<class Engine_t2>
TimeTrain<Engine_t> &TimeTrain<Engine_t>::operator=(const Engine_t2 &copy)
{
	Engine_t::operator=(copy);
	init();
	return *this;
}

template<class Engine_t>
void TimeTrain<Engine_t>::printFileTimeEngine(std::ostream &oo)
{
	oo<<std::endl<<"%Time train File..."<<std::endl;
	TimeTrain<Engine_t>::iterator MyIt(*this);
	while(MyIt){
		oo<<MyIt.t()<<" "<<MyIt.step()<<std::endl;
		++MyIt;
	}
}

template<class Engine_t>
void TimeTrain<Engine_t>::display(std::ostream &oo)
{
	oo<<"%%TimeTrain....."<<std::endl;
	TimeTrain<Engine_t>::iterator myIt(*this);
	int i=0;
	while(myIt)
	{
		oo<<"step: "<<i++<<" start: "<<myIt.beginStepTime()<<" end: "<<myIt.endStepTime()<<" substeps: "<<myIt.substep()<<std::endl;
		++myIt;
	}
}

//--------------ITerator functions-------------

template<class Engine_t>
void TimeTrainIter<Engine_t>::operator++()
{
	if(!notended_)	mye_->Additerr();
	++curpos_;
	++nextpos_;
	if(nextpos_>=mye_->size()){
		if(curloop_>=loops_)	notended_=false;
		else{	curloop_++; reset();	}
	}
}

template<class Engine_t>
void TimeTrainIter<Engine_t>::operator--()
{
	if(!notended_) mye_->Subiterr();
	--curpos_;
	--nextpos_;
	if(curpos_<0) notended_=false;
}

template<class Engine_t>
double TimeTrainIter<Engine_t>::t() const
{
	return mye_->t_(curpos_);
}

template<class Engine_t>
double TimeTrainIter<Engine_t>::dt() const
{
	return mye_->dt_(curpos_);
}

template<class Engine_t>
int TimeTrainIter<Engine_t>::substep() const
{
	return mye_->dstep_(curpos_);
}

template<class Engine_t>
double TimeTrainIter<Engine_t>::stepsize() const
{
	return mye_->dt_(curpos_)/mye_->dstep_(curpos_);
}

template<class Engine_t>
double TimeTrainIter<Engine_t>::beginStepTime() const
{
	return mye_->t_(curpos_);
}

template<class Engine_t>
double TimeTrainIter<Engine_t>::endStepTime() const
{
	return mye_->t_(nextpos_);
}

template<class Engine_t>
TimeTrainIter<Engine_t> &TimeTrainIter<Engine_t>::operator=(const TimeTrainIter<Engine_t> &rhs)
{
	if(&rhs==this) return *this;
	GenTimeIter::operator=(rhs);
	mye_=rhs.mye_;
	return *this;
}


//----------------------ASCII Ouput TO FILE------------------
template<class Engine_t>
std::ofstream &operator << (std::ofstream &OutFile,TimeTrain<Engine_t> &ObjToWrite)
{
	if(!ObjToWrite.writeASCII(OutFile))
	{
		std::cerr<<std::endl<<"Error: operator<<(TimeTrain)"<<std::endl;
		std::cerr<<" Failed to ASCII write TimeTrain"<<std::endl;
	}
	return OutFile;
}


//----------------------Ouput TO 'pretty' display ------------------
template<class Engine_t>
std::ostream &operator << (std::ostream &OutThing,TimeTrain<Engine_t> &ObjToWrite)
{
	ObjToWrite.display(OutThing);
	return OutThing;
}

//----------------------OUT BINARY ------------------

//BINARY WRITE
template<class Engine_t>
std::fstream &operator << (std::fstream &OutFile, TimeTrain<Engine_t> &ObjToWrite)
{
	if(!ObjToWrite.write(OutFile))
	{
		std::cerr<<std::endl<<"Error: operator <<(TimeTrain)"<<std::endl;
		std::cerr<<" Failed to BINARY write TimeTrain"<<std::endl;
	}
	return OutFile;
}


//----------------------INPUT BINARY ------------------
template<class Engine_t>
std::fstream &operator >> (std::fstream &inFile, TimeTrain<Engine_t> &ObjToRead)
{
	if(!ObjToRead.read(inFile))
	{
		std::cerr<<std::endl<<"Error: operator >> (TimeTrain)"<<std::endl;
		std::cerr<<" Failed to BINARY read the TimeTrain"<<std::endl;
	}
	return inFile;
}

//----------------------INPUT ASCII ------------------
template<class Engine_t>
std::ifstream &operator >> (std::ifstream &inFile, TimeTrain<Engine_t> &ObjToRead)
{
	if(!ObjToRead.readASCII(inFile))
	{
		std::cerr<<std::endl<<"Error: operator >> (TimeTrain)"<<std::endl;
		std::cerr<<" Failed to ASCII read TimeTrain"<<std::endl;
	}
	return inFile;
}

END_BL_NAMESPACE


#endif


