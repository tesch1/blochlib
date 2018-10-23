
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
	unitrain.cc--> uniform Time Train Engine
*/

#ifndef _uniTrain_Meth_h_
#define _uniTrain_Meth_h_ 1

#include "timetrain/unitrain.h"
#include "timetrain/gentrain.h"
#include "timetrain/consttrain.h"
#include "container/Vector/Vector.h"
#include "utils/utils.h"


BEGIN_BL_NAMESPACE


/******************NON constant Step Uniform Time Engine Class*****************/

int UniformTimeEngine::TimeFunc(Vector<double> &tO_, Vector<double> &dtO_, Vector<int> &dstepO_)
{
	if(nele_<=0) NstepErr();
	tO_.resize(nele_);
	dtO_.resize(nele_);
	dstepO_.resize(nele_);
	tO_(0)=beginT_;
	dt_=(endT_-beginT_)/(nele_-1);
	dtO_=dt_;
	dstepO_.fill(dstep_);
	for(int i=1;i<nele_;i++){
		tO_(i)=tO_(i-1)+dt_;
	}
	return 1;
}




UniformTimeEngine & UniformTimeEngine::operator =(const UniformTimeEngine &rhs)
{
	if(this==&rhs) return *this;
	GenTimeEngine::operator=(rhs);
	beginT_=rhs.beginT_;
	nele_=rhs.nele_;
	endT_=rhs.endT_;
	dstep_=rhs.dstep_;
	return *this;
}

template<int nste>
UniformTimeEngine & UniformTimeEngine::operator =(const ConstUniformTimeEngine<nste> &rhs)
{
	if(this==&rhs) return *this;
	beginT_=rhs.beginT_;
	nele_=nste;
	endT_=rhs.endT_;
	dstep_=rhs.dstep_;
	return *this;
}


void UniformTimeEngine::addStep()
{
	endT_+=dt_;
	nele_++;
}

void UniformTimeEngine::dropStep()
{
	endT_-=dt_;
	nele_--;
}

void  UniformTimeEngine::setStep(int st)
{
	if(st==nele_) return;
	nele_=st;
}

void UniformTimeEngine::setBeginTime(double st)
{
	if(st==beginT_) return;
	if(st<endT_) beginT_=st;
}

void UniformTimeEngine::setEndTime(double st)
{
	if(st==endT_) return;
	if(st>beginT_) endT_=st;
}


void UniformTimeEngine::setSubStep(int newSub)
{
	if(newSub!=0 && newSub>0){
		dstep_=newSub;
	}else{
		dstep_=1;
	}
}

//the Binary Writing function
//writes the following format to the file...

//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

int UniformTimeEngine::write(std::fstream &oo)
{
	static std::string sflag="UniformTimeEngine";
	static std::string esflag="endUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::write(std::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	oo<<sflag<<std::endl;
	oo.write((char *)&beginT_, sizeof(double));
	oo.write((char *)&endT_, sizeof(double));
	oo.write((char *)&nele_, sizeof(int));
	oo.write((char *)&dstep_, sizeof(int));
	oo<<std::endl<<esflag<<std::endl;
	return 1;
}

//this is the ASCII writer function
//write the following format
//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine
//but in ASCII not binary

int UniformTimeEngine::writeASCII(std::ofstream &oo)
{
	static std::string sflag="UniformTimeEngine";
	static std::string esflag="endUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::print(ostd::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	oo<<std::endl<<sflag<<std::endl;
	oo<<beginT_<<" "<<endT_<<" "<<nele_<<" "<<dstep_<<std::endl;
	oo<<esflag<<std::endl;
	return 1;
}



void UniformTimeEngine::print(std::ostream &oo)
{
	static std::string sflag="UniformTimeEngine";
	static std::string esflag="endUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::print(ostd::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return;
	}
	oo<<std::endl<<sflag<<std::endl;
	oo<<beginT_<<" "<<endT_<<" "<<nele_<<" "<<dstep_<<std::endl;
	oo<<esflag<<std::endl;
}

//reading function reads a ASCII FILE in the following format

//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

int UniformTimeEngine::readASCII(std::ifstream &in)
{
	static std::string sflag="UniformTimeEngine";
	static std::string esflag="endUniformTimeEngine";
	if(!in){
		std::cerr<<std::endl<<"Error: UniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	char tm[1000];
	std::string temp="";
	while(temp!=sflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: UniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
		std::cerr<<" No 'UniformTimeTrain' found...."<<std::endl;
		return 0;
	}
	Vector<std::string> tmm;
	in.getline(tm, 1000, '\n');
	tmm=parse_param(std::string(tm));
	if(tmm.size()!=4)
	{
		std::cerr<<std::endl<<"Error: UniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
		std::cerr<<" invalid data line in the ASCII file..."<<std::endl;
		std::cerr<<" nothing can be read..."<<std::endl;
		return 0;
	}
	beginT_=atof(tmm[0].c_str());
	endT_=atof(tmm[1].c_str());
	nele_=atoi(tmm[2].c_str());
	dstep_=atoi(tmm[3].c_str());
	while(temp!=esflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: UniformTimeEngine::readASCII(std::fstream)"<<std::endl;
		std::cerr<<" No _End_ for the 'UniformTimeTrain' found...."<<std::endl;
		return 0;
	}
	return 1;
}


//reading function reads a BINARY FILE in the following format
//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine

int UniformTimeEngine::read(std::fstream &in)
{
	static std::string sflag="UniformTimeEngine";
	static std::string esflag="endUniformTimeEngine";
	if(!in){
		std::cerr<<std::endl<<"Error: UniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	char tm[1000];
	std::string temp="";
	while(temp!=sflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: UniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" No 'UniformTimeEngine' found...."<<std::endl;
		return 0;
	}
	in.read((char *)&beginT_, sizeof(double));
	in.read((char *)&endT_, sizeof(double));
	in.read((char *)&nele_, sizeof(int));
	in.read((char *)&dstep_, sizeof(int));
	while(temp!=esflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: UniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" No _End_ for the 'UniformTimeEngine' found...."<<std::endl;
		return 0;
	}
	return 1;
}



int UniformTimeEngine::binarySize()
{
	return (17+20)*sizeof(char)+2*sizeof(double)+2*sizeof(int);
}

END_BL_NAMESPACE


#endif

