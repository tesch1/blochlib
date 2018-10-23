

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
	consttrain_meth.h--> methods for the constant uniform time train
*/

#ifndef _consttrain_meth_h_
#define _consttrain_meth_h_ 1

#ifndef _consttrain_h_
	#error "ConstTrain_meth.h Must be included from 'ConstTrain.h"
#endif

#include "timetrain/gentrain.h"
#include "container/Vector/Vector.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE


/******************Const Step Uniform Time Engine Class*****************/

template<int nsteps>
int  ConstUniformTimeEngine< nsteps>::TimeFunc(Vector<double> &tO_, Vector<double> &dtO_, Vector<int> &dstepO_)
{
	if(nsteps<=0) NstepErr();
	tO_.resize(nsteps);
	dtO_.resize(nsteps);
	dstepO_.resize(nsteps);

	tO_(0)=0;
	dt_=(endT_-beginT_)/(nsteps-1);
	dtO_(0)=dt_;
	dstepO_=dstep_;
	for(int i=1;i<nsteps;i++){
		tO_(i)=tO_(i-1)+dt_;
	}
	return 1;
}


template< int nsteps>
void ConstUniformTimeEngine<nsteps>::setSubStep(int newSub)
{
	if(newSub!=0 && newSub>0)
	{
		dstep_=newSub;
	}else{
		dstep_=1;
	}
}


template< int nsteps>
void ConstUniformTimeEngine<nsteps>::setBeginTime(double st)
{
	if(st==beginT_) return;
	if(st<endT_) beginT_=st;
}

template< int nsteps>
void ConstUniformTimeEngine<nsteps>::setEndTime(double st)
{
	if(st==endT_) return;
	if(st>beginT_) endT_=st;
}


template<int nsteps>
ConstUniformTimeEngine<nsteps> & ConstUniformTimeEngine<nsteps>::operator =(const ConstUniformTimeEngine< nsteps> &rhs)
{
	if(this==&rhs) return *this;
	beginT_=rhs.beginT_;
	endT_=rhs.endT_;
	dstep_=rhs.dstep_;
	return *this;
}

//the Binary Writing function
//writes the following format to the file...

//UniformTimeEngine
//<beingT><endT><nstep><dstep>
//endUniformTimeEngine
template<int nsteps>
int ConstUniformTimeEngine<nsteps>::write(std::fstream &oo)
{
	static std::string sflag="ConstUniformTimeEngine";
	static std::string esflag="endConstUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::write(std::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	oo<<sflag<<std::endl;
	oo.write(&beginT_, sizeof(double));
	oo.write(&endT_, sizeof(double));
	oo.write(&dstep_, sizeof(int));
	oo<<std::endl<<esflag<<std::endl;
	return 1;
}

//this is the ASCII writer function
//write the following format
//UniformTimeEngine
//<beingT><endT><dstep>
//endUniformTimeEngine
//but in ASCII not binary
template<int nsteps>
int ConstUniformTimeEngine<nsteps>::writeASCII(std::ofstream &oo)
{
	static std::string sflag="ConstUniformTimeEngine";
	static std::string esflag="endConstUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::print(ostd::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	oo<<std::endl<<sflag<<std::endl;
	oo<<beginT_<<" "<<endT_<<" "<<dstep_<<std::endl;
	oo<<esflag<<std::endl;
	return 1;
}

template<int nsteps>
void ConstUniformTimeEngine<nsteps>::print(std::ostream &out)
{
	static std::string sflag="ConstUniformTimeEngine";
		static std::string esflag="endConstUniformTimeEngine";
	if(!oo){
		std::cerr<<std::endl<<"Error: UniformTimeTrain::print(ostd::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return;
	}
	oo<<std::endl<<sflag<<std::endl;
	oo<<beginT_<<" "<<endT_<<" "<<dstep_<<std::endl;
	oo<<esflag<<std::endl;
}

//reading function reads a ASCII FILE in the following format

//ConstUniformTimeEngine
//<beingT><endT><dstep>
//endConstUniformTimeEngine
template<int nsteps>
int ConstUniformTimeEngine<nsteps>::readASCII(std::ifstream &in)
{
	static std::string sflag="ConstUniformTimeEngine";
	static std::string esflag="endConstUniformTimeEngine";
	if(!in){
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
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
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
		std::cerr<<" No 'UniformTimeTrain' found...."<<std::endl;
		return 0;
	}
	Vector<std::string> tmm;
	in.getline(tm, 1000, '\n');
	tmm=parse_param(std::string(tm));
	if(tmm.size()!=3)
	{
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::readASCII(istd::fstream)"<<std::endl;
		std::cerr<<" invalid data line in the ASCII file..."<<std::endl;
		std::cerr<<" nothing can be read..."<<std::endl;
		return 0;
	}
	beginT_=atof(tmm[0].c_str());
	endT_=atof(tmm[1].c_str());
	dstep_=atoi(tmm[2].c_str());
	while(temp!=esflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::readASCII(std::fstream)"<<std::endl;
		std::cerr<<" No _End_ for the 'ConstUniformTimeTrain' found...."<<std::endl;
		return 0;
	}
	nele_=nsteps;

	return 1;
}


//reading function reads a BINARY FILE in the following format
//UniformTimeEngine
//<beingT><endT><dstep>
//endUniformTimeEngine
template<int nsteps>
int ConstUniformTimeEngine<nsteps>::read(std::fstream &in)
{
	static std::string sflag="ConstUniformTimeEngine";
	static std::string esflag="endConstUniformTimeEngine";
	if(!in){
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" Bad file...."<<std::endl;
		return 0;
	}
	char tm[1000];
	std::string temp="";
	while(temp!=sflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
		in<<temp<<std::endl;
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" No 'ConstUniformTimeEngine' found...."<<std::endl;
		return 0;
	}
	in.read(&beginT_, sizeof(double));
	in.read(&endT_, sizeof(double));
	in.read(&dstep_, sizeof(int));
	while(temp!=esflag && !in.eof()){
		in.getline(tm, 1000, '\n');
		temp=std::string(tm);
	}
	if(in.eof())
	{
		std::cerr<<std::endl<<"Error: ConstUniformTimeEngine::read(std::fstream)"<<std::endl;
		std::cerr<<" No _End_ for the 'ConstUniformTimeEngine' found...."<<std::endl;
		return 0;
	}
	nele_=nsteps;
	return 1;
}


template<int nsteps>
int ConstUniformTimeEngine<nsteps>::binarySize()
{
	return (22+25)*sizeof(char)+2*sizeof(double)+sizeof(int);
}



END_BL_NAMESPACE
#endif



