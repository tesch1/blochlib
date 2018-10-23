

#ifndef _pulse_cc_
#define _pulse_cc_

#include "bloch/pulse.h"
#include <string.h>
#include "container/Vector/Vector.h"
#include <iostream>


BEGIN_BL_NAMESPACE


//Methods for Pulse Sclasses

SPulse &SPulse::operator=(const SPulse &rhs)
{
	if(&rhs==this) return *this;
	name_=rhs.name_;
	phase_=rhs.phase_;
	sph_=rhs.sph_;
	cph_=rhs.cph_;
	amp_=rhs.amp_;
	off_=rhs.off_;
	tb_=rhs.tb_;
	te_=rhs.te_;
	applyall_=rhs.applyall_;
	return *this;
}

double SPulse::timeForAngle(double angle)
{
	double tpulse=(angle)/(amp_);
	if(tpulse<=0)	tpulse=1e-30;
	return tpulse;
}

void SPulse::print(std::ostream &oo)
{
	oo<<"Pulse on "<<name_<<" B1: "<<amp_<<" phase: "<<phase_;
	if(off_!=0) oo<<" offset: "<<off_;
	if(!applyall_) oo<<" Starttime: "<<tb_<<" endTime: "<<te_;
	oo<<std::endl;
}

void SPulse::print()
{
	print(std::cout);
}

std::ostream &operator<<(std::ostream &oo, SPulse &out)
{
	out.print(oo);
	return oo;
}

//Methods for Pulse class

Pulse::Pulse(SPulse &in)
{
	data_.resize(1, in);
}

Pulse::Pulse(Vector<SPulse> &in)
{
	data_=in;
}

Pulse::Pulse(std::string nm, double amp, double phase, double offset, double tbegin, double tend)
{
	SPulse tm(nm, amp,phase,offset, tbegin, tend);
	data_.resize(1, tm);
}

SPulse &Pulse::operator()(int j)
{
	return data_.data()[j];
}

SPulse &Pulse::operator[](int j)
{
	return data_.data()[j];
}


double Pulse::timeForAngle(double angle,const std::string &nm)
{
	static int i=Exists(nm);
	if(i!=-1) return data_(i).timeForAngle(angle);
	return 1.e-30;
}

//Pulse Assignments
void Pulse::operator=( SPulse &rhs)
{
	data_.resize(1, rhs);
}

void Pulse::operator=( Pulse &rhs)
{
	if(&rhs==this) return;
	data_=rhs.data_;
}

//Symbol retrival
std::string Pulse::symbol(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).symbol();
	return "";
}

 void Pulse::setSymbol(int tes, std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).SetSymbol(nm);
}


//Amplitudes....
double Pulse::amplitude(const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).amplitude();
	return 0;
}

double Pulse::amplitude(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).amplitude();
	return 0;
}

double Pulse::amplitude(double time, int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).amplitude(time);
	return 0;
}

double Pulse::amplitude(double time, std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).amplitude(time);
	return 0;
}

void Pulse::setAmplitude(double amp, std::string nm)
{
	int i=Exists(nm);
	if(i!=-1)  data_(i).amplitude(amp);
}

void Pulse::setAmplitude(double amp, int nm)
{
	int i=Exists(nm);
	if(i!=-1)  data_(i).amplitude(amp);
}

//PHASES
double Pulse::phase(const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).phase();
	return 0;
}

double Pulse::phase(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).phase();
	return 0;
}

double Pulse::phase(double time, int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).phase(time);
	return 0;
}

double Pulse::phase(double time, const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).phase(time);
	return 0;
}

void Pulse::setPhase(double amp,const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1)  data_(i).phase(amp);
}

void Pulse::setPhase(double amp, int nm)
{
	int i=Exists(nm);
	if(i!=-1)  data_(i).phase(amp);
}

//OFFSETS
double Pulse::offset(const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).offset();
	return 0;
}

double Pulse::offset(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).offset();
	return 0;
}

double Pulse::offset(double time, int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).offset(time);
	return 0;
}

double Pulse::offset(double time, const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).offset(time);
	return 0;
}

void Pulse::setOffset(double off, int nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setOffset(off);
}

void Pulse::setOffset(double off, const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setOffset(off);
}

//Time functions......................
void Pulse::setTime(double tb, double te, const std::string &nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setTime(tb, te);
}

void Pulse::setTime(double tb, double te, int nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setTime(tb, te);
}

void Pulse::setBeginTime(double tb, std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setBeginTime(tb);
}

double Pulse::beginTime(std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).beginTime();
	return 0.;
}

double Pulse::beginTime(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).beginTime();
	return 0.;
}

void Pulse::setBeginTime(double tb, int nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setBeginTime(tb);
}

void Pulse::setEndTime(double te, std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setEndTime(te);
}

void Pulse::setEndTime(double te, int nm)
{
	int i=Exists(nm);
	if(i!=-1) data_(i).setEndTime(te);
}

double Pulse::endTime(std::string nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).endTime();
	return 0.;
}

double Pulse::endTime(int nm)
{
	int i=Exists(nm);
	if(i!=-1) return data_(i).endTime();
	return 0.;
}


//adding elements to the pulse list
Pulse &Pulse::operator+(const SPulse &in)
{
	data_.push_back(in);
	return *this;
}

Pulse &Pulse::operator+=(const SPulse &in)
{
	data_.push_back(in);
	return *this;
}



Pulse &Pulse::operator+=(const Pulse &in)
{
	if(&in==this) return *this;
	for(int i=0;i<in.size();i++){
		data_.push_back(in.data_.data()[i]);
	}
	return *this;
}

Pulse &Pulse::operator+(const Pulse &in)
{
	if(&in==this) return *this;
	for(int i=0;i<in.size();i++){
		data_.push_back(in.data_.data()[i]);
	}
	return *this;
}

void Pulse::add(SPulse in)
{
	data_.push_back(in);
}


void Pulse::add(std::string nm, double amp, double phase,double offset, double tbegin, double tend)
{
	static SPulse mm(nm,amp, phase,offset,tbegin,tend);
	add(mm);
}

//I/O for pulse
void Pulse::print(std::ostream &oo)
{
	for(int i=0;i<size();i++)
	{
		oo<<data_(i);
	}
}

void Pulse::print()
{
	print(std::cout);
}


std::ostream &operator<<(std::ostream &oo, Pulse &out)
{
	out.print(oo);
	return oo;
}

END_BL_NAMESPACE



#endif

