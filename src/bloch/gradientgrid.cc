
#ifndef _gradientGrid_cc_
#define  _gradientGrid_cc_


#include "bloch/gradientgrid.h"
#include "container/grids/coords.h"
#include "container/Vector/Vector.h"
#include <string>



BEGIN_BL_NAMESPACE


SGradData &SGradData::operator=(const SGradData &rhs)
{
	if(&rhs==this) return *this;
	G_=rhs.G_;
	name_=rhs.name_;
	return (*this);
}


void SGradData::print(std::ostream &oo)
{
	oo<<"Gx: "<<Gx()<<" Gy: "<<Gy()<<" Gz: "<<Gz()<<" on: "<<symbol()<<std::endl;
}

void SGradData::print(){	print(std::cout);	}


std::ostream &operator<<(std::ostream &oo, SGradData &out)
{
	out.print(oo);
	return oo;
}


GradData &GradData::operator=(const GradData &rhs)
{
	if(&rhs==this) return *this;
	data_=rhs.data_;
	return *this;
}


GradData &GradData::operator+(SGradData &rhs)
{
	data_.push_back(rhs);
	return *this;
}

GradData &GradData::operator+=(SGradData &rhs)
{
	data_.push_back(rhs);
	return *this;
}


void GradData::add(coord<> &inG, std::string inname)
{
	data_.push_back(SGradData(inG, inname));
}

void GradData::add(double Gx, double Gy, double Gz, std::string inname)
{
	data_.push_back(SGradData(Gx, Gy, Gz, inname));
}

void GradData::add(SGradData &in)
{
	data_.push_back(in);
}


void GradData::print(std::ostream &oo)
{
	for(int i=0;i<size(); i++)
	{
		oo<<(*this)(i)<<std::endl;
	}
}

void GradData::print()	{ 	print(std::cout);	}


std::ostream &operator<<(std::ostream &oo, GradData &out)
{
	out.print(oo);
	return oo;
}


END_BL_NAMESPACE


#endif


