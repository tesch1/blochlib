
/* blochParamsDensity.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-09-01
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
 	blochParamsDensity.cc-->the basic storage block for the SINGLE spin..
 */

#ifndef _bloch_params_Density_cc_
#define _bloch_params_Density_cc_ 1



#include "bloch/blochParams.h"
#include "utils/constants.h"
#include "utils/utils.h"
#include "utils/random.h"
#include "bloch/Isotope.h"




BEGIN_BL_NAMESPACE


//A small class that simply contains the 'slew' of parameters
// nessesary to perform a simple bloch simulation

//double BlochParams::Bo_=(4.7);
//double BlochParams::Ho_=(1./permVac *Bo_);
//double BlochParams::temperature_=(300);


BlochParams<BPoptions::Density | BPoptions::HighField, double>
::BlochParams():
	BasicBlochParams<double>(), offset_(0.)
{
	calcMo();
	BasicBlochParams<double>::omega()=BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma();
}

BlochParams<BPoptions::Density | BPoptions::HighField, double>
::BlochParams(std::string spin):
	BasicBlochParams<double>(spin), offset_(0.)
{
	calcMo();
	BasicBlochParams<double>::omega()=BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma();
}


BlochParams<BPoptions::Density | BPoptions::HighField, double>
::BlochParams(const BlochParams &copy):
BasicBlochParams<double>(copy)
{
	offset_=copy.offset_;
}

//Initial Condition Setter..This should be used AFTER !! the 'calcMo()' function call
//this is also the main utility for the 'Bloch' class to set the 'grand' initial
//condition from a list of BlochParams
// 'IC' is the 'Initalondion' class enum
//.. 'HalfFlag' comes from the 'BlochBasic' class where it determins the sign of the
// Mo
void
  BlochParams<BPoptions::Density | BPoptions::HighField, double>::
setInitialCondition(int IC, int HalfFlag)
{
	Random<UniformRandom<> > angleR(-1.0,1.0);
	switch(IC){
		case InitialCondition::RandomDistribution:
			if(BasicBlochParams<double>::Mo() > 0){
				angleR.high(BasicBlochParams<double>::Mo());
			}else{
				angleR.low(BasicBlochParams<double>::Mo());
			}
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomUp:
			angleR.high(abs(BasicBlochParams<double>::Mo()));
			angleR.low(0.0);
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomDown:
			angleR.high(0.0);
			angleR.low(-abs(BasicBlochParams<double>::Mo()));
			BasicBlochParams<double>::Mo()=angleR();
			break;
		case InitialCondition::RandomUpDown:
			angleR.set(0.0, 1.0);
			if(angleR()>0.5){
				BasicBlochParams<double>::Mo()=-BasicBlochParams<double>::Mo();
			}
			break;
		case InitialCondition::AllDown:
			BasicBlochParams<double>::Mo()=-BasicBlochParams<double>::Mo();
			break;
		case InitialCondition::HalfUpHalfDown:
			if(HalfFlag) BasicBlochParams<double>::Mo()=-BasicBlochParams<double>::Mo();
			break;
		case InitialCondition::AllUp: //all are up upon a CalcMo call
		default:
		break;
	}
}


void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::calcMo()
{
	//BasicBlochParams<double>::Mo()=(BasicBlochParams<double>::Bo()*
	//BasicBlochParams<double>::qn()*(BasicBlochParams<double>::qn() + 1.)*BasicBlochParams<double>::moles()*No*hbar*hbar*BasicBlochParams<double>::gamma()*
	//BasicBlochParams<double>::gamma())/3./kb/BasicBlochParams<double>::temperature();
	//gammah*hbar*tanh(hbar*pi*larmor/boltz/temp)*conc*avogd*1e3/2
	BasicBlochParams<double>::Mo()=(BasicBlochParams<double>::gamma()*hbar*
	    tanh(
		  hbar*PI*(BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma()/PI2)/kb/BasicBlochParams<double>::temperature()
		)*BasicBlochParams<double>::moles()*No*1e6/2.0);
}

//this is for use when performing 'dimensionless' integration routines
//the intial magnitization simply becomes 1...or if there are several
//species present in various concentrations
void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::calcMoNorm()
{	BasicBlochParams<double>::Mo()=BasicBlochParams<double>::moles();	}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::Mo(double inMo)
{ BasicBlochParams<double>::Mo()=inMo;	}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::Bo(double inBo){
	BasicBlochParams<double>::Bo()=inBo;
	BasicBlochParams<double>::Ho()=1/permVac*BasicBlochParams<double>::Bo();
	calcMo();
	BasicBlochParams<double>::omega()=BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma();
}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::temperature(double Tempin)
{	BasicBlochParams<double>::temperature()=Tempin;  calcMo();	}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::moles(double inmoles)
{	BasicBlochParams<double>::moles()=inmoles; calcMo();	}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::operator=(const BlochParams<BPoptions::Density | BPoptions::HighField, double>
 &rhs)
{
	if(this==&rhs) return;
	BasicBlochParams<double>::operator=(rhs);
}

BlochParams<BPoptions::Density | BPoptions::HighField, double>
 &BlochParams<BPoptions::Density | BPoptions::HighField, double>
::operator=(const std::string &rhs){
	BasicBlochParams<double>::operator=(rhs);
	calcMo();
	BasicBlochParams<double>::omega()=BasicBlochParams<double>::Bo()*BasicBlochParams<double>::gamma();
	return *this;
}


bool BlochParams<BPoptions::Density | BPoptions::HighField, double>
::operator==(const BlochParams<BPoptions::Density | BPoptions::HighField, double>
 &rhs)
{
	if(this==&rhs) return true;
	if(BasicBlochParams<double>::operator!=(rhs)) return false;
	return true;
}

bool BlochParams<BPoptions::Density | BPoptions::HighField, double>
::operator!=(const BlochParams<BPoptions::Density | BPoptions::HighField, double>
 &rhs)
{
	if(this==&rhs) return false;
	if(BasicBlochParams<double>::operator==(rhs)) return false;
	return true;
}

void BlochParams<BPoptions::Density | BPoptions::HighField, double>
::print(std::ostream &oo) const
{
	oo<<"Density BlochParam: "
	<<"Mo="<<BasicBlochParams<double>::Mo()
	<<", moles="<<BasicBlochParams<double>::moles()<<", Spin="
	  <<BasicBlochParams<double>::symbol()<<", momentum="<<BasicBlochParams<double>::momentum()<<", gamma="<<gamma();

}


int BlochParams<BPoptions::Density | BPoptions::HighField, double>
::binarySize() const
{
	static const char* tms=symbol().c_str();
	return sizeof(char)*2+4*sizeof(double)*4+sizeof(tms)+sizeof(int);
}

/* writing a bloch param goes as follows
	BP <offset> <T1> <T2> <moles> <Mo> <Bo> <temp> <Spin>
*/
bool BlochParams<BPoptions::Density | BPoptions::HighField, double>
::write(std::fstream &oo) const
{
	static const char *bpsym="BP";
	static const char* tms=symbol().c_str();
	int totalsize=binarySize();

	oo.write(&bpsym[0], sizeof(char));
	oo.write(&bpsym[1], sizeof(char));
	oo.write((char *)&totalsize, sizeof(int));
	double tmm=BasicBlochParams<double>::moles();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<double>::Mo();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<double>::Bo();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<double>::temperature();
	oo.write((const char *)&tmm, sizeof(double));
	unsigned int i=0;
	while(i<=symbol().size())
	{
		oo.write(&tms[i], sizeof(char));
		i++;
	}
	return true;
}


int BlochParams<BPoptions::Density | BPoptions::HighField, double>
::read(std::fstream &in)
{
	static char bpflag[2];
	if(!in.read(&bpflag[0], sizeof(char)) || !in.read(&bpflag[1], sizeof(char))) {
		std::cerr<<"BlochParams binary Read Fail"<<std::endl;
		return 0;
	}
	if(std::string(bpflag)!="BP"){
		return 0;
	}
	int size;
	in.read((char *)&size, sizeof(int));
	in.read((char *)&BasicBlochParams<double>::moles(), sizeof(double));
	in.read((char *)&BasicBlochParams<double>::Mo(), sizeof(double));
	in.read((char *)&BasicBlochParams<double>::Bo(), sizeof(double));
	in.read((char *)&BasicBlochParams<double>::temperature(), sizeof(double));

	int remain=size-binarySize();
	char *name_; name_=new char[remain];
	int i=0;
	while(i<remain)
	{
		in.read(&name_[i], sizeof(char));
		i++;
	}
	Isotope::operator=(std::string(name_));
	return size;
}


/************************************************************/
/*  COORD<> OFFSET TYPES *****/
/************************************************************/


BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::BlochParams():
	BasicBlochParams<coord<> >()//, offset_(0.)
{
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
}

BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::BlochParams(std::string spin):
	BasicBlochParams<coord<> >(spin)//, offset_(0.)
{
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
}


BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::BlochParams(const BlochParams &copy):
BasicBlochParams<coord<> >(copy)
{
//	offset_=copy.offset_;
}

//Initial Condition Setter..This should be used AFTER !! the 'calcMo()' function call
//this is also the main utility for the 'Bloch' class to set the 'grand' initial
//condition from a list of BlochParams
// 'IC' is the 'Initalondion' class enum
//.. 'HalfFlag' comes from the 'BlochBasic' class where it determins the sign of the
// Mo
void
  BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >::
setInitialCondition(int IC, int HalfFlag)
{
	Random<UniformRandom<> > angleR(0.0, 1.0);
	double the, phi, gam;
	switch(IC){
		case InitialCondition::RandomDistribution:
			angleR.set(0.0,PI2);
			BasicBlochParams<coord<> >::Mo().Rotate3D(angleR(),angleR(),angleR());
			break;
		case InitialCondition::RandomUp:
			angleR.set(0.0,PI2);
			phi=angleR();
			angleR.set(0.0, PI/2.0);
			the=angleR(); gam=angleR();
			BasicBlochParams<coord<> >::Mo().Rotate3D(phi,the,gam);
			break;
		case InitialCondition::RandomDown:
			angleR.set(0.0, PI2);
			phi=angleR();
			gam=angleR();
			angleR.set(PI/2.0, PI);
			the=angleR();
			BasicBlochParams<coord<> >::Mo().Rotate3D(phi,the,gam);
			break;
		case InitialCondition::RandomUpDown:
			if(angleR()>0.5){
				BasicBlochParams<coord<> >::Mo()=-BasicBlochParams<coord<> >::Mo();
			}
			break;
		case InitialCondition::AllDown:
			BasicBlochParams<coord<> >::Mo().x()=-abs(BasicBlochParams<coord<> >::Mo().x());
			BasicBlochParams<coord<> >::Mo().y()=-abs(BasicBlochParams<coord<> >::Mo().y());
			BasicBlochParams<coord<> >::Mo().z()=-abs(BasicBlochParams<coord<> >::Mo().z());
			break;
		case InitialCondition::HalfUpHalfDown:
			if(HalfFlag){
				BasicBlochParams<coord<> >::Mo().x()=-abs(BasicBlochParams<coord<> >::Mo().x());
				BasicBlochParams<coord<> >::Mo().y()=-abs(BasicBlochParams<coord<> >::Mo().y());
				BasicBlochParams<coord<> >::Mo().z()=-abs(BasicBlochParams<coord<> >::Mo().z());
			}
			break;
		case InitialCondition::AllUp: //all are up upon a CalcMo call
		default:
		break;
	}
}


void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::calcMo()
{
//	BasicBlochParams<coord<> >::Mo()=BasicBlochParams<coord<> >::moles()*No*hbar*hbar*BasicBlochParams<coord<> >::gamma()*
//	BasicBlochParams<coord<> >::gamma()/3./kb/BasicBlochParams<coord<> >::temperature()*BasicBlochParams<coord<> >::Bo()*
//	BasicBlochParams<coord<> >::qn()*(BasicBlochParams<coord<> >::qn() + 1.);
	BasicBlochParams<coord<> >::Mo()=(BasicBlochParams<coord<> >::gamma()*hbar*
	    tanh(
		  hbar*PI*(BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma()/PI2)/kb/BasicBlochParams<coord<> >::temperature()
		)*BasicBlochParams<coord<> >::moles()*No*1e6/2.0);

}

//this is for use when performing 'dimensionless' integration routines
//the intial magnitization simply becomes 1...or if there are several
//species present in various concentrations
void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::calcMoNorm()
{	BasicBlochParams<coord<> >::Mo()=BasicBlochParams<coord<> >::moles();	}

void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::Mo(const coord<> &inMo)
{ BasicBlochParams<coord<> >::Mo()=inMo;	}

void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::Bo(const coord<> &inBo){
	BasicBlochParams<coord<> >::Bo()=inBo;
	BasicBlochParams<coord<> >::Ho()=1/permVac*BasicBlochParams<coord<> >::Bo();
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
}

void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::temperature(double Tempin)
{	BasicBlochParams<coord<> >::temperature()=Tempin;  calcMo();	}

void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::moles(double inmoles)
{	BasicBlochParams<coord<> >::moles()=inmoles; calcMo();	}

/*
void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::offset(const coord<> &inoffset)
{ offset_=inoffset;	}
*/

//Isotope  BlochParams::SpinType() const { return (*this); }
//void  BlochParams::SpinType(const Isotope &iso){  Isotope::operator=(iso);	}


void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::operator=(const BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
 &rhs)
{
	if(this==&rhs) return;
	BasicBlochParams<coord<> >::operator=(rhs);
//	offset_=rhs.offset_;
}

BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
 &BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::operator=(const std::string &rhs){
	BasicBlochParams<coord<> >::operator=(rhs);
	calcMo();
	BasicBlochParams<coord<> >::omega()=BasicBlochParams<coord<> >::Bo()*BasicBlochParams<coord<> >::gamma();
	return *this;
}


bool BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::operator==(const BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
 &rhs)
{
	if(this==&rhs) return true;
	if(BasicBlochParams<coord<> >::operator!=(rhs)) return false;
//	if(offset_.x()!=rhs.offset_.x()) return false;
//	if(offset_.y()!=rhs.offset_.y()) return false;
//	if(offset_.z()!=rhs.offset_.z()) return false;
	return true;
}

bool BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::operator!=(const BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
 &rhs)
{
	if(this==&rhs) return false;
	if(BasicBlochParams<coord<> >::operator==(rhs)) return false;
//	if(offset_.x()==rhs.offset_.x()) return false;
//	if(offset_.y()==rhs.offset_.y()) return false;
//	if(offset_.z()==rhs.offset_.z()) return false;
	return true;
}

void BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::print(std::ostream &oo) const
{
	oo<<"Density BlochParam: "
	//T1="<<BasicBlochParams<coord<> >::T1()
	//  <<", T2="<<BasicBlochParams<coord<> >::T2()<<
	<<"moles="<<BasicBlochParams<coord<> >::moles()<<", Spin="
	  <<BasicBlochParams<coord<> >::symbol()<<", momentum="<<BasicBlochParams<coord<> >::momentum()<<", gamma="<<BasicBlochParams<coord<> >::gamma();

}

int BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::binarySize() const
{
	static const char* tms=symbol().c_str();
	return sizeof(char)*2+8*sizeof(double)*4+sizeof(tms)+sizeof(int);
}

/* writing a bloch param goes as follows
	BP <offset> <T1> <T2> <moles> <Mo> <Bo> <temp> <Spin>
*/
bool BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::write(std::fstream &oo) const
{
	static const char *bpsym="BP";
	static const char* tms=symbol().c_str();
	int totalsize=binarySize();

	oo.write(&bpsym[0], sizeof(char));
	oo.write(&bpsym[1], sizeof(char));
	oo.write((char *)&totalsize, sizeof(int));
	//double tmm=offset_.x();
	//oo.write((const char *)&tmm, sizeof(double));
	//tmm=offset_.y();
	//oo.write((const char *)&tmm, sizeof(double));
	//tmm=offset_.z();
	//oo.write((const char *)&tmm, sizeof(double));
	// tmm=BasicBlochParams<coord<> >::T1();
	//oo.write((const char *)&tmm, sizeof(double));
	//tmm=BasicBlochParams<coord<> >::T2();
	//oo.write((const char *)&tmm, sizeof(double));
	double tmm=BasicBlochParams<coord<> >::moles();
	oo.write((const char *)&tmm, sizeof(double));

	tmm=BasicBlochParams<coord<> >::Mo().x();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<coord<> >::Mo().y();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<coord<> >::Mo().z();
	oo.write((const char *)&tmm, sizeof(double));

	tmm=BasicBlochParams<coord<> >::Bo().x();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<coord<> >::Bo().y();
	oo.write((const char *)&tmm, sizeof(double));
	tmm=BasicBlochParams<coord<> >::Bo().z();
	oo.write((const char *)&tmm, sizeof(double));

	tmm=BasicBlochParams<coord<> >::temperature();
	oo.write((const char *)&tmm, sizeof(double));
	unsigned int i=0;
	while(i<=symbol().size())
	{
		oo.write(&tms[i], sizeof(char));
		i++;
	}
	return true;
}


int BlochParams<BPoptions::Density | BPoptions::HighField, coord<> >
::read(std::fstream &in)
{
	static char bpflag[2];
	if(!in.read(&bpflag[0], sizeof(char)) || !in.read(&bpflag[1], sizeof(char))) {
		std::cerr<<"BlochParams binary Read Fail"<<std::endl;
		return 0;
	}
	if(std::string(bpflag)!="BP"){
		return 0;
	}
	int size;
	in.read((char *)&size, sizeof(int));
	//in.read((char *)&offset_.x(), sizeof(double));
	//in.read((char *)&offset_.y(), sizeof(double));
	//in.read((char *)&offset_.z(), sizeof(double));
	//in.read((char *)&BasicBlochParams<coord<> >::T1(), sizeof(double));
	//in.read((char *)&BasicBlochParams<coord<> >::T2(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::moles(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Mo().x(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Mo().y(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Mo().z(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Bo().x(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Bo().y(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::Bo().z(), sizeof(double));
	in.read((char *)&BasicBlochParams<coord<> >::temperature(), sizeof(double));

	int remain=size-binarySize();
	char *name_; name_=new char[remain];
	int i=0;
	while(i<remain)
	{
		in.read(&name_[i], sizeof(char));
		i++;
	}
	Isotope::operator=(std::string(name_));
	return size;
}

END_BL_NAMESPACE

#endif



