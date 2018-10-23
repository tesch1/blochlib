
/* magfitter.cc ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 5.16.02
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
	magfitter.cc--> methods that aid in setting
	and adjusting the parameters of coils while
	attempting to 'fit' something about the magnetic field

*/


#include "magfitter.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/****
*  MAG FITTER  Data parts
*
****/
const std::string MagFitterData::Alterable[AltLen]={
	"R", "Z", "turns", "theta1", "theta2",
	"amps", "loops",
	"Xcenter", "Ycenter", "Zcenter"
};


MagFitterData::MagFitterData(Parameters &pset)
{	read(pset);	}

MagFitterData::MagFitterData(const Vector<std::string> &pset)
{	read(pset);	}

void MagFitterData::read(const Vector<std::string>  &pset)
{
	Parameters Ppset(pset);
	read(Ppset);
}


//need to check to make sure it is usable
void MagFitterData::setAlter(std::string in)
{
	alterInt=-1;
	alter=in;
	for(int i=0;i<AltLen;++i){
		if(in==Alterable[i]){
			alterInt=i;
			break;
		}
	}
	if(alterInt==-1){
		std::cerr<<std::endl<<"Error: MagFitterData::setAlter"<<std::endl;
		std::cerr<<" the alter type \""<<in<<"\" is NOT valid"<<std::endl;
		std::cerr<<" Allowed Types:: ";
		for(int i=0;i<AltLen;++i)	std::cerr<<" "<<Alterable[i];
		std::cerr<<std::endl;
		exit(1);
	}
}


void MagFitterData::read(Parameters &pset)
{
	use=true;
	alter=pset.getParamS("alter");
	setAlter(alter);
	upperBound=pset.getParamD("upperBound");
	lowerBound=pset.getParamD("lowerBound");
	currentValue=pset.getParamD("start");
	error=pset.getParamD("error");
	pSection=pset.getParamS("coil");
	name=pset.getParamS("name", "", false, pSection+alter);
}

//this would be a valid MINUIT string
//to be used with the function "MNPARS(char *, int err)
// format::
// <var #> '<name>' <start> <error> <lowerbound> <upperbound>
// THe <var #> must be added by 'MagFitter'
std::string MagFitterData::minuitString()
{
	return " '"+name+"' "+dbtost(currentValue)+" "+
			dbtost(error)+" "+
			dbtost(lowerBound)+" "+
			dbtost(upperBound);
}

bool  MagFitterData::setCoilParams(Parameters &pset)
{
	if(pset.section(pSection).empty()) return false; //nothing to set
	coord<> center=pset.getParamCoordD("center", pSection,  ',',false);
	switch(alterInt)
	{
		case Xcenter:
			center.x(currentValue);
			pset.setParam("center", center, pSection);
			return true;
		case Ycenter:
			center.y(currentValue);
			pset.setParam("center", center, pSection);
			return true;
		case Zcenter:
			center.z(currentValue);
			pset.setParam("center", center, pSection);
			return true;
		case R:
		case Z:
		case turns:
		case theta1:
		case theta2:
		case amps:
		case loops:
			pset.setParam(alter, currentValue,pSection);
			return true;
		default:
			std::cerr<<std::endl<<"Error: MagFittingData.setCoilParams()"<<std::endl;
			std::cerr<<std::endl<<" Bad alter Type"<<std::endl;
			return false;
	}
}


void MagFitterData::print()
{	print(std::cout);	}

void MagFitterData::print(std::ostream &oo)
{
	oo.precision(4);
	int wid=oo.width();
	oo.flags(ios::scientific | ios::left);
	oo<<"Fitter Data: name="<<name<<" alter="<<alter<<" upper="<<upperBound
	  <<" lower="<<lowerBound<<" start="<<currentValue
	  <<" error="<<error<<" coil="<<pSection<<endl;
	 oo.width(wid);
}

void MagFitterData::simplePrint(std::ostream &oo)
{
	oo.precision(4);
	int wid=oo.width();
	oo.flags(ios::scientific | ios::left);
	oo<<name<<std::string(20-name.size(), ' ');
	oo<<alter<<std::string(15-alter.size(),' ');
	oo.width(15); oo<<upperBound;
	oo.width(15); oo<<lowerBound;
	oo.width(15); oo<<currentValue;
	oo.width(15); oo<<error;
	oo<<pSection<<std::string(15-pSection.size(), ' ');
	if(!use || error==0) oo<<"NOT used";
	else oo<<"used";
	oo.width(wid);
}

std::ostream &operator<<(std::ostream &oo, MagFitterData &out)
{	out.print(oo); return oo;	}


/****
*  MAG FITTER parts
*
****/

MagFitter::MagFitter(){}

MagFitter::MagFitter(Parameters &Fpset, Parameters &Cpset)
{	read(Fpset, Cpset);	}

MagFitter::MagFitter(const Vector<std::string> &Fpset,const  Vector<std::string> &Cpset)
{	read(Fpset, Cpset);	}


void MagFitter::read(const Vector<std::string> &Fpset,const Vector<std::string> &Cpset)
{
	FitPars=Parameters(Fpset);
	CoilPars=Parameters(Cpset);
	read(FitPars, CoilPars);
}

void MagFitter::read(Parameters &Fpset,Parameters &Cpset)
{
	FitPars=Fpset;
	CoilPars=Cpset;
	subbase=FitPars.getParamS("base","",false, "par");
	int maxFit=FitPars.getParamI("numpars","",false, 1000000);
	numFits=1;
	numUseable=0;
	//while we still have parameters
	while(FitPars.addSection(subbase+itost(numFits)) && numFits<=maxFit ){
	//add a parameter to our master list
		fitparams.push_back(MagFitterData(FitPars.section(subbase+itost(numFits))));

	//now check to see if the 'coil' name is acctually
	// in the 'CoilPars' set and set the 'use' flag appropirately
		if(!CoilPars.addSection(fitparams[numFits-1].pSection)){
			std::cerr<<std::endl<<"Warning: MagFitter::read"<<std::endl;
			std::cerr<<" The 'coil' in not present...not using parameter"<<std::endl;
			fitparams[numFits-1].use=false;
			fitparams[numFits-1].error=0;
			numUseable--;
		}
		++numUseable;
		++numFits;
	}
}

void MagFitter::print(std::ostream &oo, int header)
{
	if(header){
		oo<<"Fitter Parameters"<<std::endl;
		if(fitparams.empty()){
			oo<<" Empty..."<<std::endl;
			return;
		}
		oo<<"- name ------------"
		  <<" type ---------"
		  <<" upper --------"
		  <<" lower --------"
		  <<" start --------"
		  <<" error --------"
		  <<" section ------"
		  <<" used(?) ------"<<std::endl;
		oo<<"---------------------------------------------------------------------------------------------------------------------------"<<std::endl;
	}
	if(fitparams.empty()){
		oo<<" Empty Fitter Parameters..."<<std::endl;
		return;
	}
	for(int i=0;i<fitparams.size();++i)
	{	fitparams[i].simplePrint(oo); oo<<std::endl;	}
}

//this would be a valid MINUIT string
//to be used with the function "MNPARS(char *, int err)"
std::string MagFitter::minuitString(int i)
{
	RunTimeAssert(i<fitparams.size());
	return itost(i+1)+" "+fitparams[i].minuitString();
}

//is the parameter at index 'i' used?
bool MagFitter::used(int i)
{
	if(i<fitparams.size())	return fitparams[i].use;
	return false;
}

//is the parameter name used
bool MagFitter::used(std::string name)
{
	for(int i=0;i<fitparams.size();++i){
		if(name==fitparams[i].name) return fitparams[i].use;
	}
	return false;
}

//sets the Coil Parameters from the MINUIT
// variable array
void MagFitter::setCoilParams(double *vars)
{
	for(int i=0;i<fitparams.size();++i){
		if(fitparams[i].use && fitparams[i].error!=0){
			fitparams[i].currentValue=vars[i];
			fitparams[i].setCoilParams(CoilPars);
		}
	}
}


std::ostream &operator<<(std::ostream &oo, MagFitter &out)
{	out.print(oo); return oo;	}


