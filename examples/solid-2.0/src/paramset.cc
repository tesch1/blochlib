/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-26-02
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
	paramset.h -->

	maintains a list of the normal NMR parameters
	it use the 'ParserGlobalVars' as its container for number types
	(so that they can be updated byt any assignments later in
	the sequence parser)

	strings are simply maintained statically...

	to get bools simply use 'getB("name")'
	to get ints simply use 'getI("name")'
	to get doubles simply use 'getD("name")'
	to get strings simply use 'getS("name")'
*/

#ifndef __paramset_cc__
#define __paramset_cc__ 1

#include "paramset.h"

std::string ParamSet::doubleList[10]=
{ "wr", "npts1D", "npts2D", "sw", "sw2",
"maxtstep", "mintstep", "Bfield", "rotor",
"gammaSteps"};

double ParamSet::doubleListdef[10]=
{ 0, 256, 1, 20000,0, 1e-6,
1e-7, 400e6, 0.0,1};

std::string ParamSet::stringList[4]=
{"roeq", "detect", "filesave", "multiInt"};

std::string ParamSet::stringListdef[4]=
{"Iz", "Ip", "soliddata", "no"};

std::map<std::string, std::string> ParamSet::StringList;

void ParamSet::init()
{
	std::map<std::string, double>::iterator i;

	for(int kk=0;kk<9;++kk){
		i=ParserGlobalVars.find(doubleList[kk]);
		if(i==ParserGlobalVars.end()){
			ParserGlobalVars.insert(
				std::pair<std::string, double>(doubleList[kk],doubleListdef[kk]));
		}else{
			ParserGlobalVars[doubleList[kk]]=doubleListdef[kk];
		}
	}

	std::map<std::string, std::string>::iterator j;
	for(int kk=0;kk<4;++kk){
		j=StringList.find(stringList[kk]);
		if(j==StringList.end()){
			StringList.insert(
				std::pair<std::string, std::string>(stringList[kk],stringListdef[kk]));
		}else{
			StringList[stringList[kk]]=stringListdef[kk];
		}
	}
}


void ParamSet::init(const Vector<std::string> &pset)
{
	init();
	Vector<std::string> trimmed=paramStrip(pset), tmPP;
	std::string tmp, tmp2="";

	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;

	for(int i=0;i<trimmed.size();++i){
		tmp=removeWhite(trimmed[i]);
		tmPP=parse_param(tmp,'=');

		if(tmPP.size()>=2){
			//find the name if double
			Diter=ParserGlobalVars.find(tmPP[0]);
			Siter=StringList.find(tmPP[0]);
			//put not found ones in the double
			tmp2=collapsVS(tmPP, 2, tmPP.size());
			if(Diter==ParserGlobalVars.end() && Siter==StringList.end()){
				myParse.parse(tmp2);
				ParserGlobalVars.insert(std::pair<std::string, double>(tmPP[0],myParse() ));
			}else if(Siter!=StringList.end()){ //a string param
				StringList[tmPP[0]]=tmp2;
			}else if(Diter!=ParserGlobalVars.end()){ //a double param
				myParse.parse(tmp2);
				Diter->second=myParse();
			}
		}
	}
}


ParamSet::ParamSet()
{	init();	}

ParamSet::ParamSet(Parameters &pset)
{
	parse(pset);
}

ParamSet::ParamSet(const Vector<std::string> &pset)
{
	parse(pset);
}

void ParamSet::parse(Parameters &pset)
{
	init(pset.section(""));
}

void ParamSet::parse(const Vector<std::string> &pset)
{
	init(pset);
}

bool ParamSet::getB(std::string name)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Siter!=StringList.end()){ //a string param
		if(Siter->second=="no"
			||Siter->second=="n"
			|| Siter->second=="No"
			|| Siter->second=="NO"
			|| Siter->second=="nO"
			|| Siter->second=="N"
		){
			return false;
		}else if(StringList[name]=="0" ){ return false;	}
		else{ 	return true;	}
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		if(Diter->second) return true;
		else return false;
	}else{
		BLEXCEPTION(std::string(" parameter \"")+name+"\" not found")
		return false;
	}
}

int ParamSet::getI(std::string name)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Siter!=StringList.end()){ //a string param
		return std::atoi(Siter->second.c_str());
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		return int(Diter->second);
	}else{
		std::cerr<<"Error: ParamSet::getI"<<std::endl;
		BLEXCEPTION(std::string(" parameter \"")+name+"\" not found")
		return 0;
	}
}

double ParamSet::getD(std::string name)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Siter!=StringList.end()){ //a string param
		return std::atof(Siter->second.c_str());
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		return (Diter->second);
	}else{
		std::cerr<<"Error: ParamSet::getI"<<std::endl;
		BLEXCEPTION(std::string(" parameter \"")+name+"\" not found")
		return 0;
	}
}

std::string ParamSet::getS(std::string name)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Siter!=StringList.end()){ //a string param
		return Siter->second.c_str();
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		return dbtost(Diter->second);
	}else{
		std::cerr<<"Error: ParamSet::getD"<<std::endl;
		BLEXCEPTION(std::string(" parameter \"")+name+"\" not found")
		return "";
	}
}

void ParamSet::set(std::string name, bool in)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Diter==ParserGlobalVars.end() && Siter==StringList.end()){
		ParserGlobalVars.insert(
			std::pair<std::string, double>(name,double(in)));
	}else if(Siter!=StringList.end()){ //a string param
		StringList[name]=itost(int(in));
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		ParserGlobalVars[name]=in;
	}
}

void ParamSet::set(std::string name, int in)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Diter==ParserGlobalVars.end() && Siter==StringList.end()){
		ParserGlobalVars.insert(
			std::pair<std::string, double>(name,double(in)));
	}else if(Siter!=StringList.end()){ //a string param
		StringList[name]=itost(in);
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		ParserGlobalVars[name]=in;
	}
}

void ParamSet::set(std::string name, double in)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Diter==ParserGlobalVars.end() && Siter==StringList.end()){
		ParserGlobalVars.insert(
			std::pair<std::string, double>(name,double(in)));
	}else if(Siter!=StringList.end()){ //a string param
		StringList[name]=dbtost(in);
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		ParserGlobalVars[name]=in;
	}
}

void ParamSet::set(std::string name, std::string in)
{
	std::map<std::string, double>::iterator Diter;
	std::map<std::string, std::string>::iterator Siter;
	Diter=ParserGlobalVars.find(name);
	Siter=StringList.find(name);
	//put not found ones in the double
	if(Diter==ParserGlobalVars.end() && Siter==StringList.end()){
		myParse.parse(in);
		ParserGlobalVars.insert(
			std::pair<std::string, double>(name,myParse()));
	}else if(Siter!=StringList.end()){ //a string param
		StringList[name]=in;
	}else if(Diter!=ParserGlobalVars.end()){ //a double param
		myParse.parse(in);
		ParserGlobalVars[name]=myParse();
	}
}



std::ostream &operator<<(std::ostream &oo,const ParamSet &out)
{
	std::map<std::string, double>::const_iterator Diter;
	std::map<std::string, std::string>::const_iterator Siter;
	Diter=ParserGlobalVars.begin();
	Siter=out.StringList.begin();

	oo<<std::endl<<"---Number Parameters-----"<<std::endl;
	while(Diter!=ParserGlobalVars.end()){
		oo<<Diter->first<<"="<<Diter->second<<std::endl;
		Diter++;
	}
	oo<<std::endl<<"---String Parameters-----"<<std::endl;
	while(Siter!=out.StringList.end()){
		oo<<Siter->first<<"="<<Siter->second<<std::endl;
		Siter++;
	}
	return oo;
}


#endif
