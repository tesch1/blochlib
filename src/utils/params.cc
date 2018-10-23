/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-4-01
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

/********************************

	Parameter grabber class....can suck out pararmeters from a file...
	that exsists in specific chunks

	this class has various options for turning off and on certin input
	file format bits...explained in code

	an input file can be like so... where 'moo' is a section...

	moo{
		pars1 34
		par2 54
	}

	moo2{
		par1 567
		par 34
	}


	or like so (with no section and '=' parameter separators)

	par=90
	par=78

	or combos...

********************************/


#ifndef _params_cc_
#define _params_cc_ 1


#include <list>
#include <string>
#include <map>
#include <fstream>
#include "utils/params.h"
#include "utils/utils.h"

BEGIN_BL_NAMESPACE



//assumes that 'in' is an entire 'file'
Parameters::Parameters(const Vector<std::string> &in,
	std::string open, std::string cl, char psep, bool interac)
{

	MaxLineLength=10000;
	secopen=open;
	secclose=cl;
	interactive=interac;
	paramsep=psep;
	filechop.insert(std::pair<std::string,Vector<std::string> >("", in));

}

//assumes the input vector is a section 'flag'
Parameters::Parameters(std::string flag, Vector<std::string> &in,
	std::string open, std::string cl, char psep, bool interac)
{
	MaxLineLength=10000;
	secopen=open;
	secclose=cl;
	interactive=interac;
	paramsep=psep;
	filechop.insert(std::pair<std::string,Vector<std::string> >("", in));
	addSection(flag);
}


//will read the file according to the 'rules'
Parameters::Parameters(std::string fname,
	std::string open, std::string cl, char psep, bool interac)
{
	MaxLineLength=10000;
	secopen=open;
	secclose=cl;
	paramsep=psep;
	interactive=interac;
	filechop.insert(std::pair<std::string,Vector<std::string> >("", SuckFile(fname)));
}


//will read the file according to the 'rules'
Parameters::Parameters(std::ifstream &readme,
	std::string open, std::string cl, char psep, bool interac)
{
	MaxLineLength=10000;
	secopen=open;
	secclose=cl;
	paramsep=psep;
	interactive=interac;
	filechop.insert(std::pair<std::string,Vector<std::string> >("", SuckFile(readme)));
}

//will read the file
void Parameters::operator=(std::ifstream &in)
{
	MaxLineLength=10000;
	secopen="{";
	secclose="}";
	paramsep='\0';
	interactive=true;
	filechop.insert(std::pair<std::string,Vector<std::string> >("", SuckFile(in)));
}

//will read the file
void Parameters::operator=(const Parameters &in)
{
	if(this==&in) return;
	MaxLineLength=in.MaxLineLength;
	secopen=in.secopen;
	secclose=in.secclose;
	paramsep=in.paramsep;
	interactive=in.interactive;
	filechop=in.filechop;
}



Vector<std::string> Parameters::SuckFile(std::ifstream &infile)
{
	int initmax=5000;
	Vector<std::string> tmp(initmax,"");
	char liner[10000];
	int i=0;
	while(infile.peek()!=EOF)
	{
		infile.getline(liner, MaxLineLength, '\n');
		if(i>initmax){	tmp.push_back(std::string(liner));	}
		else{	tmp[i]=std::string(liner);	}
		i++;
	}
	tmp.resizeAndPreserve(i);
	return tmp;
}

Vector<std::string> Parameters::SuckFile(std::string &infile)
{
	std::ifstream in(infile.c_str());
	if(in.fail())
	{
		BLEXCEPTION("Error: cannot load input file...")
	}
	return SuckFile(in);
}

void Parameters::read(std::string &infile){
	filechop.insert(std::pair<std::string,Vector<std::string> >("", SuckFile(infile)));
}

void Parameters::read(std::ifstream &infile){
	filechop.insert(std::pair<std::string,Vector<std::string> >("", SuckFile(infile)));
}

//gets the insides of a 'secopen' 'seclose' section
Vector<std::string> Parameters::GetInside(std::string flag, Vector<std::string> &from, std::string seco, std::string secc)
{
	int i, gotr=-1, gotl=-1, gotsp=-1;
	int len=from.size(), Lct=0;
	Vector<std::string> out;
	std::string holdd, nowhite;
	//find the flag, find the '{' line, and find the '}' line
	for(i=0;i<len;i++){
		nowhite=removeWhite(from[i]);
		nowhite=nowhite.substr(0, nowhite.find("#"));
		if(gotsp==-1  && nowhite.substr(0, flag.size())==flag){
			Vector<std::string> st=parse_param(from[i], paramsep);
			if(st.size()>=1){
				if(st[0].find(flag)<st[0].size()){ gotsp=i+1;	}
			}
		}
		if(from[i].find(seco)<from[i].size() && gotsp!=-1){
			if(Lct==0){gotl=i; }
			++Lct;
		}
		if(from[i].find(secc)<from[i].size() && gotsp!=-1 && gotl>=0)
		{
			--Lct;
			if(Lct==0){
				gotr=i;
				break;
			}
		}
	}
	//check that we got everything
	if(gotsp==-1 || gotr==-1 || gotl==-1){
		/*cerr<<endl<<"Error Parameters::GetInside()"<<endl;
		cerr<<" could not find a valid set"<<endl;
		cerr<<" syntax should be"<<endl;
		cerr<<"--> "<<flag<<seco<<endl;
		cerr<<"--> ..."<<endl<<"--> stuff"<<endl<<"--> ..."<<endl;
		cerr<<"--> "<<secc<<endl;*/
		return Vector<std::string>();
	}
	//make sure we get everything past the initi l '{'
	if(from[gotl].find(seco)<(from[gotl].size()-1) && gotsp!=-1 && from[gotsp-1].find(flag)){
		if(from[gotl].find(seco)-from[gotl].size() > 0)
			out.push_back(from[gotl].substr(from[gotl].find(seco)+1, from[gotl].size()));
	}

	//get the stuff in the middle
	for(i=gotl+1;i<gotr;i++){
		out.push_back(from[i]);
	}

	//get everything before the closing '}'
	if(from[gotr].find(secc)>0){
		if(from[gotr].find(secc)-1 > 0)
			out.push_back(from[gotr].substr(0, from[gotr].find(secc)-1));
	}
	return out;
}


//get a double parameter 'par' from the section 'sec'
double Parameters::getParamD(std::string par,std::string sec, bool required, double def)
{
	Vector<std::string> tmp, vsec=section(sec);
	for(int i=0;i<vsec.size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec[i]);	}
		else{ tmp=parse_param(vsec[i], paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par){	return std::atof(tmp[1].c_str());	}
		}
	}
	if(interactive && required){
		double out;
		std::cout<<std::endl<<"Parameter: "<<par<<" not found in section "<<sec<<" Please Enter: ";
		std::cout.flush();
		std::cin>>out;
		return out;
	}else{
		return def;
	}
}

//get a int parameter 'par' from the section 'sec'
int Parameters::getParamI( std::string par,std::string sec, bool required, int def)
{
	Vector<std::string> tmp, vsec=section(sec);
	for(int i=0;i<vsec.size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec[i]);	}
		else{ tmp=parse_param(vsec[i], paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par){	return std::atoi(tmp[1].c_str());	}
		}
	}
	if(interactive && required){
		int out;
		std::cout<<std::endl<<"Parameter: "<<par<<" not found in section "<<sec<<" Please Enter: ";
		std::cout.flush();
		std::cin>>out;
		return out;
	}else{
		return def;
	}

}

//get a std::string parameter 'par' from the section 'sec'
std::string Parameters::getParamS( std::string par,std::string sec, bool required, std::string def)
{
	Vector<std::string> tmp, vsec=section(sec);
	for(int i=0;i<vsec.size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec[i]);	}
		else{ tmp=parse_param(vsec[i], paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par){	return collapsVS(tmp, 2, tmp.size());	}
		}
	}
	if(interactive && required){
		std::string out;
		std::cout<<std::endl<<"Parameter: '"<<par<<"' not found in section '"<<sec<<"' Please Enter: ";
		std::cout.flush();
		std::cin>>out;
		return out;
	}else{
		return def;
	}

}

//get a char parameter 'par' from the section 'sec'
char Parameters::getParamC(std::string par, std::string sec, bool required, char def)
{
	Vector<std::string> tmp, vsec=section(sec);
	for(int i=0;i<vsec.size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec[i]);	}
		else{ tmp=parse_param(vsec[i], paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par){	return tmp[1][0];	}
		}
	}
	if(interactive && required){
		char out;
		std::cout<<std::endl<<"Parameter: "<<par<<" not found in section "<<sec<<" Please Enter: ";
		std::cout.flush();
		std::cin>>out;
		return out;
	}else{
		return def;
	}

}

//get a coord<double, 3> parameter 'par' from the section 'sec'
coord<> Parameters::getParamCoordD(std::string par, std::string sec, char sep, bool required, coord<> def)
{
	std::string maxx=getParamS(par, sec, required,"");
	if(maxx==""){
		return def;
	}
	Vector<std::string> maxs=parse_param(maxx, sep);
	coord<> mymax(0);
	if(maxs.size()==1){
		mymax(std::atof(maxs[0].c_str()), std::atof(maxs[0].c_str()),std::atof(maxs[0].c_str()));
	}else if(maxs.size()==3){
		mymax(std::atof(maxs[0].c_str()), std::atof(maxs[1].c_str()),std::atof(maxs[2].c_str()));
	}else{
		std::string mess=" Coords inputs must be written as ";
		mess+="\n \""+par+paramsep+"<num>"+sep+"<num>"+sep+"<num>\" or";
		mess+="\n \""+par+paramsep+"<num>\"";
		BLEXCEPTION(mess)
	}
	return mymax;
}

//get a coord<int, 3> parameter 'par' from the section 'sec'
coord<int,3> Parameters::getParamCoordI(std::string par, std::string sec, char sep, bool required, coord<int,3> def)
{
	std::string maxx=getParamS(par, sec, required,"");
	if(maxx==""){
		return def;
	}
	Vector<std::string> maxs=parse_param(maxx, sep);
	coord<int, 3> mymax;
	if(maxs.size()==1){
		mymax(std::atoi(maxs[0].c_str()), std::atoi(maxs[0].c_str()),std::atoi(maxs[0].c_str()));
	}else if(maxs.size()==3){
		mymax(std::atoi(maxs[0].c_str()), std::atoi(maxs[1].c_str()),std::atoi(maxs[2].c_str()));
	}else{
		std::string mess=" Coords inputs must be written as ";
		mess+="\n \""+par+paramsep+"<int>"+sep+"<int>"+sep+"<int>\" or";
		mess+="\n \""+par+paramsep+"<int>\"";
		BLEXCEPTION(mess)
	}
	return mymax;
}

//get a Vector<double> parameter 'par' from the section 'sec'
Vector<double> Parameters::getParamVectorD(std::string par, std::string sec, char sep, bool required, Vector<double> def)
{
	std::string maxx=getParamS(par, sec, required,"");
	if(maxx==""){
		return def;
	}
	Vector<std::string> maxs=parse_param(maxx, sep);
	Vector<double> mymax;
	if(maxs.size()==1){
		mymax.push_back(std::atof(maxs[0].c_str()));
	}else{
		for(int i=0;i<maxs.size();++i){
			mymax.push_back(std::atof(maxs[i].c_str()));
		}
	}
	return mymax;
}

//get a Vector<int> parameter 'par' from the section 'sec'
Vector<int> Parameters::getParamVectorI(std::string par, std::string sec, char sep, bool required, Vector<int> def)
{
	std::string maxx=getParamS(par, sec, required,"");
	if(maxx==""){
		return def;
	}
	Vector<std::string> maxs=parse_param(maxx, sep);
	Vector<int> mymax;
	if(maxs.size()==1){
		mymax.push_back(std::atoi(maxs[0].c_str()));
	}else{
		for(int i=0;i<maxs.size();++i){
			mymax.push_back(std::atoi(maxs[i].c_str()));
		}
	}
	return mymax;
}


//get a Vector<coord<> > parameter 'par' from the section 'sec'
// syntax must be "num, num, num | num, num, num |...
Vector<coord<> > Parameters::getParamVectorCoordD(std::string par, std::string sec, char sep, bool required, Vector<coord<> > def)
{
	std::string maxx=getParamS(par, sec, required,"");
	if(maxx==""){
		return def;
	}
	Vector<std::string> maxs=parse_param(maxx, '|');
	Vector<coord<> > mymax;
	for(int i=0;i<maxs.size();++i){
		Vector<std::string> inter=parse_param(maxs[i], sep);
		if(maxs.size()==1){
			mymax.push_back(coord<>(std::atof(inter[0].c_str()),std::atof(inter[0].c_str()),std::atof(inter[0].c_str())));
		}else if(maxs.size()==3){
			mymax.push_back(coord<>(std::atof(inter[0].c_str()), std::atof(inter[1].c_str()),std::atof(inter[2].c_str())));
		}else{
			std::cerr<<std::endl<<"Error: Parameters::getParamVectorCoordD()"<<std::endl;
			std::string mess=" VectorCoords inputs must be written as  ";
			mess+="\n \""+par+paramsep+"<num>"+sep+"<num>"+sep+"<num>|<num>"+sep+"<num>"+sep+"<num>|...\" or";
			mess+="\n \""+par+paramsep+"<num> | <num> |...\"";
			BLEXCEPTION(mess)
		}
	}
	return mymax;
}


//returns the section 'sec' from the std::map if it exsits.
Vector<std::string> &Parameters::section(std::string sec)
{
	std::map<std::string,Vector<std::string> >::iterator i;
	i = filechop.find(sec);
	if(i==filechop.end())
	{
		std::cerr<<std::endl<<"Warning:: Pararmeters::section(sec)"<<std::endl;
		std::cerr<<" No File has been loaded yet..."<<std::endl;
		std::cerr<<" cannot find add any secs yet..."<<std::endl;
		return ZeroType<Vector<std::string> >::zero();
	}
	return i->second;
}

//Adds the section 'sec' from the std::map if it exsits.
//return false if not found...
bool Parameters::addSection(std::string sec)
{
	std::map<std::string,Vector<std::string> >::iterator i;
	//std::map<std::string,Vector<std::string> >::iterator j;
	i = filechop.find("");
	if(i==filechop.end())
	{
		std::cerr<<std::endl<<"Warning:: Pararmeters::addSection(sec)"<<std::endl;
		std::cerr<<" No File has been loaded yet..."<<std::endl;
		std::cerr<<" cannot add any sections yet..."<<std::endl;
		return false;
	}else if(sec=="" && i!=filechop.end()){
		std::cerr<<std::endl<<"Warning:: Pararmeters::addSection(sec)"<<std::endl;
		std::cerr<<" cannot 're-add' the main file section.."<<std::endl;
		std::cerr<<" performing no actions.."<<std::endl;
		return true;
	}

	//j = filechop.find(sec);
	//if(j != filechop.end()){
	//	std::cerr<<std::endl<<"Warning:: Pararmeters::addSection(sec)"<<std::endl;
	//	std::cerr<<" section "<<sec<<" already added...overwriting old one"<<std::endl;
	//}
/*	std::map<std::string,Vector<std::string> >::iterator j=filechop.find(sec);
	if(j!=filechop.end() && sec!=""){
		//filechop.erase(j);
		filechop[j->first]=GetInside(sec, i->second, secopen, secclose) ;
	}
	else{*/
	if(i==filechop.end()){
		filechop.insert(
				std::pair<std::string, Vector<std::string> >(
					sec, GetInside(sec, i->second, secopen, secclose) ) );
	}else{
		filechop[sec]=GetInside(sec, i->second, secopen, secclose);
	}
	i = filechop.find(sec);
	if(i==filechop.end() || (i->second).empty()) return false;
	return true;
}



/****** Parmeter Sett *****/
//this gives one the ability to set (i.e. change or create) a parameter in a 'seciton'
void Parameters::setParam(std::string par,  int toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				tmp[1]+=itost(toset);
				(*vsec)[i]=collapsVS(tmp);
				return;
			}
		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;
	in+=itost(toset);
	vsec->push_back(in);
}

void Parameters::setParam(std::string par,  double toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				tmp[1]+=dbtost(toset);
				(*vsec)[i]=collapsVS(tmp);
				return;
			}

		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;
	in+=dbtost(toset);
	vsec->push_back(in);

}
void Parameters::setParam(std::string par,  std::string toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				tmp[1]+=toset;
				(*vsec)[i]=collapsVS(tmp);
				return;
			}

		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;
	in+=toset;
	vsec->push_back(in);
}

void Parameters::setParam(std::string par,  char toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				tmp[1]+=std::string(&toset, 1);
				(*vsec)[i]=collapsVS(tmp);
				return;
			}

		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;
	in+=std::string(&toset,1);
	vsec->push_back(in);
}

void Parameters::setParam(std::string par,  Vector<double> toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				for(int j=0;j<toset.size();++j){
					tmp[1]+=dbtost(toset[j]);
					if(j!=toset.size()-1) tmp[1]+=",";
				}
				(*vsec)[i]=collapsVS(tmp);
				return;
			}

		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;

	for(int j=0;j<toset.size();++j){
		in+=dbtost(toset[j]);
		if(j!=toset.size()-1) in+=",";
	}
	vsec->push_back(in);

}
void Parameters::setParam(std::string par,  Vector<int> toset, std::string sec)
{
	Vector<std::string> tmp;
	Vector<std::string> *vsec=&section(sec);
	for(int i=0;i<vsec->size();i++)
	{
		if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
		else{ tmp=parse_param(vsec->get(i), paramsep); }
		if(tmp.size()>0){
			if(tmp[0]==par)
			{
				if(paramsep=='\0')	tmp[1]=" ";
				else tmp[1]=paramsep;
				for(int j=0;j<toset.size();++j){
					tmp[1]+=itost(toset[j]);
					if(j!=toset.size()-1) tmp[1]+=",";
				}
				(*vsec)[i]=collapsVS(tmp);
				return;
			}

		}
	}
	//the param as not in the list..so add it
	std::string in=par;
	if(paramsep=='\0')	in+=" ";
	else in+=paramsep;

	for(int j=0;j<toset.size();++j){
		in+=itost(toset[j]);
		if(j!=toset.size()-1) in+=",";
	}
	vsec->push_back(in);
}

//writes the parameters to a file
bool Parameters::print(std::ostream &out)
{
	std::map<std::string,Vector<std::string> >::iterator i;
	i=filechop.begin();
	while(i!=filechop.end())
	{
		if(i->first!=""){
			out<<"#Sub Sections"<<std::endl;
			out<<i->first<<secopen<<std::endl;
		}else{
			out<<"#Main Section"<<std::endl;
		}
		Vector<std::string> tmv=paramStrip(i->second);
		if(!tmv.empty()){
			for(int j=0;j<tmv.size();++j){
				out<<"\t"<<tmv.get(j)<<std::endl;
			}
		}
		if(i->first!=""){
			out<<secclose<<std::endl;
		}
		i++;
	}
	return true;
}

//writes the parameters to a file
bool Parameters::print(std::string out)
{
	std::ofstream outf(out.c_str());
	if(outf.fail())
	{
		std::cerr<<std::endl<<"Error: Parameters::print()"<<std::endl;
		std::cerr<<" could not open file to print"<<std::endl;
		return false;
	}
	return print(outf);
}

bool Parameters::print(const char *out)
{	return print(std::string(out));	}


std::ostream &operator<<(std::ostream &oo, Parameters &out)
{
	out.print(oo);
	return oo;
}

END_BL_NAMESPACE



#endif


