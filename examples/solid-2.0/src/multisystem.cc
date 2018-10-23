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
	multisystem.cc -->
	maintains a list to a bunch of SolidSystems

	the one requirement is that the System does NOT change
	is Hilbert Space dimensions...the things that should change are
	the specific interactions...

	it parses a parameter set like
	-----

	#these are the global (or intial parameters)
	# the numspins and the type cannot change in the sub sections
	# if no subsections are present, the the global is the default
	# spin system returned
	numspins 2
	T 1H 0
	T 1H 1
	D 2000 0 1

	# the first spin section (will overwrite any global interactions)
	spin1{
		C 5000 -2556 0 1
		C -5000 -2556 0 0
		#D 2134 0 1
	}

	# the second spin section (will overwrite any global interactions)
	spin2{
		D 2534 0 1
		J 500 0 1
	}



*/


#ifndef __multisystem_cc__
#define __multisystem_cc__ 1

#include "multisystem.h"


MultiSystem::MultiSystem(Parameters &pset)
{	parse(pset);	}

//cleans the input parameter set for easy and correct reading
// remvove extra sections so it does not read the bad ones
Vector<std::string> MultiSystem::preParse(const Vector<std::string> &pset)
{
	Vector<std::string> out;
	std::string loopArg;
	for(int i=0;i<pset.size();++i)
	{
		loopArg=removeWhite(pset[i]);
		if(loopArg.size()>0 && loopArg[0]!='#' && loopArg[0]!='\n')
		{
			loopArg=removeWhite(pset[i]);
			if(loopArg.find("{")<loopArg.size()){
				while( (i-1)<pset.size()){
					if(loopArg.find("}")<loopArg.size()) break;
					++i;
				}
				if(i==pset.size())
				{
					BLEXCEPTION(" Bad bunch of '}' '{' ")
				}

			}else{
				out.push_back(pset[i].substr(0,pset[i].find("#")));
			}
		}
	}
	return out;
}


//this will grab all the
//subsections named 'powder1...powderN' and add them
// to the std::map
void MultiSystem::parse(Parameters &pset)
{
	std::string base="spin";
	global=preParse(pset.section(""));
//if there is only one section 'named' powders
// then there is no list, only one powder....

	int maxFit=pset.getParamI("numsys","",false, 100);
	int numPows=0;
//count the number of params present
	while(pset.addSection(base+itost(numPows+1)) && numPows<=maxFit )
	{	numPows++;		}

	if(numPows==0)
	{
		Vector<std::string> empt(1," ");
		Parameters pp(empt);
		addSystem(base+itost(1),pp);
		return;
	}

//add a parameter to our master list
	int i=0;
	while(i<numPows)
	{
		Parameters tmPP=pset.section(base+itost(i+1));
		addSystem(base+itost(i+1), tmPP);
		++i;
	}
}
//add a powder to the list from a parameter set list
void MultiSystem::addSystem(std::string name,Parameters &pset)
{
	Vector<std::string> spV=global, curChunk=pset.section("");

	for(int i=0;i<curChunk.size();++i){
		spV.push_back(curChunk[i]);
	}
	SolidSys onSys;
	onSys.read(spV);
	SysMapIter j=Systems.find(name);
	if(j==Systems.end()){
		Systems.insert(std::pair<std::string, SolidSys>(name, onSys));
	}else{
		Systems[name]=onSys;
	}
}


//add a powder to the list from a SolidSys
void MultiSystem::addSystem(std::string name,const SolidSys &onSys)
{
	SysMapIter j=Systems.find(name);
	if(j==Systems.end()){
		Systems.insert(std::pair<std::string, SolidSys>(name, onSys));
	}else{
		Systems[name]=onSys;
	}
}

//returns the pointer to the powder in the list
SolidSys *MultiSystem::getSystem(std::string use)
{
	SysMapIter i;
	if(use=="default")
	{
		i=Systems.begin();
		return &(i->second);
	}
	i=Systems.find(use);
	if(i == Systems.end())
	{
		std::cerr<<"Error: MultiSystem::getSystem"<<std::endl;
		std::cerr<<" SoildSys '"<<use<<"' does not exsist"<<std::endl;
		BLEXCEPTION(std::string("SoildSys '")+use+"' does not exsist")
	}
	return &(i->second);
}

//sets the Bfield in all the systems
void MultiSystem::setBfield(double BF)
{
	SysMapIter i=Systems.begin();
	while(i!=Systems.end())
	{
		i->second.setBfield(BF);
		i++;
	}
}


std::ostream &operator<<(std::ostream &oo,const MultiSystem &out)
{
	MultiSystem::SysMap::const_iterator i=out.Systems.begin();
	while(i != out.Systems.end()){
		oo<<"--Spin System--"<<std::endl;
		oo<<"NAME: "<< i->first<<std::endl;
		SolidSys *tmm=const_cast<SolidSys *>(&(i->second));
		oo<<(*tmm)<<std::endl;
		++i;
	}
	return oo;
}

#endif
