
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
	multipowder.h -->
	This has several classes associated with it

	its main purpose is to allow the easy changeing between
	powder averages in the middle of performing a simulation

	however, to properly perform the powder integration using different
	powder sets, at different times, a double (or triple, etc) integral is
	needed in order to avoid powder angle correlation.  THus making the
	simulation quite long, and not only that we have an
	unknown number of loops to perform (it depends on howmany 'powder()' functions
	are called in the main sequence) so there are a few helper classes to handle
	this odd problem...they basically allow one iterator for all the loops


*/

#ifndef __multipowder_cc__
#define __multipowder_cc__ 1

#include "multipowder.h"


/***** THe multi Powder Class ****/

MultiPowder::MultiPowder(Parameters &pset)
{	parse(pset);	}

//this will grab all the
//subsections named 'powder1...powderN' and add them
// to the std::map
void MultiPowder::parse(Parameters &pset)
{
	std::string base="powder";

	int maxFit=pset.getParamI("numpows","",false, 100);
	int numPows=0;
//count the number of params present
	while(pset.addSection(base+itost(numPows+1)) && numPows<=maxFit )
	{	 numPows++;	}

//if there is only one section 'named' powders
// then there is no list, only one powder....
	if(numPows==0 )
	{
		if(pset.addSection(base))
		{
			Parameters tmpP(pset.section(base));
			addPowder(base,tmpP);
			return;
		}
	}

//add a parameter to our master list
	int i=0;
	while(i<numPows)
	{
		Parameters tmPP=Parameters(pset.section(base+itost(i+1)));
		addPowder(base+itost(i+1), tmPP);
		++i;
	}
}


//add a powder from a parameter set
void MultiPowder::addPowder(std::string name,Parameters &PowderPset)
{
	std::string ptype=PowderPset.getParamS("aveType");
	int thst=PowderPset.getParamI("thetaStep", "", false, 0);
	int phst=PowderPset.getParamI("phiStep", "", false, 0);
	int gst=PowderPset.getParamI("gammaStep", "", false, 0);
	powder onPow(ptype, thst, phst,gst);
	PowMapIter j=Powders.find(name);
	if(j==Powders.end()){
		Powders.insert(std::pair<std::string, powder>(name, onPow));
	}else{
		Powders[name]=onPow;
	}
}

//add a powder from an already delcared powder
void MultiPowder::addPowder(std::string name,const powder &powd)
{
	PowMapIter j=Powders.find(name);
	if(j==Powders.end()){
		Powders.insert(std::pair<std::string, powder>(name, powd));
	}else{
		Powders[name]=powd;
	}
}

//finds the powders used from the list inside a pulsedata
// and returns the 'iterator' object
Vector<nonUniqueMapEle<std::string,powder> >
    MultiPowder::findUsedPowders(const Vector<PulseData> &pdata)
{
	Vector<nonUniqueMapEle<std::string,powder> > out;
	std::string last;
	typedef nonUniqueMapEle<std::string,powder> EleMent;

	if(pdata.size()>=1){
		last=pdata[0].usepow;
		out.push_back(EleMent(pdata[0].usepow, Powders[pdata[0].usepow]));
	}
	for(int i=1;i<pdata.size();++i)
	{
		if(last!=pdata[i].usepow)
		{
			out.push_back(EleMent(pdata[i].usepow, Powders[pdata[i].usepow]));
		}
	}
	return out;
}

powder *MultiPowder::getPowder(std::string use)
{
	PowMapIter i;
	if(use=="default")
	{
		i=Powders.begin();
		return &(i->second);
	}

	i=Powders.find(use);
	if(i == Powders.end())
	{
		BLEXCEPTION(std::string(" Powder '")+use+"' does not exsist")
	}
	return &(i->second);
}

ostream &operator<<(ostream &oo,const MultiPowder &out)
{
	MultiPowder::PowMap::const_iterator j=out.Powders.begin();
	while(j !=out.Powders.end())
	{
		oo<<"--Powder Angles--"<<std::endl;
		oo<<"NAME: "<< j->first<<std::endl;
		oo<<"size: "<<j->second.size()<<std::endl;
		++j;
	}
	return oo;
}
#endif


