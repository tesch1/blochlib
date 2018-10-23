

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-25-02
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
 	sequenceparse.h --> this takes in a 'vector string' (usually from
 	a Parameter set, and will create a 'PulseData' vector  from the input

 	the input allows for variable setting, and looping over parameters

 	variables are simply declared via a simple syntax

 	A=B (a gets set from B) where B can be some expression

 	loops are performed via the

 	loop i=1:40

 	...do things

 	end

 	syntax where the i=1:40 means loop from 1 to 40

 	to set a Pulse the function 'pulse' is used and can be

 	pulse(time) used if amplitude is set already via the 'amplitude' function
 	pulse(time, phase)
 	pulse(time, phase, amplitude)
 	pulse(time, phase, amplitude, offset)

 	the amplitude, is in Hz, phase in degreed, and offset in Hz

 	delays are set via

 	delay(time) --> NONE incremented delay
 	delay(time, dt) --> incremented delay (a 2D) type step

 	rotor anlges are set via

 	rotor(angle) --> angles in degrees

 	cycler traces are set via

 	cycler(SpinOp)

	To change the spin system use

	spinsys(sysname)

	To change the powder system use

	powder(powsec)
*/

#ifndef __sequenceparse_cc__
#define __sequenceparse_cc__ 1

#include "sequenceparse.h"
#include "pulsedata.h"


/**************** Sequence Parser bits *************/
void SequenceParse::addGs()
{
	//myParse.addGlobalVar("rotor", acos(1.0/sqrt(3.0)));
	//myParse.addGlobalVar("wr",0);
	curRot=acos(1.0/sqrt(3.0));
	curAmp=0;
	curOff=0;
	curSys="default";
	curPow="default";
	curOn="1H";
	curCycler="Ie";
	rWhite="";
	currentTime=0;
}

SequenceParse::SequenceParse()
{
	addGs();
}

SequenceParse::SequenceParse(bool gv)
 :addGvars(gv)
{
	addGs();
}


SequenceParse::SequenceParse(Parameters &pset,bool addGv)
 :addGvars(addGv)
{
	addGs();
	pulses.resize(0);
	data=pset.section("");
	ScriptParse::parse(data);
}

SequenceParse::SequenceParse(const Vector<std::string> &pset,bool addGv)
 :addGvars(addGv)
{
	addGs();
	pulses.resize(0);
	data=pset;
	ScriptParse::parse(data);
}

SequenceParse::SequenceParse(const SequenceParse &rhs):
	ScriptParse(rhs), addGvars(rhs.addGvars), pulses(rhs.pulses)
{

	if(this==&rhs) return;
	*this=rhs;
}

void SequenceParse::operator=(const SequenceParse &rhs)
{
	if(this==&rhs) return;
	ScriptParse::operator=(rhs);
	addGvars=rhs.addGvars;
	sequence=rhs.sequence;
	pulses=rhs.pulses;

}

//this adds the proper line to the pulseData list
//proper line to the sequence vector
void SequenceParse::addLine(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps, tmP;
	std::string tmS=removeWhite(inSqe);
	//a special change for a 'rotor' assignment
	if(tmS.substr(0,6)=="rotor="){
		tmP=parse_param(tmS, '=');
		sequence.push_back(std::string("rotor(")+collapsVS(tmP, 2, tmP.size())+")");
	}else{
		sequence.push_back(tmS);
	}
}

void SequenceParse::addPulses(std::string curLin)
{
	int curAdd=0;
	bool added=false;

	curRot=myParse.getVar("rotor")*DEG2RAD;

	PulseData tmPD;
	tmPD.myParse=&myParse;
	rWhite=removeWhite(curLin);
	//these are the global parmeter functions

	if(rWhite.find("amplitude(")<rWhite.size()){
		inside=getInside(rWhite);
		myParse.parse(inside);
		curAmp=myParse();

	}else if(rWhite.find("offset(")<rWhite.size()){
		inside=getInside(rWhite);
		myParse.parse(inside);
		curOff=myParse();

	}else if(rWhite.find("cycler(")<rWhite.size()){
		inside=getInside(rWhite);
		if(inside.size()==0) inside="default";
		curCycler=inside;
	}else if(rWhite.find("=")<rWhite.size() &&
			 rWhite.find("==")!=rWhite.find("=") &&
			 rWhite.find("!=")!=rWhite.find("=")+1 ){
		Vector<std::string> ctL=parse_param(rWhite, '=');
		if(ctL.size()>=2){
			std::string lhs=ctL[0], rhs=collapsVS(ctL, 2, ctL.size());;
			myParse.parse(rhs);
			myParse.addVar(lhs, myParse());
		}

	}else if(rWhite.find("pulse(")<rWhite.size() ||rWhite.find("delay(")<rWhite.size()){
		Vector<std::string> ctL=parse_param(rWhite, '|');
		tmPD.rotorangle=curRot;
		tmPD.amp.resize(ctL.size());
		tmPD.offset.resize(ctL.size());
		for(int kk=0;kk<ctL.size();++kk){
			tmPD.amp[kk]=curAmp;
			tmPD.offset[kk]=curOff;
		}
		tmPD.usesys=curSys;
		tmPD.usepow=curPow;
		tmPD.cycler=curCycler;
		added=(tmPD.read(rWhite));
	//the Cs and Rs
	}
	curCycler="Ie";
	if(added){
		double dt=tmPD.t2-tmPD.t1;
		tmPD.t1=currentTime;
		tmPD.t2=currentTime+dt;
		currentTime=tmPD.t2;

		pulses.push_back(tmPD);
		//cout<<pulses.size()<<endl;
		curAdd++;
		added=false;
	}
	//std::cout<<curLin<<std::endl<<tmPD;
}
//this is the master decider function
// and snaps it to the proper above function
bool SequenceParse::decide(std::string inSqe)
{
	if(inSqe.find("pulse(")<inSqe.size() ||
		inSqe.find("delay(")<inSqe.size() ||
		inSqe.find("amplitude(")<inSqe.size() ||
		inSqe.find("offset(")<inSqe.size() ||
		inSqe.find("rotor(")<inSqe.size() ||
		inSqe.find("spinsys(")<inSqe.size() ||
		inSqe.find("powder(")<inSqe.size() ||
		inSqe.find("cycler(")<inSqe.size() )
	{
			addPulses(inSqe);
	}else if(addGvars &&
	      inSqe.find("=")<inSqe.size() &&
		   inSqe.find("==")!=inSqe.find("=") &&
		   inSqe.find("!=")!=inSqe.find("=")+1)
	{
		addGlobalVar(inSqe);
	}else{
		ScriptParse::decide(inSqe);
	}
	return true;
}

bool SequenceParse::decideSeP(std::string inSqe)
{
	if(inSqe.find("pulse(")<inSqe.size() ||
		inSqe.find("delay(")<inSqe.size() ||
		inSqe.find("amplitude(")<inSqe.size() ||
		inSqe.find("offset(")<inSqe.size() ||
		inSqe.find("rotor(")<inSqe.size() ||
		inSqe.find("rotor=")<inSqe.size() ||
		inSqe.find("spinsys(")<inSqe.size() ||
		inSqe.find("powder(")<inSqe.size() ||
		inSqe.find("cycler(")<inSqe.size()){
			addLine(inSqe);
/*	if(inSqe.find("pulse(")<inSqe.size()){
		addPulse(inSqe);
	}else if(inSqe.find("delay(")<inSqe.size()){
		addDelay(inSqe);
	}else if(inSqe.find("on(")<inSqe.size()){
		addOn(inSqe);
	}else if(inSqe.find("amplitude(")<inSqe.size()){
		addAmplitude(inSqe);
	}else if(inSqe.find("offset(")<inSqe.size()){
		addOffset(inSqe);
	}else if(inSqe.find("rotor(")<inSqe.size()){
		addRotor(inSqe);
	}else if(inSqe.find("powder(")<inSqe.size()){
		addPowder(inSqe);
	}else if(inSqe.find("spinsys(")<inSqe.size()){
		addSpinSys(inSqe);
	}else if(inSqe.find("cycler(")<inSqe.size()){
		addCycler(inSqe);*/
	}else{
		ScriptParse::decide(inSqe);
	}

	return true;
}



#endif
