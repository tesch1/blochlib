

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



/**************** Sequence Parser bits *************/
SequenceParse::SequenceParse(Parameters &pset)
{
	data=pset.section("");
	ScriptParse::parse(data);
}

SequenceParse::SequenceParse(const Vector<std::string> &pset)
{
	data=pset;
	ScriptParse::parse(data);
}

SequenceParse::SequenceParse(const SequenceParse &rhs):
	ScriptParse(rhs)
{

	if(this==&rhs) return;
	*this=rhs;
}

void SequenceParse::operator=(const SequenceParse &rhs)
{
	if(this==&rhs) return;
	ScriptParse::operator=(rhs);
	sequence=rhs.sequence;
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

//this should take in a string like "pulse(#,#)" and add the
//proper line to the sequence vector
void SequenceParse::addPulse(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	//std::cout<<"screw"<<inSqe<<std::endl;
	if(tmS.find("pulse(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'pulse' usage for \"")+inSqe+"\""+
			"\n should be \"pulse(time)\" \"pulse(time, phase)\", "+
			"\n \"pulse(time, phase, amp)\" or \"pulse(time, phase, amp,offset)\"")			
		}

		out+="pulse ";
		if(ps.size()>=1){
			myParse.parse(ps[0]);
			out+=" "+dbtost(myParse(),"%.9f");
		}

		if(ps.size()>=2){
			myParse.parse(ps[1]);
			out+=" "+dbtost(myParse(),"%.9f");
		}

		if(ps.size()>=3){
			myParse.parse(ps[2]);
			out+=" amplitude "+dbtost(myParse(),"%.9f");
		}

		if(ps.size()>=4){
			myParse.parse(ps[3]);
			out+=" offset "+dbtost(myParse(),"%.9f");
		}
	}
	sequence.push_back(out);
}



//this should take in a string like "delay(#,#)" and add the
//proper line to the sequence vector
void SequenceParse::addDelay(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("delay(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'delay' usage for \"")+inSqe+"\""+
			"\n should be \"delay(time)\" \"delay(time, dt)\", ")
		}

		out+="delay";
		if(ps.size()>=2){
			myParse.parse(ps[1]);
			out+="step ";
		}

		if(ps.size()>=1){
			myParse.parse(ps[0]);
			out+=" "+dbtost(myParse(),"%.9f");
		}

		if(ps.size()>=2){
			myParse.parse(ps[1]);
			out+=" dt "+dbtost(myParse(),"%.9f");
		}
	}
	sequence.push_back(out);
}


//this should take in a string like "rotor(#)" and add the
//proper line to the sequence vector
void SequenceParse::addRotor(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("rotor(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'rotor' usage for \"")+inSqe+"\""+
			"\n should be \"rotor(angle)\" (angle in DEGREES), ")
		}

		out+="rotor";
		if(ps.size()>=1){
			myParse.parse(ps[0]);
			out+=" "+dbtost(myParse(),"%.9f");
		}
	}
	sequence.push_back(out);
}

//this should take in a string like "cycler(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addCycler(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("cycler(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'cycler' usage for \"")+inSqe+"\""+
			"\n should be \"cycler(SpinOp)\" (a HamiltonianGen Type)  ")
		}

		out+="cycler";
		if(ps.size()>=1){
			out+=" "+ps[0];
		}
	}
	sequence.push_back(out);
}

//this should take in a string like "amplitude(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addAmplitude(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("amplitude(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'amplitude' usage for \"")+inSqe+"\""+
			"\n should be \"amplitude(amp)\" (in Hz) ");
		}

		out+="amplitude";
		if(ps.size()>=1){
			myParse.parse(ps[0]);
			out+=" "+dbtost(myParse(),"%.9f");
		}
	}
	sequence.push_back(out);
}

//this should take in a string like "offset(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addOffset(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("offset(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'offset' usage for \"")+inSqe+"\""+
			"\n should be \"offset(offs)\" (in Hz) ")
		}

		out+="offset";
		if(ps.size()>=1){
			myParse.parse(ps[0]);
			out+=" "+dbtost( myParse(),"%.9f");
		}
	}
	sequence.push_back(out);
}

//this should take in a string like "spinsys(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addSpinSys(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("spinsys(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'spinsys' usage for \"")+inSqe+"\""+
			"\n should be \"spinsys(SysName)\" (in the 'spins' section) ")
		}

		out+="useSys";
		if(ps.size()>=1){
			out+=" "+ps[0];
		}
	}
	sequence.push_back(out);
}

//this should take in a string like "powder(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addPowder(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("powder(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'powder' usage for \"")+inSqe+"\""+
			"\n should be \"powder(PowName)\" (in the 'powders' section) ")
		}

		out+="usePowder";
		if(ps.size()>=1){
			out+=" "+ps[0];
		}
	}
	sequence.push_back(out);
}
//this should take in a string like "on(spinop)" and add the
//proper line to the sequence vector
void SequenceParse::addOn(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("on(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'on' usage for \"")+inSqe+"\""+
			"\n should be \"on(SpinName)\" (like '1H') ")
		}

		out+="on";
		if(ps.size()>=1){
			out+=" "+ps[0];
		}
	}
	sequence.push_back(out);
}

//this is the master decider function
// and snaps it to the proper above function
bool SequenceParse::decide(std::string inSqe)
{
	if(inSqe.find("pulse(")<inSqe.size() ||
		inSqe.find("delay(")<inSqe.size() ||
		inSqe.find("amplitude(")<inSqe.size() ||
		inSqe.find("offset(")<inSqe.size() ||
		inSqe.find("rotor=")<inSqe.size() ||
		inSqe.find("wr=")<inSqe.size() ||
		inSqe.find("rotor(")<inSqe.size() ||
		inSqe.find("spinsys(")<inSqe.size() ||
		inSqe.find("powder(")<inSqe.size() ||
	inSqe.find("U(")<inSqe.size() ||
	inSqe.find("cycler(")<inSqe.size() ||
		(inSqe.find("=")<inSqe.size() && inSqe.find("==")!=inSqe.find("=") && inSqe.find("!=")!=inSqe.find("=")+1 )
		){
			addLine(inSqe);
	/*
	if(inSqe.find("pulse(")<inSqe.size()){
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
	}else if(inSqe.find("spinsys(")<inSqe.size()){
		addSpinSys(inSqe);
	}else if(inSqe.find("powder(")<inSqe.size()){
		addPowder(inSqe);
	}else if(inSqe.find("cycler(")<inSqe.size()){
		addCycler(inSqe);
	*/
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
	inSqe.find("U(")<inSqe.size() ||
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
