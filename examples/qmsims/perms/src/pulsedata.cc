#include "blochlib.h"

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

#ifndef __Pulse_Data_cc__
#define __Pulse_Data_cc__ 1

#include "pulsedata.h"

using namespace BlochLib;

//the empty constructor
PulseData::PulseData():
	rotorangle(acos(1./sqrt(3.0))),
	t1(0.0), t2(0.0),
	cycler("Fe"), myParse(NULL)
{
	calcTime();
}

//a basic pulse constructor
PulseData::PulseData(double ina, double inp, double ino,double inamp, std::string inon):
	rotorangle(acos(1./sqrt(3.0))),
	angle(1,ina), phase(1,inp), offset(1,ino),amp(1,inamp),
	t1(0), t2(ina/360.0*1.0/inamp),
	delay(1,t2),dtDelay(0),incDelay(false),
	on(1,inon),cycler("Fe"), myParse(NULL)
{}

//create a Pulse data from a delay
PulseData::PulseData(double del):
	rotorangle(acos(1./sqrt(3.0))),
	angle(1,0), phase(1,0), offset(1,0), amp(1,0),
	t1(0), t2(del),
	delay(1,del),on(1,"1H"),cycler("Fe"), myParse(NULL)
{}


//the copy contructor
PulseData::PulseData(const PulseData &rhs)
{
	*this=rhs;
}

//the copy contructor
PulseData::~PulseData()
{
	myParse=NULL;
}


//the assignment operator
void PulseData::operator=(const PulseData &rhs)
{
	rotorangle=rhs.rotorangle;
	angle=rhs.angle;
	phase=rhs.phase;
	offset=rhs.offset;
	amp=rhs.amp;
	t1=rhs.t1;
	t2=rhs.t2;
	delay=rhs.delay;
	dtDelay=rhs.dtDelay;
	incDelay=rhs.incDelay;
	on=rhs.on;
	usesys=rhs.usesys;
	usepow=rhs.usepow;
	cycler=rhs.cycler;
	myParse=rhs.myParse;
}

//calculates the time for a pulse..
// Used internally to create a time list based on
// the input of an angle....
void PulseData::calcTime(double start)
{
	t1=start;
	for(int i=0;i<amp.size();++i){
		if(delay[i]==0.0 && angle[i]!=0 && amp[i]!=0){
			t2=t1+angle[i]/360.0/amp[i];
		}else if(delay[i]!=0){
			t2=t1+delay[i];
		}
	}
}



// THe pulse matrix
hmatrix PulseData::Pulse(SolidSys &A)
{
	matrix tm=A.F0();
	for(int j=0;j<amp.size();++j){
		if(amp[j]!=0){
			for(int i=0;i<A.size();++i){
				if(A.symbol(i)==on[j]){
					tm+= (amp[j]*(A.Ix(i)*cos(phase[j])
								+ A.Iy(i)*sin(phase[j]))
								+offset[j]*A.Iz(i));
				}
			}
		}
	}
	return tm;
}

//Simple text based output
std::ostream &operator<<(std::ostream &oo,const PulseData &out)
{
	out.print(oo, false);
	return oo;
}

//the pusle error function
void PulseData::PulseErr(std::string mySec, const char * file, int line)
{
	std::cerr<<"error:  PulseData::read"<<std::endl;
	std::cerr<<" Bad 'pulse' usage for \""<<mySec<<"\""<<std::endl;
	std::cerr<<" should be \"spin:pulse(time)\" \"spin:pulse(time, phase)\", "<<std::endl;
	std::cerr<<" \"spin:pulse(time, phase, amp)\" or \"spin:pulse(pulse(time, phase, amp,offset)\""<<std::endl;
	throw BL_exception(file, line);
}

//the dealy error function
void PulseData::DelayErr(std::string mySec,const char * file, int line)
{
	std::cerr<<"error:  PulseData::read"<<std::endl;
	std::cerr<<" Bad 'delay' usage for \""<<mySec<<"\""<<std::endl;
	std::cerr<<" should be \"spin:delay(time)\" " <<std::endl;
	throw BL_exception(file, line);
}

/**** REads a pulse line ***/
//a pulse parameter in its most generic form
// looks like this
//   spin:pulse(time, phase, amp, offset) | spin2:pulse(time, phase, amp, offset) ...
// where spin=="1H" or "13C" or similar
// the '|' means that both things are beign pulsed at the same time
// it can also look like this
//   spin:pulse(time, phase, amp, offset) | spin2:delay(time)
// where instead of a pulse on the second spin there is 'nothing'
// but a delay


//any number parameter (time, amp, phase, offset) can be
// a numeric expression parable by 'Parser'
// any global variables should be set in 'ParserGlobalVars'
// for this particular app, the 'SeqeunceRun' class handles this

bool PulseData::read(std::string inSec)
{
	Vector<std::string> multiP; //this holds the split '|' list
	Vector<std::string> spinP; //holds the split from ':'
	Vector<std::string> commP; //hold the spilt from ','
	std::string inside; //gets the inside of '(' ')'

	bool got=false;

	//first spin the '|' list
	multiP=parse_param(inSec, '|');
	amp.resize(multiP.size());
	delay.resize(multiP.size());
	offset.resize(multiP.size());
	angle.resize(multiP.size());
	phase.resize(multiP.size());
	incDelay.resize(multiP.size());
	dtDelay.resize(multiP.size());
	on.resize(multiP.size());

	for(int j=0;j<multiP.size();++j){
		spinP=parse_param(multiP[j], ':'); //get the spin
		if(spinP.size()<=1 || spinP.size()>2) PulseErr(inSec,__FILE__, __LINE__);
		//set the current On
		on[j]=spinP[0];

		//now see what the action is 'pulse' or 'delay'
		if(spinP[1].find("pulse(")<spinP[1].size()){
			inside=getInside(spinP[1]); //get the inside
			if(inside.size()==0) PulseErr(inSec,__FILE__, __LINE__);

			commP=parseComma(inside); //split the comma list
			if(commP.size()==0) PulseErr(inSec,__FILE__, __LINE__);

			double dt=0;
			if(commP.size()>=1){ //the time bit
				if(myParse!=NULL){
					myParse->parse(commP[0]);
					dt=(*myParse)();
				}else{
					dt=std::atof(commP[0].c_str());
				}
			}

			if(commP.size()>=2){ //the phase bit
				if(myParse!=NULL){
					myParse->parse(commP[1]);
					phase[j]=(*myParse)()*DEG2RAD;
				}else{
					phase[j]=std::atof(commP[1].c_str())*DEG2RAD;
				}
			}

			if(commP.size()>=3){ //the amplitude bit
				if(myParse!=NULL){
					myParse->parse(commP[2]);
					amp[j]=(*myParse)();
				}else{
					amp[j]=std::atof(commP[2].c_str());
				}
			}

			if(commP.size()>=4){ //the offset bit
				if(myParse!=NULL){
					myParse->parse(commP[3]);
					offset[j]=(*myParse)();
				}else{
					offset[j]=std::atof(commP[3].c_str());
				}
			}
			if(commP.size()>=1){
				t1=0;
				t2=dt+t1;
				delay[j]=0;
				incDelay[j]=false;
				dtDelay[j]=0;
			//	std::cout<<(dt*amp[j])*360.0<<" "<<dt<<" "<<1.0/amp[j]<<std::endl;
				angle[j]=(dt*amp[j])*360.0;
				got=true;
			}

		}else if(spinP[1].find("delay(")<spinP[1].size()){
			inside=getInside(spinP[1]); //get the inside
			if(inside.size()==0) DelayErr(inSec,__FILE__, __LINE__);

			commP=parseComma(inside); //split the comma list
			if(commP.size()==0) DelayErr(inSec,__FILE__, __LINE__);

			double dt=0;
			if(commP.size()>=1){ //the time bit
				if(myParse!=NULL){
					myParse->parse(commP[0]);
					dt=(*myParse)();

				}else{
					dt=std::atof(commP[0].c_str());
				}
				delay[j]=dt;
				t1=0;
				t2=dt;
				incDelay[j]=false;
				dtDelay[j]=dt;
				amp[j]=0;
				phase[j]=0;
				angle[j]=0;
				got=true;
			}
		}
	}
	return got;
}






/**** THe Main Pulse Data Vector Reader....****/
Vector<PulseData> readPulses(Parameters &pset, Parser *myP)
{	return readPulses(pset.section("pulses"), myP);	}

Vector<PulseData> readPulses(const Vector<std::string> &mySecin, Parser *myP)
{
	Vector<std::string>  curLin;
	double curRot=acos(1.0/sqrt(3.0)), curAmp=0, curOff=0;
	std::string curSys="default", curPow="default", curOn="1H", curCycler="Ie",
	            rWhite="", inside;
	double currentTime=0;
	Vector<PulseData> outdata;
	int i, curAdd=0;
	bool added=false;


	Vector<std::string> mySec=paramStrip(mySecin); //clean up the input
	for(i=0;i<mySec.size();++i){
		PulseData tmPD;
		tmPD.myParse=myP;
		rWhite=removeWhite(mySec[i]);
		//these are the global parmeter functions
		if(rWhite.find("rotor(")<rWhite.size()){
			inside=getInside(rWhite);
			if(myP!=NULL){
				myP->parse(inside);
				curRot=(*myP)()*DEG2RAD;
			}else{
				curRot=std::atof(inside.c_str())*DEG2RAD;
			}
		}else if(rWhite.find("spinsys(")<rWhite.size()){
			inside=getInside(rWhite);
			if(inside.size()==0) inside="default";
			curSys=inside;
		}else if(rWhite.find("powder(")<rWhite.size()){
			inside=getInside(rWhite);
			if(inside.size()==0) inside="default";
			curPow=inside;
		}else if(rWhite.find("amplitude(")<rWhite.size()){
			inside=getInside(rWhite);
			if(myP!=NULL){
				myP->parse(inside);
				curAmp=(*myP)();
			}else{
				curAmp=std::atof(inside.c_str());
			}
		}else if(rWhite.find("offset(")<rWhite.size()){
			inside=getInside(rWhite);
			if(myP!=NULL){
				myP->parse(inside);
				curOff=(*myP)();
			}else{
				curOff=std::atof(inside.c_str());
			}
		}else if(rWhite.find("cycler(")<rWhite.size()){
			inside=getInside(rWhite);
			if(inside.size()==0) inside="default";
			curCycler=inside;
		}else if(rWhite.find("=")<rWhite.size() &&
		         rWhite.find("==")!=rWhite.find("=") &&
		         rWhite.find("!=")!=rWhite.find("=")+1 ){
			if(myP!=NULL){
				Vector<std::string> ctL=parse_param(rWhite, '=');
				if(ctL.size()>=2){
					std::string lhs=ctL[0], rhs=collapsVS(ctL, 2, ctL.size());;
					myP->parse(rhs);
					myP->addVar(lhs, (*myP)());
				}
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
			outdata.push_back(tmPD);
			curAdd++;
			added=false;
		}
		//std::cout<<curLin<<std::endl<<tmPD;
	}
	return outdata;
}

void PulseData::printHead(std::ostream &oo)
{
	std::string top;
	static std::string rottr="rotor(deg) ";
	static std::string systr="spinsys   ";
	static std::string powtr="powder   ";
	static std::string statr="start T(s)   ";
	static std::string endtr="end T(s)     ";
	static std::string instr="input...  ";
	top=rottr+systr+powtr+statr+endtr+instr;
	oo<<top<<std::endl;
	oo<<std::string(top.size(), '-')<<std::endl;
}

void PulseData::print(std::ostream &oo, bool printhead) const
{
	std::string tmp;
	std::string top;
	static std::string rottr="rotor(deg) ";
	static std::string systr="spinsys   ";
	static std::string powtr="powder   ";
	static std::string statr="start T(s)   ";
	static std::string endtr="end T(s)     ";
	static std::string instr="input...  ";
	top=rottr+systr+powtr+statr+endtr+instr;
	if(printhead) printHead(oo);

	tmp=dbtost(rotorangle*RAD2DEG, "%3.2f");
	oo<<tmp<<std::string(rottr.size()-tmp.size(), ' ');

	tmp=usesys;
	oo<<tmp<<std::string(systr.size()-tmp.size(), ' ');

	tmp=usepow;
	oo<<tmp<<std::string(powtr.size()-tmp.size(), ' ');

	tmp=dbtost(t1, "%1.4g");
	oo<<tmp<<std::string(statr.size()-tmp.size(), ' ');

	tmp=dbtost(t2, "%1.4g");
	oo<<tmp<<std::string(endtr.size()-tmp.size(), ' ');

	for(int i=0;i<amp.size();++i){
		if(i!=0 )	oo<<" | ";
		if(amp[i]!=0){
			oo<<on[i]<<":"<<"pulse(angle="<<(t2-t1)*(amp[i])*360.0<<","
			              <<"phase="<<phase[i]*RAD2DEG<<","
			              <<"amp="<<amp[i]<<","
			              <<"offset="<<offset[i]<<") ";
		}else{
			oo<<on[i]<<":"<<"delay(dt="<<(t2-t1)<<") ";
		}
	}
	oo<<std::endl;

}

std::ostream &operator<<(std::ostream &oo, Vector<PulseData> &pdata)
{
	for(int i=0;i<pdata.size();++i){
		if(i==0){
			pdata[i].print(oo, true);
		}else{
			pdata[i].print(oo);
		}
	}
	return oo;
}
#endif

