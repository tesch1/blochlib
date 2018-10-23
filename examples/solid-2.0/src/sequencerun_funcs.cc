

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08.5.02
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
 	sequencerun_funcs.cc -->
 	this maintains the pulse sections and performs the
 	propogation and fid collection (the main runner class)

 	it is public to sequence parse so it can also parse
 	'non-section sequences'

 	it adds  these function for scripting
 	(it also has the sequenceparse ones well)

 	use(seqname)
 	  --> use this subsection
 	use(seqname, repeat)
 	  --> use the subsection reapeating it n times
 	use(seqname, repeat, hold)
 	  --> use the subsection reapeating it n times
 	  --> and 'held' so that it is only appiled once in 2D cycles

 	fid()
 	  --> collect an fid at this poin in the prop run
 	  --> if the sequence is ptop it will collect one point

 	fid(index)
 	  --> collect an fid at this poin in the prop run
 	  --> if the sequence is ptop it will collect one point
 	  --> and ADDS it to an fid at index

 	savefid()
 	  --> saves an fid into a file...using the default out name

 	savefid(name)
 	  --> saves an fid into a file of name 'name'

 	ptop()
 	  --> spcifies that the sequence is point to point

	rotor(angle)
	  --> set the collection rotor angle

	spinsys(sys)
	  --> set the collection spin system

	powder(angle)
	  --> set the collection power average

NOTE:: any variables decalred outside the 'subsections'
in the pulses parameter set will be global to all the
subsections....

it also hold the main ParamSet

*/

#ifndef __sequencerun_funcs_cc__
#define __sequencerun_funcs_cc__ 1

#include "sequencerun.h"
#include <signal.h>

extern "C"{
void SC__sigdie(int);
void SC__sigdie(int)
{
	signal(SIGHUP, SC__sigdie); /* set function calls */
	signal(SIGQUIT ,SC__sigdie);
	signal(SIGKILL  , SC__sigdie);
	signal(SIGINT  , SC__sigdie);
	signal(SIGSEGV  , SC__sigdie);
	signal(SIGABRT , SC__sigdie);
}
}


/******************* THe MAIN FUNCTIONS... *****/
void SequenceRun::doPulseData(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	Vector<PulseData> pdata(1);
	//double dt=0;

	if(tmS.find("pulse(")<tmS.size() ||tmS.find("delay(")<tmS.size()  ){

		pdata[0].rotorangle=Params.getD("rotor")*DEG2RAD;
		Vector<std::string> ctL=parse_param(tmS, '|');
		pdata[0].amp.resize(ctL.size(),amp);
		pdata[0].offset.resize(ctL.size(),offset);
		pdata[0].usesys=usesys;
		pdata[0].usepow=usepow;
		pdata[0].cycler=cycler;
		pdata[0].myParse=&myParse;
		bool added=(pdata[0].read(tmS));
		if(added){
			pdata[0].t1+=currentTime;
			pdata[0].t2+=currentTime;
			if(!toRun){
				pdata[0].print(std::cout);
			}else{
				myProp.pdata=pdata;
				PropList[numProp].docycler=false;
				PropList[numProp].hold=false;
				PropList[numProp].repeat=1;
				PropList[numProp].U=myProp.propogator(powpt);
				numProp++;
			}
			currentTime=pdata[0].t2;
		}
	}
}

void SequenceRun::setAmplitude(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("amplitude(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'amplitude' usage for \"")+inSqe+"\""+
			"\n should be \"amplitude(Hz)\" (number in Hz) ")
		}

		if(ps.size()>=1){
			myParse.parse(ps[0]);
			amp=myParse();
		}
	}
}

void SequenceRun::setOffset(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("offset(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'offset' usage for \"")+inSqe+"\""+
			"\n should be \"offset(Hz)\" (num in Hz) ")
		}

		if(ps.size()>=1){
			myParse.parse(ps[0]);
			offset=myParse();
		}
	}
}

void SequenceRun::setSystem(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("spinsys(")<tmS.size())
	{
		tmS=getInside(tmS);
		if(tmS.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'spinsys' usage for \"")+inSqe+"\""+
			"\n should be \"spinsys(sys)\" ")
		}

		usesys=tmS;
		myProp.setSystem(usesys);
	}
}

void SequenceRun::setOn(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("on(")<tmS.size())
	{
		tmS=getInside(tmS);
		if(tmS.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'on' usage for \"")+inSqe+"\""+
			"\n should be \"on(SpinLabel)\" ")
		}

		on=tmS;
	}
}

void SequenceRun::setPowder(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("powder(")<tmS.size())
	{
		tmS=getInside(tmS);
		if(tmS.size()==0 )
		{
			BLEXCEPTION(std::string(" Bad 'powder' usage for \"")+inSqe+"\""+
			"\n should be \"powder(pow)\" ")
		}
		if(usepow!=tmS){
			usepow=tmS;
			myProp.setPowder(usepow);
			fid1D.fill(0);
			powpt=0;
		}
	}
}

void SequenceRun::doCycler(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("cycler(")<tmS.size())
	{
		tmS=getInside(tmS);
		if(tmS.size()==0)
		{
			BLEXCEPTION(std::string(" Bad 'cycler' usage for \"")+inSqe+"\""+
			"\n should be \"cycler(SpinOp)\" ")
		}

		cycler=tmS;

		if(!toRun){
			std::cout<<"Ro will be 'cycled' with \""<<cycler<<"\""<<std::endl;
		}else{
			PropList[numProp].props=&myProp;

			PropList[numProp].cycler=cycler;
			PropList[numProp].docycler=true;
			numProp++;

		//scomplex norm=trace(trac);
		//	ro=trace(ro,adjoint(trac))*trac;
		}
		cycler="Ie";
	}
}



//this should take in a string like "ptop()" and add the
//proper line to the sequence vector
void SequenceRun::addPtop(std::string inSqe)
{
	std::string out;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("ptop(")<tmS.size())
	{
		tmS=getInside(tmS);
		isPtoP=true;
	}
	if(!toRun && isPtoP){
		std::cout<<"Point-TO-Point Spectrum"<<std::endl;
	}
}

//this should take in a string like "A=B" and add the
//proper line to the sequence vector
void SequenceRun::addGVar(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("=")<tmS.size()
	   && inSqe.find("==")!=inSqe.find("=")
	   && inSqe.find("!=")!=inSqe.find("=")-1 )
	{
		ps=parse_param(tmS, '=');
		if(ps.size()<=1)
		{
			BLEXCEPTION(std::string(" Bad 'variable' usage for \"")+inSqe+"\""+
			"\n should be \"A=(b*...)\" ")
		}
		std::string lhs=ps[0];
		std::string rhs=collapsVS(ps,2, ps.size());

		if(ps.size()>=1){
			Params.set(lhs, rhs);
		}
	}
}

//Shows the current pulse sequence
void SequenceRun::Shows(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("show()")<tmS.size())
		toRun=false;
}

//sets ro to the roeq if nothing is in the argument,
//is something is pressent is sets ro to that value...
void SequenceRun::setRo(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("ro(")<tmS.size()){
		tmS=getInside(tmS);
		if(tmS.size()!=0){
			if(tmS=="roeq") tmS=Params.getS("roeq");
			HamiltonianGen myGen;
			ro=myGen.Hamiltonian((*myProp.curSys), tmS,
							 myProp.curPow->theta(powpt),
							 myProp.curPow->phi(powpt),
							 myProp.curPow->gamma(powpt));
		}else{
			ro=roeq;
		}
	}
}

//sets a new detection input string...
// the detection is set in the 'doFID' function
// so all we need to do here is reset the input string
// this could be done with a 'detect=SPinOp" command as well
// but ust to make the sysntax similar to 'ro(SPinOp)'
void SequenceRun::setDetect(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("detect(")<tmS.size()){
		tmS=getInside(tmS);
		if(tmS.size()!=0){
			Params.set("detect", tmS);
		}
	}
}

//resets all the main parameters to their 0 values
void SequenceRun::ResetPars(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("reset()")<tmS.size()){
		fid2D.fill(0);
		fid1D.fill(0);
		currentTime=0;
		numProp=0;
		ro=roeq;
		pt2D=0;
		canSave=false;
		notSaved=true;
	}
}

//sets the 2D flag
void SequenceRun::to2D(std::string inSqe)
{
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("2D(")<tmS.size())
	{
		is2D=true;
		if(pt2D==-1) pt2D=0;
		tmS=getInside(tmS);
		if(tmS.size()>0){
			myParse.parse(tmS);
			Params.set("npts2D", myParse());
			fid2D.resizeAndPreserve(Params.getI("npts2D"), Params.getI("npts1D"), 0.0);
		}
	}
}

//this functions alters (adds, or changes) an interaction in the
// current spin system
// uses the syntax
//  alterSys(what, number)
// where 'what' is the input syntax for a 'setSpinParam' function
// in the SolidSys class
//  'C1iso' means the second spins isotropic shift
//  'C1del' means the second spins anisotropic shift
//  'C1eta' means the second spins assymetry shift
//  'D01' means the dipole coupling between spin 0 and 1
//  'J01' means the jcoup between spin 0 and 1
//  'Q1' means the qaud coupling on spin 1
//  'Q1eta' means the assymetry of the qud on spin 1
// it will look for the interaction
void SequenceRun::doAlterSys(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("alterSys(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()<2)
		{
			BLEXCEPTION(std::string(" Bad 'use' usage for \"")+inSqe+"\""+
			"\n should be \"alterSys(what, number)\" "+
			"\n where what is something like 'C1eta', 'D01', 'C1iso', etc ")
		}
		myParse.parse(ps[1]);
		if(toRun){
			myProp.curSys->setSpinParam(ps[0], myParse());
		}else{
			std::cout<<" Spin System Alteration: "<<ps[0]<<" to "<<myParse()<<std::endl;
		}
	}
}

//dumps out the current propogators
void SequenceRun::showProps(std::string inSqe)
{
	Vector<std::string> ps;
std::string tmS=removeWhite(inSqe);
	if(tmS.find("U(")<tmS.size())
	{
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()<1){
			for(int i=0;i<numProp;i++){
			std::cout<<"---Propagator "<<i+1<<" : "<<std::endl;
			std::cout<<PropList[i]<<std::endl;
			}
		}else{
			for(int i=0;i<ps.size();i++){
				myParse.parse(ps[i]);
				int indx=int(myParse());
				if(numProp < indx){
std::cout<<"---Propagator "<<indx<<" : "<<std::endl;
std::cout<<PropList[indx]<<std::endl;
				}else{
					std::cout<<" Propogator index '"<<indx<<"'does not exsist"<<std::endl;
				}
			}
		}
	}
}

//if somethign happens, like a [ctlr]-C or an error
// setting this flag will dump the current state
// of the system for looking at
void SequenceRun::setDumpOnDie(std::string inSqe)
{
	Vector<std::string> ps;
std::string tmS=removeWhite(inSqe);
	if(tmS.find("dumpOnDie()")<tmS.size())
	{
		dumpOnDie = true;
	}
}


//insteead of performing an acctuall seuqnce simply
//dumps the 'PulseData vector generated to std::cout
// parses 'use(seqname, repeta, hold)'
void SequenceRun::showUsed(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	bool held=false;
	int repeat=1;
	bool rehash=false;
	if(tmS.find("use(")<tmS.size())
	{

		rehash=true;
		if(tmS.find("reuse(")<tmS.size()) rehash=false;
		tmS=getInside(tmS);
		ps=parseComma(tmS);
		if(tmS.size()==0 || ps.size()<1)
		{
			BLEXCEPTION(std::string(" Bad 'use' usage for \"")+inSqe+"\""+
			"\n should be \"use(seq)\", \"use(seq, repeat)\", \"use(seq, repeat, hold)\" ")
		}

		if(ps.size()>=2){
			myParse.parse(ps[1]);
			repeat=int(myParse());
			if(repeat<0){
				BLEXCEPTION(std::string(" Bad 'repeat' usage for \"")+inSqe+"\""+
				"\n repeat should be >=0 (use none or more...) ")
			}
		}

		if(ps.size()>=3){
			if(ps[2]=="hold") held=true;
			else{
				myParse.parse(ps[2]);
				held=bool(myParse());
			}
		}



		SubPulseMapIter jj=subProc.find(ps[0]);

		if(jj!=subProc.end()){
			jj->second.hold=held;
			jj->second.repeat=repeat;
			if(!toRun){
					//set the Bfield all the systems
					myProp.Systems->setBfield(Params.getD("Bfield"));
					std::cout<<jj->first<<"---------------------------"<<std::endl;
					if(held) std::cout<<"Sequence will be held (valid for 2Ds only)"<<std::endl;
					std::cout<<"Sequence will repeated \""<<repeat<<"\" times"<<std::endl;
					if(jj->second.haveParsed && rehash) jj->second.haveParsed=false;
					Vector<PulseData> tmm=jj->second.parse();
					if(!jj->second.hold)
						if(jj->second.parse().size()>0)
							currentTime=jj->second.parse()[jj->second.parse().size()-1].t2*double(repeat);

					for(int i=0;i<tmm.size();++i){
						if(i==0) tmm[i].print(std::cout, true);
						else tmm[i].print(std::cout);
					}

			}else{
				//set the Bfield all the systems
				myProp.Systems->setBfield(Params.getD("Bfield"));
				if(jj->second.haveParsed && rehash) jj->second.haveParsed=false;
				try{
					if(jj->second.repeat>0) jj->second.propogator(powpt, pt2D);
				}catch(BL_exception e){
					e.print(std::cerr);
					BLEXCEPTION("leaving evaluation");
				}
				if(!jj->second.hold)
					if(jj->second.parse().size()>0)
						currentTime=jj->second.parse()[jj->second.parse().size()-1].t2*double(repeat);

				if(jj->second.repeat>0){
					PropList[numProp]=(jj->second);
					numProp++;
				}
			}
			
		}else{
			BLEXCEPTION(std::string(" Sub Pulse section does NOT exsist use... \"")+inSqe+"\"")
		}
	}
}

void SequenceRun::doFID(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	int idx=-1;

	//set the Bfield all the systems
	myProp.Systems->setBfield(Params.getD("Bfield"));

	if(tmS.find("fid(")<tmS.size())
	{
		tmS=getInside(tmS);
		fid1D.resizeAndPreserve(Params.getI("npts1D"));

		if(tmS.size()>=1){
			myParse.parse(tmS);
			idx=int(myParse());
			if(!is2D){
				if(idx>=fid1D.size()){
					fid1D.resizeAndPreserve(idx);
					fid1D[idx]=0.0;
					Params.set("npts1D", idx+1);
				}
			}else{
				if(idx>=fid2D.rows()){
					fid2D.resizeAndPreserve(idx+1, Params.getI("npts1D"), 0.0);
					Params.set("npts2D", idx+1);
				}
			}
		}
	}else{
		return;
	}

	if(!toRun){
		std::cout<<" collect FID ";
		if(idx!=-1) std::cout<<" at point "<<idx;
		std::cout<<std::endl;
		return;
	}

	myProp.setPowder(usepow);
	myProp.setSystem(usesys);
	myProp.curSys->setRotorAngles(0.0,Params.getD("rotor")*DEG2RAD);

//the FID collectors
	oneFID<SolidSys, hmatrixs> hFID;
	oneFID<SolidSys, dmatrixs> dFID;

	if(myProp.curSys->isDiagonal()){
		dFID.setFunction(myProp.curSys);
		dFID.setSize(Params.getI("npts1D"));
		dFID.setSw(Params.getD("sw"));
		dFID.setWr(Params.getD("wr"));
		dFID.setRotorAngle(Params.getD("rotor"));
		if(myProp.pdata.size()>0) dFID.setStartTime(currentTime);
	}else{
		hFID.setFunction(myProp.curSys);
		hFID.setSize(Params.getI("npts1D"));
		hFID.setSw(Params.getD("sw"));
		hFID.setWr(Params.getD("wr"));
		hFID.setRotorAngle(Params.getD("rotor"));
		if(myProp.pdata.size()>0) hFID.setStartTime(currentTime);
	}

	static HamiltonianGen myGen;
	detect=myGen.Hamiltonian((*myProp.curSys), Params.getS("detect"),
							 myProp.curPow->theta(powpt),
							 myProp.curPow->phi(powpt),
							 myProp.curPow->gamma(powpt));

	if(isPtoP){
		if(!is2D){
			if(idx!=-1){
				for(int jj=0;jj<numProp;++jj){
					try{
						PropList[jj].propogate(ro, powpt, 0);
					}catch(BL_exception e){
						e.print(std::cerr);
						BLEXCEPTION("cannot evalutate point-to-point sectrum");
					}
				}
				fid1D[idx]+=myProp.curPow->weight(powpt)*trace(ro, detect)/trace(detect, adjoint(detect));

			//	if(MPIworld.master()) std::cout<<"    on Point: "<<idx<<"     ";
				if(idx==0 && currentTime!=0) Params.set("sw", 1.0/currentTime);
			}else{
				for(int i=0;i<fid1D.size();++i){
					fid1D[i]+=myProp.curPow->weight(powpt)*trace(ro, detect)/trace(detect, adjoint(detect));
					if(i!=fid1D.size()-1){
						for(int jj=0;jj<numProp;++jj){
							try{
								PropList[jj].propogate(ro, powpt, i);
							}catch(BL_exception e){
								e.print(std::cerr);
								BLEXCEPTION("cannot evalutate point-to-point sectrum");
							}
						}
					}
				}
				//if( currentTime!=0) Params.set("sw", 1.0/currentTime);
			}
			numProp=0;
		}else{
			if(idx!=-1){
				//if(MPIworld.master()) std::cout<<"    on FID: "<<idx<<"     ";
				for(int i=0;i<fid1D.size();++i){
					for(int jj=0;jj<numProp;++jj){
						try{
							PropList[jj].propogate(ro, powpt, i);
						}catch(BL_exception e){
							e.print(std::cerr);
							BLEXCEPTION("cannot evalutate point-to-point sectrum");
						}
					}
					fid1D[i]=myProp.curPow->weight(powpt)*trace(ro, detect);
				}
				fid1D/=trace(detect, adjoint(detect));
				//if( currentTime!=0) Params.set("sw", 1.0/currentTime);

				Range all(Range::Start, Range::End);
				fid2D.put(idx,all,fid2D(idx, all)+fid1D(all));
				numProp=0;
			}
		}
	}else{
		if(!is2D){
			if(idx!=-1){
				for(int jj=0;jj<numProp;++jj){
					PropList[jj].propogate(ro, powpt, 0);
				}
				fid1D[idx]+=myProp.curPow->weight(powpt)*trace(ro, detect)/trace(detect, adjoint(detect));
			}else{
				for(int jj=0;jj<numProp;++jj){
					PropList[jj].propogate(ro, powpt, 0);
					//std::cout<<PropList[jj].holdU;
				}
				if(myProp.curSys->isDiagonal()){
					fid1D+=myProp.curPow->weight(powpt)*
						   dFID.FID(ro, detect,
							 myProp.curPow->theta(powpt),
							myProp.curPow->phi(powpt),
							 myProp.curPow->gamma(powpt))/trace(detect, adjoint(detect));
					Params.set("sw", dFID.sw());
				}else{
					fid1D+=myProp.curPow->weight(powpt)*
						   hFID.FID(ro, detect,
							 myProp.curPow->theta(powpt),
							myProp.curPow->phi(powpt),
							 myProp.curPow->gamma(powpt))/trace(detect, adjoint(detect));
					Params.set("sw", hFID.sw());
				}
			}
		}else{
			for(int jj=0;jj<numProp;++jj){
				PropList[jj].propogate(ro, powpt, pt2D);
			}
			Vector<scomplex> tmFID(Params.getI("npts1D"), 0);
			if(myProp.curSys->isDiagonal()){
				tmFID=myProp.curPow->weight(powpt)*
							dFID.FID(ro, detect,
							myProp.curPow->theta(powpt),
							myProp.curPow->phi(powpt),
							myProp.curPow->gamma(powpt))/trace(detect, adjoint(detect));
				Params.set("sw", dFID.sw());
			}else{
				tmFID=myProp.curPow->weight(powpt)*
							hFID.FID(ro, detect,
							myProp.curPow->theta(powpt),
							myProp.curPow->phi(powpt),
							myProp.curPow->gamma(powpt))/trace(detect, adjoint(detect));
				Params.set("sw", hFID.sw());
			}
			Range all(Range::Start, Range::End);

			if(idx!=-1){
				//if(idx==0) Params.set("sw2", 1.0/currentTime);

				//if(MPIworld.master()) std::cout<<"    on FID: "<<idx<<"     ";
				fid2D.put(idx,all,fid2D(idx, all)+tmFID(all));

			}else{
				fid2D.put(pt2D,all, fid2D(pt2D,all)+tmFID(all));
				pt2D++; //advance a point if not doing an index value
			}
			numProp=0; //reset the list to zero
			//currentTime=0; //reset the time list back to 0
		}
	}
	if(MPIworld.master()) std::cout<<"\r";

}

//can take 0 to 2 args
// savefid() --> simply dumps the last fid out (ither 1D or 2D)
// savefid(name) -> dumps thefid with the name
// savefid(name, i) --> if a 2 exp, will print out the i'th fid
// we can only do these if 'canSave' has been set

void SequenceRun::saveFID(std::string inSqe)
{
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);

	//recollect the fids
	if(canSave && toRun){
		MPIworld.reduce(fid1D, Reduce::Add);
		MPIworld.reduce(fid2D, Reduce::Add);
	}
	if(MPIworld.master()){
		int saveWhat=1;
		if(tmS.find("savefidmatlab(")<tmS.size()){
			saveWhat=10;
		}else if(tmS.find("savefidbinary(")<tmS.size()){
			saveWhat=20;
		}else if(tmS.find("savefid(")<tmS.size() || tmS.find("savefidtext(")<tmS.size()){
			saveWhat=1;
		}else{
			return;
		}


		tmS=getInside(tmS);
		ps=parseComma(tmS);
		std::string fname=Params.getS("filesave");
		if(ps.size()>=1) fname=ps[0];

		Vector<scomplex> curFID=fid1D;
		if(ps.size()<=1){
			if(canSave && toRun){
				if(pt2D==0)	curFID=fid2D.row(0);
				if(is2D){
					if(saveWhat==10){
						matstream matout(fname, std::ios::out | std::ios::binary);
						matout.put("vdat", fid2D);
						matout.close();
					}else if(saveWhat==20){
						std::fstream out(fname.c_str(), std::ios::binary | std::ios::out);
						if(out.fail()){
							std::cerr<<"Error: SaveFid \""<<Params.getS("filesave")<<"\" file cannot be opened for writing"<<std::endl;
							return;
						}

						out<<"npts1 ="<<fid2D.rows()<<std::endl;
						out<<"npts2 ="<<fid2D.cols()<<std::endl;
						double tmR, tmI;
						for(int p=0;p<fid2D.rows();++p){
							for(int h=0;h<fid2D.cols();++h){
								tmR=Re(fid2D(p,h));
								tmI=Im(fid2D(p,h));
								out.write((char *)&tmR, sizeof(double));
								out.write((char *)&tmI, sizeof(double));
							}
						}
					}else{
						std::ofstream out(fname.c_str());
						if(out.fail()){
							std::cerr<<"Error: SaveFid \""<<fname<<"\" file cannot be opened for writing"<<std::endl;
							return;
						}
						out<<"npts1 ="<<fid2D.rows()<<std::endl;
						out<<"npts2 ="<<fid2D.cols()<<std::endl;
						for(int p=0;p<fid2D.rows();++p){
							for(int h=0;h<fid2D.cols();++h){
								out<<p<<" "<<h<<" "<<Re(fid2D(p,h))<<" "<<Im(fid2D(p,h))<<std::endl;
							}
						}
					}

				}else{
					double dt=1/Params.getD("sw")/2;
					if(saveWhat==10){
						matstream matout(fname, std::ios::out | std::ios::binary);
						matout.put("vdat", fid1D);
						matout.close();
					}else if(saveWhat==20){
						std::fstream out(fname.c_str(), std::ios::binary | std::ios::out);
						if(out.fail()){
							std::cerr<<"Error: SaveFid \""<<Params.getS("filesave")<<"\" file cannot be opened for writing"<<std::endl;
							return;
						}

						out<<"npts1 ="<<fid1D.size()<<std::endl;
						double tmR, tmI;
						for(int p=0;p<fid1D.size();++p){
							tmR=Re(fid1D(p));
							tmI=Im(fid1D(p));
							out.write((char *)&tmR, sizeof(double));
							out.write((char *)&tmI, sizeof(double));
						}
						out.close();
					}else{
						plotterFID(curFID, fname, dt);
					}
				}
				notSaved=false;
			}
		}else if(ps.size()>=2){
			std::string fname=ps[0];
			myParse.parse(ps[1].c_str());
			int idx=int(myParse());
			if(idx<0 || idx>fid2D.rows()){
				BLEXCEPTION(std::string(" Bad 'savefid' usage for \"")+inSqe+"\""+
				"\n the index counter should be >=0 and <npts2D")
			}
			fname=fname+itost(idx);

			if(canSave && toRun){
				if(is2D) curFID=fid2D.row(idx);
				double dt=1/Params.getD("sw")/2;
				if(isPtoP) dt=currentTime/2;
				plotterFID(curFID, fname, dt);
				if(saveWhat==10){
					matstream matout(fname, std::ios::out | std::ios::binary);
					matout.put("vdat", fid1D);
					matout.put("sw1", Params.getD("sw"));
					matout.close();
				}else if(saveWhat==20){
					std::fstream out(fname.c_str(), std::ios::binary | std::ios::out);
					if(out.fail()){
						std::cerr<<"Error: SaveFid \""<<Params.getS("filesave")<<"\" file cannot be opened for writing"<<std::endl;
						return;
					}
					out<<"sw1 ="<<Params.getS("sw")<<std::endl;
					out<<"npts1 ="<<fid1D.size()<<std::endl;
					double tmR, tmI;
					for(int p=0;p<fid1D.size();++p){
						tmR=Re(fid1D(p));
						tmI=Im(fid1D(p));
						out.write((char *)&tmR, sizeof(double));
						out.write((char *)&tmI, sizeof(double));
					}
					out.close();
				}else{
					plotterFID(curFID, fname, dt);
				}

				notSaved=false;
			}
		}
	}
}


#endif
