

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
 	sequencerun.h -->
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

#ifndef __sequencerun_cc__
#define __sequencerun_cc__ 1

#include "sequencerun.h"

#include <signal.h>

subSequence::subSequence():
	haveParsed(false), props(NULL),hold(false), repeat(1),
	lastPow(-1), lastT(-1),cycler("Fe"), docycler(false)
{}

subSequence::subSequence(const Vector<std::string> &indata, Propogation &inp):
	data(indata),haveParsed(false), props(&inp),
	hold(false), repeat(1), lastPow(-1), lastT(-1),
cycler("Fe"), docycler(false)
{
}
subSequence::~subSequence(){	props=NULL;	}

subSequence::subSequence(const subSequence &cp)
{
	*this=cp;
}

void subSequence::operator=(const subSequence &cp)
{
	if(this==&cp) return;
	haveParsed=cp.haveParsed;
	data=cp.data;
	props=cp.props;
	hold=cp.hold;
	repeat=cp.repeat;
	lastPow=cp.lastPow;
	lastT=cp.lastT;
	cycler=cp.cycler;
	U=cp.U;
	holdU=cp.holdU;
	myPulses=cp.myPulses;
	docycler=cp.docycler;
}

void subSequence::propogator(int powpt, int pt, double startT)
{
	if(haveParsed && lastPow==powpt && startT==lastT){	return ;	}

	lastPow=powpt; lastT=startT;


	if(data.empty()){
		U=props->curSys->Fe();
		return;
	}

	if(!haveParsed){
		SequenceParse myP;
		myP.parse(data);
		myPulses=readPulses(myP.sequence, &myP.myParse);
		haveParsed=true;
	}

	if(startT!=0){
		if(myPulses.size()>0){
			double dt=myPulses[0].t2-myPulses[0].t1;
			myPulses[0].t1=startT+dt;
			myPulses[0].t2=myPulses[0].t1+dt;
			for(int i=1;i<myPulses.size();++i){
				dt=myPulses[i].t2-myPulses[i].t1;
				myPulses[i].t1=myPulses[i-1].t2;
				myPulses[i].t2=myPulses[i].t1+dt;
			}
		}
	}
	/*
		for(int i=0;i<myPulses.size();++i){
			if(i==0) myPulses[i].print(std::cout, true);
			else myPulses[i].print(std::cout);
		
	}*/
	props->pdata=myPulses;
	U=props->propogator(powpt);

}

void subSequence::propogate(matrixs &ro, int powpt, int fidpt)
{

	if(docycler){
		HamiltonianGen myGen;
		if(cycler!="Fe" && cycler!="Ie" && cycler!=""){
			matrixs trac=myGen.Hamiltonian((*(props->curSys)), cycler,
									 props->curPow->theta(powpt),
									 props->curPow->phi(powpt),
									 props->curPow->gamma(powpt));
			scomplex norm=trace(trac, adjoint(trac));
			ro=trace(ro,adjoint(trac))*trac/norm;
			//std::cout<<ro<<std::endl;
		}
		return;
	}

	if(fidpt==0) holdU=U;

	if(hold){
		if(fidpt==0){
			for(int i=0;i<repeat;++i){
				ro.prop(holdU);
			}
		}


	}else{
		for(int i=0;i<repeat;++i){
			ro.prop(holdU);
		}
	}
}

Vector<PulseData> &subSequence::parse()
{

	if(!haveParsed){
		SequenceParse myP;
		myP.parse(data);
		myPulses=readPulses(myP.sequence, &myP.myParse);
		haveParsed=true;
	}

	return myPulses;
}

std::ostream &operator<<(std::ostream &otr, subSequence oo){
	//otr<<"**Pulse Data**"<<std::endl;
	//otr<<"\t"<<oo.parse();
	otr<<"**Current Propogator**"<<std::endl;
	otr<<"\t"<<oo.U<<std::endl;
	return otr;
}



/*************
*********************** SubSequence FUnction *********
************************/

/**************** Sequence Parser bits *************/
SequenceRun::SequenceRun():
	//SequenceParse(),
	amp(0), phase(0), offset(0),currentTime(0),
	on("1H"), usesys("default"), usepow("default"), cycler("Ie"),
	isPtoP(false), is2D(false), powpt(0),pt2D(-1), toRun(false),
	dumpOnDie(false)
{
	myProp.Systems=&Systems;
	myProp.Powders=&Powders;
	myProp.Params=&Params;
	numProp=0;
	PropList.resize(1000);
}


SequenceRun::SequenceRun(Parameters &pset) :
	//SequenceParse(),
	amp(0), phase(0), offset(0),currentTime(0),
	on("1H"), usesys("default"), usepow("default"), cycler("Ie"),
	isPtoP(false), is2D(false),  powpt(0),pt2D(-1),toRun(false),
	dumpOnDie(false)
{
	myProp.Systems=&Systems;
	myProp.Powders=&Powders;
	myProp.Params=&Params;
	SequenceRun::parsePulse(pset);
	numProp=0;
	PropList.resize(1000);
}

SequenceRun::SequenceRun(Vector<std::string> &pset):
	//SequenceParse(),
	amp(0), phase(0),  offset(0),currentTime(0),
	on("1H"), usesys("default"), usepow("default"), cycler("Ie"),
	isPtoP(false), is2D(false),  powpt(0),pt2D(-1), toRun(false),
	dumpOnDie(false)
{
	myProp.Systems=&Systems;
	myProp.Powders=&Powders;
	myProp.Params=&Params;
	numProp=0;
	PropList.resize(1000);
}




//this is the master decider function
// and snaps it to the proper above function
bool SequenceRun::decide(std::string inSqe)
{
	//need to override these in the SequenceParse as
	// here they acctually calculate propogators
	//if(inSqe.find("pulse(")<inSqe.size()){
	//	 doPulse(inSqe);
	//}else if(inSqe.find("delay(")<inSqe.size()){
	//	 doDelay(inSqe);
	//}else if(inSqe.find("on(")<inSqe.size()){
	//	 setOn(inSqe);
	try{
		if(inSqe.find("pulse(")<inSqe.size() || inSqe.find("delay(")<inSqe.size()){
			doPulseData(inSqe);
		}else if(inSqe.find("amplitude(")<inSqe.size()){
			setAmplitude(inSqe);
		}else if(inSqe.find("offset(")<inSqe.size()){
			setOffset(inSqe);
		}else if(inSqe.find("spinsys(")<inSqe.size()){
			setSystem(inSqe);
		}else if(inSqe.find("powder(")<inSqe.size()){
			setPowder(inSqe);
		}else if(inSqe.find("cycler(")<inSqe.size()){
			doCycler(inSqe);
		}else if(inSqe.find("dumpState()")<inSqe.size()){
			dumpState();
		}
	
		else if(inSqe.find("=")<inSqe.size()
			&& inSqe.find("==")!=inSqe.find("=")
			&& inSqe.find("!=")!=inSqe.find("=")-1 ){
			addGVar(inSqe);
		}else if(inSqe.find("ptop()")<inSqe.size()){
			addPtop(inSqe);
		}else if(inSqe.find("show()")<inSqe.size()){
			Shows(inSqe);
		}else if(inSqe.find("2D()")<inSqe.size()){
			to2D(inSqe);
		}else if(inSqe.find("ro(")<inSqe.size()){
			setRo(inSqe);
		}else if(inSqe.find("detect(")<inSqe.size()){
			setDetect(inSqe);
		}else if(inSqe.find("U(")<inSqe.size()){
			showProps(inSqe);
		}else if(inSqe.find("use(")<inSqe.size()){
			showUsed(inSqe);
		}else if(inSqe.find("reset(")<inSqe.size()){
			ResetPars(inSqe);
		}else if(inSqe.find("alterSys(")<inSqe.size()){
			doAlterSys(inSqe);
		}else if(inSqe.find("savefid(")<inSqe.size()){
			saveFID(inSqe);
		}else if(inSqe.find("savefidmatlab(")<inSqe.size()){
			saveFID(inSqe);
		}else if(inSqe.find("savefidtext(")<inSqe.size()){
			saveFID(inSqe);
		}else if(inSqe.find("savefidbinary(")<inSqe.size()){
			saveFID(inSqe);
		}else if(inSqe.find("fid(")<inSqe.size()){
			doFID(inSqe);
		}else{
			decideSeP(inSqe);
		}
	}catch(BL_exception e){
		e.print(std::cerr);
		BLEXCEPTION("could not evaluate art of script");
	}
	return true;
}

void SequenceRun::setParams(Parameters &pset)
{
	setParams(pset.section(""));
}

void SequenceRun::setParams(const Vector<std::string> &pset)
{
	Params.parse(pset);
}

void SequenceRun::setSystems(Parameters &pset)
{
	Systems.parse(pset);
}

void SequenceRun::setSystems(const Vector<std::string> &pset)
{
	Parameters pt(pset);
	Systems.parse(pt);
}

void SequenceRun::setPowders(Parameters &pset)
{
	Powders.parse(pset);
}

void SequenceRun::setPowders(const Vector<std::string> &pset)
{
	Parameters pt(pset);
	Powders.parse(pt);
}


//Paraser the pulses ection into the 'global' items and the
// subPulse sections
void SequenceRun::parsePulse(Parameters &pset)
{
	std::string base="sub";
	mainProc=paramStrip(pset.section(""));
	data=mainProc;
	//mainProc.print(std::cout, "\n");
//if there is only one section 'named' powders
// then there is no list, only one powder....

	int maxFit=pset.getParamI("numsub","",false, 100);
	int numPows=0;
//count the number of params present
	while(pset.addSection(base+itost(numPows+1)) && numPows<=maxFit )
	{	numPows++;		}

	if(numPows==0)	return;

//add a parameter to our master list
	int i=0;
	while(i<numPows)
	{
		if(!pset.section(base+itost(i+1)).empty()){
			subProc.insert(
				std::pair<std::string, subSequence >(
					base+itost(i+1), subSequence(pset.section(base+itost(i+1)), myProp)
				)
			);
		}
		++i;
	}
}


int SequenceRun::run()
{
	fid1D.resize(Params.getI("npts1D"));
	fid1D=0;
	fid2D.resize(Params.getI("npts2D"), Params.getI("npts1D"));
	fid2D.fill(scomplex(0.,0.));
	notSaved=true;
	//if(Params.getI("npts2D")>1){	is2D=true; pt2D=0;	}

	//for(int pp=0;pp<myProp.setPowder(usepow);
	myProp.setSystem(usesys);
	myProp.setPowder(usepow);
	HamiltonianGen myGen;
	roeq=myGen.Hamiltonian(*(myProp.curSys), Params.getS("roeq"));
	if(MPIworld.master() && myProp.curSys->isDiagonal())
	  std::cout<<"Diagonal Hamiltonian"<<std::endl;
	if(MPIworld.master() )
	  std::cout<<"-----(theta, phi  ,gamma)----point/total---"<<std::endl;

	canSave=false;


	//do a master-slave parallel model
	if(MPIworld.parallel()){
		int done=-1, curP;
		powpt=0;
		ro=roeq;
		numProp=0;
		U=myProp.curSys->Fe();
		currentTime=0;
		if(MPIworld.master()){

			//the initialization
			for(int proc=1;proc<MPIworld.size();++proc){
				if(powpt>=(myProp.curPow->size()-1)) break;
				std::string powPtPrint="    ";
				std::cout<<"     ("<<dbtost(myProp.curPow->theta(powpt), "%1.3f")<<","
						 <<dbtost(myProp.curPow->phi(powpt), "%1.3f")<<", "
						 <<dbtost(myProp.curPow->gamma(powpt), "%1.3f")<<")";

				powPtPrint=itost(powpt, 5)+"/"+itost(myProp.curPow->size(),5)+"    ";
				std::cout<<"   "<<powPtPrint<<"\r";
				MPIworld.put(powpt, proc); //std::cout.flush();
				powpt++;
			}

			//the main loop
			//the last powder point must be done on the master
			// as we need to savefid from the master
			while(powpt<(myProp.curPow->size()-1)){
				int curProc=MPIworld.getAny(curP);
				MPIworld.put(powpt, curProc);
				std::string powPtPrint="    ";
				std::cout<<"     ("<<dbtost(myProp.curPow->theta(powpt), "%1.3f")<<","
						 <<dbtost(myProp.curPow->phi(powpt), "%1.3f")<<", "
						 <<dbtost(myProp.curPow->gamma(powpt), "%1.3f")<<")";

				powPtPrint=itost(powpt, 5)+"/"+itost(myProp.curPow->size(),5)+"    ";
				std::cout<<"   "<<powPtPrint<<"\r"; //std::cout.flush();
				++powpt;
			}

			//send termination
			for(int proc=1;proc<MPIworld.size();++proc) MPIworld.put(done, proc);

			//set the save flag when we are done
			canSave=true;
			//do the final powder point
			powpt++;
		//set the Bfield all the systems
			myProp.Systems->setBfield(Params.getD("Bfield"));
			parse(mainProc); //do it and hopefully save
			//save it if we have not yet
			if(notSaved){
				saveFID(std::string("savefid()"));
			}


		}else{
			while(1){
				MPIworld.get(powpt, 0);
				if(powpt==done) break;
			//set the Bfield all the systems
				myProp.Systems->setBfield(Params.getD("Bfield"));
				ro=roeq;
				numProp=0;
				U=myProp.curSys->Fe();
				currentTime=0;
				parse(mainProc);
				MPIworld.put(powpt, 0);
			}
			canSave=true;
			//isue a save command so that the 'reduce' will work
			if(notSaved){
				saveFID("");
			}
		}
	}else{ //serial mode
	//set the Bfield all the systems
		myProp.Systems->setBfield(Params.getD("Bfield"));
		for(powpt=0;powpt<myProp.curPow->size();++powpt){
			ro=roeq;
			numProp=0;
			U=myProp.curSys->Fe();
			currentTime=0;
			std::string powPtPrint="    ";
			std::cout<<"     ("<<dbtost(myProp.curPow->theta(powpt), "%1.3f")<<","
					 <<dbtost(myProp.curPow->phi(powpt), "%1.3f")<<", "
					 <<dbtost(myProp.curPow->gamma(powpt), "%1.3f")<<")";

			powPtPrint=itost(powpt, 5)+"/"+itost(myProp.curPow->size(),5)+"    ";
			std::cout<<"   "<<powPtPrint<<"\r";

			//if we get a 'savefid comman' we only want to save after the powder is done
			if(powpt==myProp.curPow->size()-1) canSave=true;
			else canSave=false;

			parse(mainProc);
			//std::cout.flush();

		}
		//save it if we have not yet
		if(notSaved){
			saveFID(std::string("savefid()"));
		}
	}
	return 1;
}


void SequenceRun::printSeq(std::ostream &oo)
{

	bool oldRun=toRun;
	toRun=false;
	currentTime=0;

	//first parse the main to see if it sets any vars
	oo<<"******Main Pulse Process "<<std::endl;
	PulseData::printHead(oo);
	if(!mainProc.empty()){
		parse(mainProc);
	}

	Vector<PulseData> tmPPulse=readPulses(sequence);
	oo<<tmPPulse<<std::endl;


	toRun=oldRun;
}

//The sig die procedure
void SequenceRun::dumpState()
{
	//Vector<PulseData> tmPPulse=readPulses(sequence);
	std::cerr<<(this->Systems)<<std::endl;
	std::cerr<<(this->Powders)<<std::endl;
	std::cerr<<(this->Params)<<std::endl;
std::cerr<<sequence<<std::endl;

std::cerr<<"**Powder Angle:::"<<"point: "<<powpt<<"/"<<myProp.curPow->size()<<" ("<<dbtost(myProp.curPow->theta(powpt), "%1.3f")<<","
	<<dbtost(myProp.curPow->phi(powpt), "%1.3f")<<", "
	<<dbtost(myProp.curPow->gamma(powpt), "%1.3f")<<")"<<std::endl;
	std::cerr<<"**Initial density operator::: "<<std::endl;
	std::cerr<<roeq<<std::endl;
	std::cerr<<"**Current Propogators in sequences:::"<<std::endl;
	showProps("U()");
std::cerr<<"**Current density opertor (rho):::"<<std::endl;
	std::cerr<<ro<<std::endl;
}

std::ostream &operator<<(std::ostream &oo, SequenceRun &out)
{
	oo<<(out.Systems)<<std::endl;
	oo<<(out.Powders)<<std::endl;
	oo<<(out.Params)<<std::endl;
	out.printSeq(oo);

	return oo;
}

#endif
