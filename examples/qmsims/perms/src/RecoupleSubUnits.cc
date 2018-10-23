//#include "fullbinarytree.h"
#include <string>
#include "blochlib.h"
//#include "RecoupleSequence.h"
#include "RecoupleSubUnits.h"
#include "spindex.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/******************* RECOUPLE SUB UNIT ***********/


RecoupleSubUnits::
  RecoupleSubUnits()
{
	subs_=parseSubUnits("o");
	types_=parseTypes("Cp");
	units=parseUnits("o");
	debugflag=0;
}


RecoupleSubUnits::
  RecoupleSubUnits(std::string intype,  std::string inunits, std::string insubs)
{
	setParams(intype, inunits, insubs);
	debugflag=0;
}

RecoupleSubUnits::
 RecoupleSubUnits(Vector<std::string> intype, Vector<std::string> inunits, Vector<std::string> insubs)
{
	setParams(intype, inunits, insubs);
	debugflag=0;
}

bool RecoupleSubUnits::checkType(std::string in)
{
	if(in != "Cp" && in != "C" && in!="Cpp" && in!="R" && in!="Rp" &&
	   in !="CppBar" &&  in !="CpBar" && in !="RBar" && in !="RpBar" && in!="CBar")
	{
		std::cerr<<std::endl<<"Error: RecoupleSubUnits"<<std::endl;
		std::cerr<<" there can only be C, Cp, Cpp, R, or Rp "<<std::endl;
		std::cerr<<"  CBar, CpBar, CppBar, RBar, or RpBar "<<std::endl;
		std::cerr<<" for 'type' your input " <<in<<std::endl;
		return false;
	}
	return true;
}

//looks at the input string "o,O,W"
//parses them into sparate chars, then checks to make sure there are
// no dups
Vector<char> RecoupleSubUnits::parseSubUnits(std::string in)
{
	Vector<std::string> tmv=parse_param(in, ',');
	Vector<char> outV;

	//check to make sure only one char is in the list
	for(int i=0;i<tmv.size();++i){
		if(tmv[i].size()>1){
			std::cerr<<std::endl<<"Error: RecoupleSubUnits::parseSubUnits"<<std::endl;
			std::cerr<<" there can only be 1 char per SubUnit id"<<std::endl;
			std::cerr<<" for 'type' your input " <<in<<std::endl;
			std::cerr<<" id tag #" <<i<<" id: "<<tmv<<std::endl;
			exit(1);
		}else{
			outV.push_back(tmv[i][0]);
		}
	}

	//check for dupes
	for(int i=0;i<outV.size();++i){
		for(int j=0;j<outV.size();++j){
			if(i!=j && outV[i]==outV[j]){
				std::cerr<<std::endl<<"Error: RecoupleSubUnits::parseSubUnits"<<std::endl;
				std::cerr<<" you have Duplicate entries in you id "<<std::endl;
				std::cerr<<" sub unit string, this cannot be " <<in<<std::endl;
				std::cerr<<" you input string: "<<in<<" the dupe: "<<outV[i]<<std::endl;
				exit(1);
			}
		}
	}
	return outV;
}

//will take in a string like "R, Rp, Cpp" (type units) and
//parses them into individual pieces
// it also checks that they are valid types
Vector<std::string> RecoupleSubUnits::parseTypes(std::string in)
{
	Vector<std::string> tmv=parse_param(in, ',');
	/*for(int i=0;i<tmv.size();++i){
		if(!checkType(tmv[i])){
			std::cerr<<std::endl<<"Error: RecoupleSubUnits::parseTypes"<<std::endl;
			std::cerr<<" you have an invalid type in your input string "<<std::endl;
			std::cerr<<" this cannot be--input::" <<in<<std::endl;
			exit(1);
		}
	}*/
	return tmv;
}

//private function that parses the input
// like "oowwee, ooOOEE, " into 'oowwee' 'ooOOEE'
Vector<std::string> RecoupleSubUnits::parseUnits(std::string inunit)
{	return parse_param(inunit, ',');	}


//chops a string like 'ooWW' in to 'o' 'o' 'W' 'W'
//and checks that each char is in the 'subs' vector
Vector<char> RecoupleSubUnits::chopUnits(std::string inunit)
{
	Vector<char> subunit;
	//chop up the string
	for(unsigned int i=0;i< inunit.length();++i){
		if(!isspace(inunit[i])) subunit.push_back(inunit[i]);
	}

	//check that they are in fact in the 'subs' vector
	for(int i=0;i< subunit.size();++i){
		bool got=false;
		for(int j=0;j<subs_.size();++j)	if(subunit[i]==subs_[j]) got=true;
		if(!got){
			std::cerr<<std::endl<<"Error: RecoupleSubUnits::parseUnits"<<std::endl;
			std::cerr<<" characters in your Units string are not'"<<std::endl;
			std::cerr<<" in your SubUnits character IDs "<<std::endl;
			std::cerr<<" Valid character::"<<subs_<<std::endl;
			std::cerr<<" your input::" <<subunit<<std::endl;
			exit(1);
		}
	}
	return subunit;
}


//set the internal parameters from input strings
void
	RecoupleSubUnits::setParams(
		std::string insubIDs,
		std::string intype,
		std::string inunits)
{
	subs_=parseSubUnits(insubIDs);
	types_=parseTypes(intype);
	if(subs_.size() != types_.size()){
		std::cerr<<std::endl<<"Error: RecoupleSubUnits::setParams"<<std::endl;
		std::cerr<<" There must the same number of Types defined"<<std::endl;
		std::cerr<<" for each SubUnit Id specfied"<<std::endl;
		std::cerr<<" Your subunit Ids: "<<insubIDs<<std::endl;
		std::cerr<<" Your Types Ids: "<<intype<<std::endl;
		exit(1);
	}
	units=parseUnits(inunits);
/*	if(units.size()<2){
		std::cerr<<std::endl<<"Error: RecoupleSubUnits::setParams"<<std::endl;
		std::cerr<<" There must be at least '2' Unit Strings in order"<<std::endl;
		std::cerr<<" to do anything"<<std::endl;
		std::cerr<<" Your Unit String: "<<inunits<<std::endl;
		exit(0);
	}*/
}

//set the internal parameters from input strings
void
	RecoupleSubUnits::setParams(
		Vector<std::string> insubIDs,
		Vector<std::string> intype,
		Vector<std::string> inunits
		)
{
	std::string tmm1="";
	std::string tmm2="";
	std::string tmm3="";
	for(int i=0;i<insubIDs.size();++i)
	{
		if(i!=insubIDs.size()-1) tmm1+=insubIDs[i]+",";
		else tmm1+=insubIDs[i];
	}
	for(int i=0;i<intype.size();++i)
	{
		if(i!=intype.size()-1) tmm2+=intype[i]+",";
		else tmm2+=intype[i];
	}
	for(int i=0;i<inunits.size();++i)
	{
		if(i!=inunits.size()-1) tmm3+=inunits[i]+",";
		else tmm3+=inunits[i];
	}
	setParams(tmm1, tmm2, tmm3);
}

//calculates the number of times step to satify the maxtstep condition
int RecoupleSubUnits::calcPropSteps(double maxtstep, double bt, double et)
{
	int st=1;
	st=int(std::ceil((et-bt)/maxtstep));
	if(st<=0){ st=1; }
	return st;
}

//this generates the Pulse Trains for
// each character in 'subUnits'
// it returns the rotor speed!!!
// which should be defined in the top of the pulses parameter section
double RecoupleSubUnits::generateTrains(Parameters &in)
{
	if(types_.size()<1) return 0.0;
	trains.resize(types_.size());
	subPropTime.resize(types_.size());

//find and evaulate any GLOABL variables
	SequenceParse gVar(true);
	gVar.parse(in);

//now look for each subunit..
	for(int i=0;i<types_.size();++i){

		if(!in.addSection(types_[i])){
			BLEXCEPTION(std::string(" the sub Unit pulse section ")+types_[i]+" cannot be found")
		}

	//find the sub section and parse it up to a
	//valid pulse sequence
		SequenceParse subSeq;
		subSeq.parse(in.section(types_[i]));

	//assigning the pulsese quence to the proper entry in the trains
		trains[i]=subSeq.pulses;
		subPropTime[i]=trains[i][trains[i].size()-1].t2-trains[i][0].t1;
		if(debugflag>0) std::cout<<"RecoupleSubUnits info--"<<trains[i]<<std::endl;
	}
	return gVar.myParse.getVar("wr");
}



//you need to sys.setPowderAngles before you call this
// AND You Must call 'generateTrains!!!'
// sys--> the Hamiltonian generator
// dtmin --> the desired dt max time step
void
	RecoupleSubUnits::
	 generateSubProps(SolidSys &sys,
	 	                double wr,
					 	double dtmin)
{
	if(sys.size()<1 || types_.size()<1) return;
	subProps.resize(types_.size());
	subProps.fill(sys.Fe());

//make sure the trains have been generated
	if(types_.size() != trains.size())
	{
		BLEXCEPTION(std::string("You Must generate the Reoupling Trains before calling")+
		"\n this function (i.e. a 'generateTrains(Parameters)");
	}

	//one prop for each type
	for(int i=0;i<trains.size();++i){
		hmatrix H;
		double dt, ton;

		int end=trains[i].size();

		for(int k=0;k<end;++k){
			int propSteps=calcPropSteps(dtmin, trains[i][k].t1,trains[i][k].t2);
			sys.setRotorAngles(1.0, trains[i][k].rotorangle);
			dt=(trains[i][k].t2-trains[i][k].t1)/double(propSteps);
			ton=trains[i][k].t1+dt/2.0;
			if(debugflag>=3)
				std::cout<<"calculating propogator for: "
					     <<trains[i][k]<<std::endl;
			for(int j=0;j<propSteps;++j){
				H=(sys.H(ton-dt/2.0, ton+dt/2.0, wr)+trains[i][k].Pulse(sys));
				subProps[i]*=Mexp(H, -complexi*dt*PI2);
				ton+=dt;
			}
		}
	}
}

//you need to sys.setPowderAngles before you call this
// this function will call 'generateSubProps'
void
	RecoupleSubUnits::
	 generateProps(SolidSys &sys,	// sys--> the Hamiltonian generator
				double wr,
				double dtmin	// the desired dt
				)
{
	generateSubProps(sys, wr, dtmin);
	generateProps();
}

//you need to sys.setPowderAngles before you call this
// this function will NOT call 'generateSubProps'
void RecoupleSubUnits::generateProps()
{
	if(units.size()<1) return;


	if(subs_.size() != types_.size()){
		std::cerr<<std::endl<<"Error: RecoupleSubUnits::generateProps"<<std::endl;
		std::cerr<<" There must the same number of Types defined"<<std::endl;
		std::cerr<<" for each SubUnit Id specfied"<<std::endl;
		std::cerr<<" Your subunit Ids: "<<subs_<<std::endl;
		std::cerr<<" Your Types Ids: "<<types_<<std::endl;
		exit(1);
	}

	if(subProps.size() != subs_.size() || subProps.size() != types_.size()){
		std::cerr<<std::endl<<"Error: RecoupleSubUnits::generateProps"<<std::endl;
		std::cerr<<" must call 'generateSubProps' first"<<std::endl;
		exit(1);
	}

	Props.resize(units.size());
	propTime.resize(units.size());

	for(int i=0;i<units.size();++i){ //go through each 'unit'
//cout<<endl<<"UNITS: "<<units.size()<<"|"<<units<<"| index: "<<i<<endl;
		Vector<char> subunit=chopUnits(units[i]);
		if(debugflag>0)
				std::cout<<"RecoupleSubUnits info--"<<"Total Unit: "<<subunit<<std::endl;
		for(int j=0;j<subunit.size();++j){ //go through each 'sub unit'
			if(debugflag>=2)
				std::cout<<"RecoupleSubUnits info----"<<"Props: "<<j<<" which: "<<subunit[j]<<endl;
			if(j==0){	//the inital U
				for(int k=0;k<subs_.size();++k){
					if(subs_[k]==subunit[j]){
						Props[i]=subProps[k];
						propTime[i]=subPropTime[k];
						break;
					}
				}
			}else{
				for(int k=0;k<subs_.size();++k){
					if(subs_[k]==subunit[j]){
						Props[i]*=subProps[k];
						propTime[i]+=subPropTime[k];
						break;
					}
				}
			}
		}
	}
}




