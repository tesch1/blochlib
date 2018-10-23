
#include "RecoupleContent.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(int argc,char* argv[])
{
	std::string fname;
	query_parameter(argc, argv, 1,"Enter File Name to parse ", fname);
	std::cout<<std::endl<<"Parameter File \""<<fname<<"\""<<std::endl;

	Parameters pset(fname);
	pset.addSection("pulses");
	pset.addSection("params");
	pset.addSection("spins");

	SolidSys sys(pset.section("spins"));


	double wr=pset.getParamD("wr", "params");
	double beta=pset.getParamD("rotor", "params");

	int steps=10;
	int baseSym=pset.getParamI("baseSym", "pulses");
	int factorSym=pset.getParamI("factorSym", "pulses");
	int wrfactor =pset.getParamI("wrfactor", "pulses");

	double amp=baseSym*wr;

//	rcS.generateSubProps(sys,baseSym, factorSym, amp, wr, beta, steps);
//	rcS.generateProps(sys,baseSym, factorSym, amp, wr, beta, steps);

	RecoupleContent rc;
	rc.SubUnits.setParams(pset.getParamS("subUnitIds", "pulses"),
				  pset.getParamS("typeIds", "pulses"),
				  pset.getParamS("units", "pulses"));

	rc.SubUnits.debugflag=pset.getParamI("debuglevel", "params", false, 0);

	int order=pset.getParamI("permutations", "pulses");

	rc.generateTrains(order, baseSym, factorSym, amp);


}

