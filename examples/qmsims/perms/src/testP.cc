

#include "sequenceparse.h"

using namespace BlochLib;
using namespace std;


int main(int argc,char* argv[])
{
	int q=1, len;
	std::string fname;

	query_parameter(argc, argv, q++,"Enter File Name to parse ", fname);
	std::cout<<std::endl<<"Parameter File "<<fname<<std::endl;

	Parameters pset(fname);
	pset.addSection("spins");
	pset.addSection("params");
	pset.addSection("pulses");
	pset.addSection("permutate");

	SequenceParse gSeq(true);
	gSeq.parse(pset.section("pulses"));
	cout<<gSeq.myParse.getGlobalVar("wr")<<endl;

	Vector<string> subP=parse_param(pset.getParamS("typeIds", "permutate"),',');

	Parameters Pulses=pset.section("pulses");
	Vector<Parameters> subPulses(subP.size());
	for(int i=0;i<subP.size();++i){
		Pulses.addSection(subP[i]);
		subPulses[i]=Pulses.section(subP[i]);

		SequenceParse mySeq;
		mySeq.parse(paramStrip(Pulses.section(subP[i])));
		cout<<mySeq.pulses<<endl;
		cout<<"--------------------------------------------"<<endl;
	}
}


