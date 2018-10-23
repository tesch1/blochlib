

#include "magfitter.h"

/* a program to test the  MagFitter class */
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

int main(int argc, char **argv)
{
	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);

	Parameters pset(fn);

	pset.addSection("fitter");
	std::string CoiSec=pset.getParamS("section");
	pset.addSection(CoiSec);

	MagFitter mf(pset.section("fitter"), pset.section(CoiSec));
	cout<<mf<<endl;
		mf.coilParams().print(cout);

	double *vals; vals =new double[mf.size()];
	for(int i=0;i<5;++i){
		for(int j=0;j<mf.size();++j){
			vals[j]=Rand()*5.0;
		}
		mf.setCoilParams(vals);
		cout<<mf<<endl;
		mf.coilParams().print(cout);
	}
}

