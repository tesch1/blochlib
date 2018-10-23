

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/* a simple 'hardcoded' 90 pulse on the spins system below */


int main(int argc,char* argv[]){

	Vector<std::string> spinl;
	spinl.push_back(" numspin 2\n");
	spinl.push_back("T 1H 0\n");
	spinl.push_back("T 1H 1\n");
	spinl.push_back("C 1243 1548 0 0\n");
	spinl.push_back("C -4561 3215 0 1\n");
	//spinl.push_back("D -2120 0 1\n");

	SolidSys mysys(spinl);
	int q=1;
	int npts=512;
	query_parameter(argc,argv,q++, "Enter Collection points: ", npts);
	double sw=20000;
	query_parameter(argc,argv,q++, "Enter sweep width(Hz): ", sw);
	double wr=200;
	double beta=54.7;
	query_parameter(argc,argv,q++, "Enter rotor speed: ", wr);
	query_parameter(argc,argv,q++, "Enter rotor angle: ", beta);
	int thd=144, phs=89;
	query_parameter(argc,argv,q++, "theta points: ", thd);
	query_parameter(argc,argv,q++, "phi points: ", phs);

	cout<<mysys<<endl;

	powder zcw(powder::zcw, thd, phs);

	matrix detect=mysys.Fp();
	hmatrix roeq=mysys.Fx()+mysys.Fz();

	Vector<complex> fid(npts, 0);

	oneFID<SolidSys> myfid(mysys, npts, sw, wr, beta);
	fid=myfid.FID(zcw,roeq, detect);

	plotterFID(fid, "fid", 1.0/myfid.sw());
}

