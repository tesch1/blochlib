

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;


int main(int argc,char* argv[]){

	Vector<std::string> spinl;
	spinl.push_back(" numspin 1\n");
	spinl.push_back("T 2H 0\n");
	spinl.push_back("Q 100e3 0 0\n");

	SolidSys mysys(spinl);
	int q=1;
	int npts=512;
	query_parameter(argc,argv,q++, "Enter Collection points: ", npts);
	double wr=200;
	double beta=54.7;
	query_parameter(argc,argv,q++, "Enter rotor speed: ", wr);
	query_parameter(argc,argv,q++, "Enter rotor angle: ", beta);
	int thd=144, phs=89;
	query_parameter(argc,argv,q++, "theta points: ", thd);
	query_parameter(argc,argv,q++, "phi points: ", phs);
	double bf=0;
	query_parameter(argc,argv,q++, "B field: ", bf);
	double wq=0;
	query_parameter(argc,argv,q++, "Quad Freq (Hz): ", wq);
	mysys.qua[0].Q(wq);
	double sw=600000;
	query_parameter(argc,argv,q++, "Sweep width (Hz): ", sw);
	double lb=50;
	query_parameter(argc,argv,q++, "Line Broadening (Hz): ", lb);
	mysys.SetBfield(bf);
	//mysys.qua[0].order=2;
	cout<<mysys<<endl;

	powder zcw(powder::alderman, thd, phs);

	matrix detect=mysys.Fp();
	hmatrix roeq=mysys.Fx();

	Vector<complex> fid(npts, 0);


	oneFID<SolidSys> myfid(mysys, npts, sw, wr, beta);
	fid=myfid.FID(zcw, roeq, detect);
	complex sm=0;
	sm=sum(fid(Range(npts/2, npts)))/double(npts/2.0);
	fid-=sm;
	plotterFID(fid, "fid", 1.0/myfid.sw(),lb, 1,1);
}

