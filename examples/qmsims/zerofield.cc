#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

class ZeroFieldDipole :
	public SolidSys
{
	public:
		double D;
		ZeroFieldDipole():SolidSys(), D(0){}
		ZeroFieldDipole(double d):SolidSys(), D(d){}
		ZeroFieldDipole(SolidSys &sys, double d):SolidSys(sys),  D(d){}

		hmatrix Hamiltonian(double t1, double t2, double wr=0.0)
		{
			 double tmid=((t2-t1)/2.0+t1)*wr*PI2;
             return D*(
				 A2(tmid,theRotations.beta,-2)*T22(*this, 0,1)
				 -A2(tmid,theRotations.beta,-1)*T21(*this, 0,1)
				 +A2(tmid,theRotations.beta,0)*T20(*this, 0,1)
				 -A2(tmid,theRotations.beta,1)*T2m1(*this, 0,1)
				 -A2(tmid,theRotations.beta,2)*T2m2(*this, 0,1));
		}

		hmatrix Hamiltonian(double wr, double rot, double alpha, double beta, double t1, double t2)
		{
			return sqrt(6.0)*D*(
				 A2(alpha,beta,-2)*T22(*this, 0,1)
				 -A2(alpha,beta,-1)*T21(*this, 0,1)
				 +A2(alpha,beta,0)*T20(*this, 0,1)
				 -A2(alpha,beta,1)*T2m1(*this, 0,1)
				 +A2(alpha,beta,2)*T2m2(*this, 0,1));
		}


};


int main(int argc,char* argv[]){

	Vector<std::string> spinl;
	spinl.push_back(" numspin 2\n");
	spinl.push_back("T 1H 0\n");
	spinl.push_back("T 1H 1\n");
	spinl.push_back("D 2000 0 1\n");

	SolidSys mysys(spinl);
	int q=1;
	double d=512;
	query_parameter(argc,argv,q++, "Enter dipole coupling: ", d);
	ZeroFieldDipole myD(mysys, d);
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
	matrix roeq=mysys.Fx();

	Vector<complex> fid(npts, 0);

	oneFID<ZeroFieldDipole> myfid(myD, npts, sw, wr, beta);
	fid=myfid.FID(zcw, roeq, detect);

	plotterFID(fid, "fid", 1.0/myfid.sw());
}

