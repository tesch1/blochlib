

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*
	the second order quadrople EXPLICITLY calculated in the
	'high field' approx to second order


*/

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


class mCos
{
	double a, b;
	public:
		mCos():a(0), b(0){}
		mCos(double ai, double bi): a(ai), b(bi){}
		double operator()(double t)
		{	return cos(PI2*a*t+PI2*b);	}
};

class mSin
{
	double a, b;
	public:
		mSin():a(0), b(0){}
		mSin(double ai, double bi): a(ai), b(bi){}
		double operator()(double t)
		{	return sin(PI2*a*t+PI2*b);	}
};

template<class M1, class M2>
MUL_STRUCTURE(M1, M2) commutator(M1 mat1, M2 mat2)
{
	return mat1*mat2-mat2*mat1;
}


int main(int argc,char* argv[]){

#ifdef HAVE_MPI
	int size;
	int rank;
	int err;
	//MPI_Status *stat;
	err=MPI_Init(&argc, &argv);
	err=MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* Get my rank   */
	err=MPI_Comm_size(MPI_COMM_WORLD, &size);   /* Get the total */
#endif
	int q=1;
	std::string fn;
	query_parameter(argc,argv,q++, "Enter file to parse: ", fn);
	Parameters pset(fn);
	pset.AddSection("spins");
	pset.AddSection("pars");

	SolidSys mysys(pset.Section("spins"));
	int npts=512;
	double wr=200;
	double  wq, wo, J;
	double thetap=144, phip=89;
	int sample;
	wr=pset.GetParamD("wr","pars");
	npts=pset.GetParamI("npts","pars");
	sample=pset.GetParamI("sample","pars");
	double tm=pset.GetParamD("timestep","pars");

	std::string avetype=pset.GetParamS("avetype","pars");
	powder mypow;
	int gammast=0;
	if(avetype=="liquid"){
		thetap=pset.GetParamD("theta","pars");
		thetap*=Pi/180.0;
		phip=pset.GetParamD("phi","pars");
		phip*=Pi/180.0;
		mypow=powder(thetap, phip);
	}else{
		int thetast=pset.GetParamI("thetastep","pars");
		int phist=pset.GetParamI("phistep","pars");
		gammast=pset.GetParamI("gammastep","pars");
		mypow=powder(avetype, thetast, phist, gammast);
	}

	wq=pset.GetParamD("wq","pars")*1.e6*PI2;
	wo=pset.GetParamD("wo","pars")*1.e6*PI2;
	J=pset.GetParamD("J","pars")*PI2;
	double etap=pset.GetParamD("eta","pars")*PI2;
	int dhq=pset.GetParamI("secondorder","pars");
	std::string fname=pset.GetParamS("outfile","pars");

	int wrsteps=pset.GetParamI("wrsteps", "pars");
	double wrstep=pset.GetParamD("wrstep", "pars");

	int quiet=pset.GetParamI("verbose", "pars", false, 0);
	cout<<mysys<<endl;

	matrix detect=mysys.Ip(0);
	hmatrix roeq=mysys.Ix(0);

	Vector<complex>  A(5,0);
	Vector<double> xo(5,0);


	xo(0)=tm;
	hmatrix H=mysys.Fe();
	matrix U=mysys.Fe();
	complex A20, A21, A21m, A22, A22m;

	matrix T20=mysys.Iz(1)*mysys.Iz(1),

	T21T2m1=commutator(mysys.Iz(1)*mysys.Imi(1)+ mysys.Imi(1)*mysys.Iz(1),mysys.Iz(1)*mysys.Ip(1) +mysys.Ip(1)*mysys.Iz(1)),

	T22T2m2=commutator(mysys.Imi(1)*mysys.Imi(1),mysys.Ip(1)*mysys.Ip(1)),
	//[IZ*IZ,IZI+ + I+IZ]
	T20T21=commutator(T20, mysys.Iz(1)*mysys.Ip(1)+mysys.Ip(1)*mysys.Iz(1)),
	//[IZ*IZ,IZI- + I-IZ]
	T20T2m1=commutator(T20, mysys.Iz(1)*mysys.Imi(1)+mysys.Imi(1)*mysys.Iz(1)),
	//[IZ*IZ,(I+)^2]
	T20T22=commutator(T20, mysys.Ip(1)*mysys.Ip(1)),
	//[IZ*IZ,(I-)^2]
	T20T2m2=commutator(T20, mysys.Imi(1)*mysys.Imi(1));

	if(wrsteps==0) wrsteps=1;
	double wrstart=wr;
	std::string fbase=fname;

#ifdef HAVE_MPI
	if(wrsteps%size !=0)
	{
		std::cerr<<"ERROR: ***"<<endl;
		std::cerr<<" Number of rotor steps must be a "<<endl
		<<" multiple of the number of processors..."<<endl
		<<" THe end.."<<endl;
		exit(0);
	}

	int start=rank*int(wrsteps/size);
	int end=(rank+1)*int(wrsteps/size);
	for(int i=start;i<end;++i)
	{
		fname=fbase+itost(i);
		wr=wrstart+double(i)*wrstep;
		if(quiet>=1) cout<<"PROC:: "<<i<<" "<<fname<<" "<<wr<<endl;

#else
	for(int wrss=0; wrss<wrsteps;wrss++)
	{
		wr=wrstart+double(wrss)*wrstep;
		if(wrsteps==1) fname=fbase;
		else fname=fbase+itost(wrss);
#endif
		Vector<complex> fid(npts,0);
		powder::iterator myit(mypow);
		while(myit){

			mCos x1(wr, myit.gamma());
			mSin x2(wr, myit.gamma());
			mCos x3(2.0*wr, 2.0*myit.gamma());
			mSin x4(2.0*wr, 2.0*myit.gamma());

			Integrate<mCos> intC1(x1);
			Integrate<mSin> intS1(x2);
			Integrate<mCos> intC2(x3);
			Integrate<mSin> intS2(x4);

			double Co=1./8.*sqrt(3.)/sqrt(2.)*(3.*cos(myit.theta())*cos(myit.theta())-1.+etap*sin(myit.theta())*sin(myit.theta())*cos(2.*myit.phi()));
			double C1=1./4.*sqrt(3./2.)*sin(2.*myit.theta())*(1.-1./3.*etap*cos(2.*myit.phi()));
			double C2=1./8.*sqrt(3./2.)*(sin(myit.theta())*sin(myit.theta())+etap/3.*(cos(thetap)*cos(myit.theta())+1.)*cos(2.*myit.phi()));
			double S1=1./2./sqrt(6.)*etap*sin(myit.theta())*sin(2.*myit.phi());
			double S2=1./4./sqrt(6.)*etap*cos(myit.theta())*sin(2.*myit.phi());

			for(double k=0;k<npts;++k)
			{
				if(quiet>1) cout<<"wr: "<<wr<<"/"<<wrstart+double(wrsteps)*wrstep
				    <<" Powder: ("<<myit.theta()<<","<<myit.phi()
				    <<","<<myit.gamma()<<") "<<"on Pt: "<<k<<"\r"; cout.flush();
				for(double kk=1;kk<=sample;++kk)
				{

					xo(1)=intC1((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
					xo(2)=intS1((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
					xo(3)=intC2((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
					xo(4)=intS2((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);

					A20=-C1*xo(1)*sqrt(3./2.)*2.*sqrt(2.)/3.-
					S1*xo(2)*sqrt(3./2.)*2.*sqrt(2.)/3.+
					C2*xo(3)*sqrt(6.)*2./3.-
					S2*xo(4)*sqrt(6.)*2./3.;

					A21=Co*2.*sqrt(2.)/3.+
					C1*(-1./3.*xo(1)-complexi/sqrt(3)*xo(2))+
					S1*(-1/3*xo(2)+complexi*1/sqrt(3.)*xo(1))-
					C2*(2.*sqrt(2.)/3.*xo(3)-2.*complexi*sqrt(2./3.)*xo(4))+
					S2*(2.*sqrt(2.)/3.*xo(4)+2.*complexi*sqrt(2./3.)*xo(3));

					A21m=-Co*2.*sqrt(2.)/3.-
					C1*(-1./3.*xo(1)+complexi/sqrt(3.)*xo(2))-
					S1*(-1./3.*xo(2)-complexi*1./sqrt(3.)*xo(1))+
					C2*(2.*sqrt(2.)/3.*xo(3)+2.*complexi*sqrt(2./3.)*xo(4))-
					S2*(2.*sqrt(2.)/3.*xo(4)-2.*complexi*sqrt(2./3.)*xo(3));

					A22=Co*2./3.+
					C1*(sqrt(2.)/3.*xo(1)-complexi*sqrt(2./3.)*xo(2))+
					S1*(sqrt(2.)/3.*xo(2)+complexi*sqrt(2./3.)*xo(1))+
					C2*(4./3.*xo(3)-2*complexi*1/sqrt(3.)*xo(4))-
					S2*(4./3.*xo(4)+2.*complexi*1/sqrt(3.)*xo(3));

					A22m=Co*2./3.+
					C1*(sqrt(2.)/3.*xo(1)+complexi*sqrt(2./3.)*xo(2.))+
					S1*(sqrt(2.)/3*xo(2.)-complexi*sqrt(2./3.)*xo(1))+
					C2*(4./3.*xo(3)+2*complexi*1./sqrt(3.)*xo(4))-
					S2*(4./3.*xo(4)-2.*complexi*1./sqrt(3.)*xo(3));


					//H=	wo*tm*mysys.Iz(1)+
					if(dhq)
					{
						H=	wq*A20*(3.*T20)
						-3./2.*wq*wq/wo*A21*A21m*T21T2m1
						-3./2.*wq*wq/wo*A22*A22m*T22T2m2
						-sqrt(3./2.)*wq*wq/wo*3.*
						(
							T20T21*A20*A21 -
							A21m*A20*T20T2m1 -
							A20*A22/2.0*T20T22 +
							A20*A22m/2.0*T20T2m2
						)
						+J*tm*mysys.Iz(0)*mysys.Iz(1);
					}else{
						H=	wq*A20*(3.*T20)
						+J*tm*mysys.Iz(0)*mysys.Iz(1);
					}

					U=Mexp(H, -complexi)*U;
				}

				//printTime(int(k+1));


				fid(int(k))+=myit.weight()*trace(detect*U*roeq*adjoint(U));

			}
			++myit;
		}
		plotterFID(fid, fname, tm*sample);
	}

#ifdef HAVE_MPI
	MPI::Finalize();
#endif

}

