

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*function x=cosI(x,a,b);
m=length(x);

for l=1:m;
x(1,l)=cos(2*pi*a*x(1,l)+2*pi*b);

end;
return x
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

	Vector<std::string> spinl;
	spinl.push_back(" numspin 2\n");
	spinl.push_back("T 13C 0\n");
	spinl.push_back("T 14N 1\n");

	SYS mysys(spinl);
	int q=1;
	int npts=512;
	double wr=200;
	double wrs=54.7, wq, wo, J;
	double thetap=144, phip=89;
	int sample, Nwr;
	query_parameter(argc,argv,q++, "Enter initial rotor speed: ", wr);
#ifdef HAVE_MPI
	query_parameter(argc,argv,q++, "Enter rotor speed step size: ", wrs);
	query_parameter(argc,argv,q++, "Enter number of rotor steps: ", Nwr);

#endif
	query_parameter(argc,argv,q++, "Enter Collection points: ", npts);
	query_parameter(argc,argv,q++, "Samples..: ", sample);
	query_parameter(argc,argv,q++, "Theta angle (deg): ", thetap);
	thetap*=Pi/180.0;
	query_parameter(argc,argv,q++, "phi angle(deg): ", phip);
	phip*=Pi/180.0;
	query_parameter(argc,argv,q++, "Enter wq (MHz): ", wq);
	wq*=1.e6*PI2;
	query_parameter(argc,argv,q++, "Enter wo (MHz): ", wo);
	wo*=1.e6*PI2;
	query_parameter(argc,argv,q++, "Enter J(Hz): ", J);
	J*=PI2;
	double etap;
	query_parameter(argc,argv,q++, "Enter eta(0..1): ", etap);
	double rotp;
	query_parameter(argc,argv,q++, "Inital Rotor Phase(deg): ", rotp);
	rotp*=Pi/180.0;
	double tm;
	query_parameter(argc,argv,q++, "Time Step...(s): ", tm);
	std::string fname;
	query_parameter(argc,argv,q++, "File Output: ", fname);

	double Co=1./8.*sqrt(3.)/sqrt(2.)*(3.*cos(thetap)*cos(thetap)-1.+etap*sin(thetap)*sin(thetap)*cos(2.*phip));
	double C1=1./4.*sqrt(3./2.)*sin(2.*thetap)*(1.-1./3.*etap*cos(2.*phip));
	double C2=1./8.*sqrt(3./2.)*(sin(thetap)*sin(thetap)+etap/3.*(cos(thetap)*cos(thetap)+1.)*cos(2.*phip));
	double S1=1./2./sqrt(6.)*etap*sin(thetap)*sin(2.*phip);
	double S2=1./4./sqrt(6.)*etap*cos(thetap)*sin(2.*phip);

	matrix matro(5,5);

	matro(0,0)=0;
	matro(0,1)=-sqrt(3./2.)*C1*2*sqrt(2.)/3.;
	matro(0,2)=-sqrt(3./2.)*S1*2*sqrt(2.)/3.;
	matro(0,3)=sqrt(6.)*C2*2./3.;
	matro(0,4)=-sqrt(6.)*S2*2./3.;
	matro(1,0)=Co*2.*sqrt(2.)/3.;
	matro(1,1)=C1*(-1./3.)+S1*complexi/sqrt(3.);
	matro(1,2)=-C1*complexi*1/sqrt(3.)+S1*(-1./3.);
	matro(1,3)=-C2*(2.*sqrt(2.)/3.)+S2*2*complexi*sqrt(2./3.);
	matro(1,4)=C2*2.*complexi*sqrt(2./3.)+S2*(2.*sqrt(2.)/3.);
	matro(2,0)=-Co*2.*sqrt(2.)/3.;
	matro(2,1)=-C1*(-1./3.)+complexi*1/sqrt(3.)*S1;
	matro(2,2)=-C1*complexi*1/sqrt(3.)-S1*(-1./3.);
	matro(2,3)=C2*(2.*sqrt(2.)/3.)+2*S2*complexi*sqrt(2./3.);
	matro(2,4)=C2*2.*complexi*sqrt(2./3.)-S2*(2.*sqrt(2.)/3.);
	matro(3,0)=Co*2./3.;
	matro(3,1)=C1*(1./2.*2.*sqrt(2.)/3.)+S1*complexi*sqrt(2./3.);
	matro(3,2)=-C1*complexi*sqrt(2./3.)+S1*(sqrt(2.)/3.);
	matro(3,3)=C2*((1.+1./3.))-S2*2.*complexi*1./sqrt(3.);
	matro(3,4)=-C2*2.*complexi*1./sqrt(3.)-S2*((1.+1./3.));
	matro(4,0)=Co*2./3.;
	matro(4,1)=C1*(1./2.*2.*sqrt(2.)/3.)-complexi*sqrt(2./3.)*S1;
	matro(4,2)=C1*complexi*sqrt(2./3.)+S1*(sqrt(2.)/3.);
	matro(4,3)=C2*((1.+1./3.))+S2*complexi*1./sqrt(3.);
	matro(4,4)=C2*2.*complexi*1./sqrt(3.)-S2*((1.+1./3.));

	//cout<<mysys<<endl;
	//cout<<matro<<endl;

	matrix detect=mysys.Ip(0);
	hmatrix roeq=mysys.Ix(0);

	Vector<complex>  A(5,0);
	Vector<double> xo(5,0);


	Vector<complex> fid(npts, 0);
	xo(0)=tm;
	hmatrix H=mysys.Fe();
	matrix U=mysys.Fe();

	matrix T22=T2(mysys.A,1,1,2),
	T2m2=T2(mysys.A,1,1,-2),
	T21=T2(mysys.A,1,1,1),
	T2m1=T2(mysys.A,1,1,-1),
	T20=T2(mysys.A,1,1,0);


	matrix Iz1Iz2=mysys.Iz(0)*mysys.Iz(1);

#ifdef HAVE_MPI
	if(Nwr%size !=0)
	{
		std::cerr<<"ERROR: ***"<<endl;
		std::cerr<<" Number of rotor steps must be a "<<endl
		<<" multiple of the number of processors..."<<endl
		<<" THe end.."<<endl;
		exit(0);
	}

	int start=rank*int(Nwr/size);
	int end=(rank+1)*int(Nwr/size);
	string fbase=fname;
	double wrb=wr;
	for(int i=start;i<end;++i)
	{
		fname=fbase+itost(i);
		wr=wrb+double(i)*wrs;
		cout<<"PROC:: "<<i<<" "<<fname<<" "<<wr<<endl;

#endif
		mCos x1(wr, rotp);
		mSin x2(wr, rotp);
		mCos x3(2.0*wr, 2.0*rotp);
		mSin x4(2.0*wr, 2.0*rotp);

		Integrate<mCos> intC1(x1);
		Integrate<mSin> intS1(x2);
		Integrate<mCos> intC2(x3);
		Integrate<mSin> intS2(x4);



	for(double k=0;k<npts;++k)
	{
		cout<<"on Pt: "<<k<<endl;
		for(double kk=1;kk<=sample;++kk)
		{
			xo(1)=intC1((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
			xo(2)=intS1((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
			xo(3)=intC2((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);
			xo(4)=intS2((sample*(k)+(kk-1.0))*tm, (sample*k+kk)*tm);

			A=matro*xo;
			H=	wo*tm*mysys.Iz(1)+
				wq*sqrt(3.0/2.0)*
				(
					A(4)*T22
					+A(3)*T2m2
					-A(2)*T21
					+A(1)*T2m1
					+A(0)*sqrt(2./3.)*T20
				)
				+J*tm*Iz1Iz2;
			/*cout<<xo<<endl<<A*1.e9<<endl<<wq*sqrt(3.0/2.0)*
				(
					A(4)*T22
					+A(3)*T2m2
					-A(2)*T21
					+A(1)*T2m1
					+A(0)*sqrt(2./3.)*T20
				)<<endl;
			*/
			U=Mexp(H, -complexi)*U;
		}

		//printTime(k+1);
		fid(k)=trace(detect*U*roeq*adjoint(U));

	}
	plotterFID(fid, fname, tm*sample);
#ifdef HAVE_MPI
}
 MPI_Finalize();
#endif

}

