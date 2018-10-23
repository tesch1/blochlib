

/* and intregration method tester */

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


class WillRos{
	public:
		rmatrix Jacobi;
		double S, R, B,C;

		WillRos(){
			Jacobi.resize(3,3);
			S=30.0; R=0.415; B=16.5; C=10.0;
		}

		void Jacobian(Vector<coord<> > &iny)
		{

			/* Will Ross Jacboi*/
			Jacobi(0,0)=S-R*iny[0].x()-iny[0].y()-iny[0].z();
			Jacobi(0,1)=-iny[0].x();
			Jacobi(0,2)=-iny[0].x();

			Jacobi(1,0)=iny[0].y();
			Jacobi(1,1)=iny[0].x()-C;
			Jacobi(1,2)=0;

			Jacobi(2,0)=-iny[0].z();
			Jacobi(2,1)=0;
			Jacobi(2,2)=B-iny[0].x()-iny[0].z();
		}

		void jacobian(double t, Vector<coord<> > &iny, rmatrix &dfdt){
			Jacobian(iny);
			dfdt=Jacobi;
		}

		void function(double t, Vector<coord<> > &iny, Vector<coord<> > &dydt){

		/* normal Will Ros*/

			dydt[0].x()=S*iny[0].x()-R*iny[0].x()*iny[0].x()-iny[0].x()*iny[0].y()-iny[0].x()*iny[0].z();
			dydt[0].y()=iny[0].x()*iny[0].y()-C*iny[0].y();
			dydt[0].z()=B*iny[0].z()-iny[0].x()*iny[0].z()-0.5*iny[0].z()*iny[0].z();

			if(iny.size()>1){
				Jacobian(iny);
				dydt.put(1, Jacobi*iny(1));
				dydt.put(2, Jacobi*iny(2));
				dydt.put(3, Jacobi*iny(3));
			}
		}
};


class Lorentz{
	public:
		rmatrix Jacobi;
		double S, R, B;

		Lorentz(){
			Jacobi.resize(3,3);
			S=10; R=28; B=8./3.;
		}

		void Jacobian(Vector<coord<> > &iny)
		{
		/* Lorentz jacobi  */
			Jacobi(0,0)=-S;
			Jacobi(0,1)=S;
			Jacobi(0,2)=0;

			Jacobi(1,0)=R-iny[0].z();
			Jacobi(1,1)=-1;
			Jacobi(1,2)=-iny[0].x();

			Jacobi(2,0)=iny[0].y();
			Jacobi(2,1)=iny[0].x();
			Jacobi(2,2)=-B;
		}

		void jacobian(double t, Vector<coord<> > &iny,  rmatrix &dfdy){
			Jacobian(iny);
			dfdy=Jacobi;
		}


		void function(double t, Vector<coord<> > &iny, Vector<coord<> > &dydt)
		{
			dydt[0].x()=S*(iny[0].y()-iny[0].x());
			dydt[0].y()=R*iny[0].x()-iny[0].y()-iny[0].x()*iny[0].z();
			dydt[0].z()=iny[0].x()*iny[0].y()-B*iny[0].z();

			if(iny.size()>1){
				Jacobian(iny);
				dydt.put(1, Jacobi*iny(1));
				dydt.put(2, Jacobi*iny(2));
				dydt.put(3, Jacobi*iny(3));
			}
		}
};

class LorentzD{
	public:
		rmatrix Jacobi;
		double S, R, B;

		LorentzD(){
			Jacobi.resize(3,3);
			S=10; R=28; B=8./3.;
		}

		void Jacobian(Vector<double > &iny)
		{
		/* Lorentz jacobi  */
			Jacobi(0,0)=-S;
			Jacobi(0,1)=S;
			Jacobi(0,2)=0;

			Jacobi(1,0)=R-iny[2];
			Jacobi(1,1)=-1;
			Jacobi(1,2)=-iny[0];

			Jacobi(2,0)=iny[1];
			Jacobi(2,1)=iny[0];
			Jacobi(2,2)=-B;
		}

		void jacobian(double t, Vector<double > &iny,  rmatrix &dfdy){
			dfdy.resize(3,3);
			/* Lorentz jacobi  */
			dfdy(0,0)=-S;
			dfdy(0,1)=S;
			dfdy(0,2)=0;

			dfdy(1,0)=R-iny[2];
			dfdy(1,1)=-1;
			dfdy(1,2)=-iny[0];

			dfdy(2,0)=iny[1];
			dfdy(2,1)=iny[0];
			dfdy(2,2)=-B;
		}


		void function(double t, Vector<double > &iny, Vector<double > &dydt)
		{
			dydt[0]=S*(iny[1]-iny[0]);
			dydt[1]=R*iny[0]-iny[1]-iny[0]*iny[2];
			dydt[2]=iny[0]*iny[1]-B*iny[2];

		/*	if(iny.size()>1){
				Jacobian(iny);
				dydt.put(1, Jacobi*iny(1));
				dydt.put(2, Jacobi*iny(2));
				dydt.put(3, Jacobi*iny(3));
			}*/
		}
};

//a 'stiff' set of equations
class VanDerPol{
	public:
		rmatrix Jacobi;
		double S, R, B;

		VanDerPol(){
			Jacobi.resize(2,2);
			S=0.013; R=1000.0; B=2500.0;
		}

		void Jacobian(Vector<double > &iny)
		{
		/* jacobi  */
			Jacobi(0,0)=0.0;
			Jacobi(0,1)=1.0;

			Jacobi(1,0)=R*2*iny[0]*iny[1]-1.0;
			Jacobi(1,1)=R*(1.0-iny[0]*iny[0]);
		}

		void jacobian(double t, Vector<double > &iny,  rmatrix &dfdy){
			dfdy.resize(2,2);
			dfdy(0,0)=0.0;
			dfdy(0,1)=1.0;

			dfdy(1,0)=R*2.0*iny[0]*iny[1]-1.0;
			dfdy(1,1)=R*(1.0-iny[0]*iny[0]);
		}


		void function(double t, Vector<double > &iny, Vector<double > &dydt)
		{
			dydt[0]=iny[1];
			dydt[1]=R*(1.0-iny[0]*iny[0])*iny[1]-iny[0];
		}
};


int main(int argc,char* argv[])
{
	typedef WillRos MyODE ;
	typedef coord<> T;
	MyODE MyDiffs;
	Vector<T > IC(1, 1),start(1, 0);

	//set the variation bits to the identity matrix to begin with
	//start[1].x()=1;
	//start[2].y()=1;
	//start[3].z()=1;

	//time info
	double tstep=.01;
	double startT=0.;
	double endT=100.;
	query_parameter(argc,argv,1, "end time: ", endT);
	int nsteps=int((endT-startT)/tstep)+1;
	query_parameter(argc,argv,2, "dt: ", tstep);
	double subst=tstep/10;

	//an integration object
	bs<MyODE, T > odes(MyDiffs);

	odes.relToler=1e-6; //set the relative precision to be 'higher' then normal (1e-6)
	odes.absToler=1e-3; //set the absolute precision to be 'higher' then normal (1e-3)
	//output data...
	string fname="data";
	ofstream oo(fname.c_str());
	ofstream lyo("lyp");
	Vector<T > *out=odes.get_out();
	Range R(0,1);
	cout<<start(R)<<endl;

	/** Initial conditions for 'Vector<coord<> >' types (LORETNZ and WILLROSS)  **/
	start[0].x()=10.;
	start[0].y()=5.;
	start[0].z()=6.;

	/** Initial conditions for 'Vector<double>' types (VANDERPOL)  **/
	//start[0]=2;
	//start[1]=0.;

	//A Lyapunov object (there is only 1 coord<> we care about)
	//Lyapunov<coord<> > myLyps(1, out);

	double tmS=startT;
	IC=start;
	odes.setInitialCondition(IC);
	int pt=0;
	Vector<double> times=(Spread<double>(0.0, endT, tstep));
	Vector<Vector<T > > data=odes.solve(times);

	/** Uncomment for 'Vector<double>' types (VANDERPOL)  **/
//	for(int i=0;i<data.size();i++)
	//	oo<<times[i]<<" "<<data[i][0]<<" "<<data[i][1]<<endl;//<<" "<<data[i][2]<<endl;

	/** Uncomment for 'Vector<coord<> >' types  **/
	for(int i=0;i<data.size();i++)
		oo<<times[i]<<" "<<data[i]<<endl;

	odes.print(cout);
	/** Uncomment for 'Vector<coord<> >' solving Lypapunov equations also types  **/
	/*for(int k=0;k<times.size()-1;++k)
	{
		odes.solve(times[k], times[k+1]);
		oo<<times[k]<<" "<<(*out)<<endl;
		//myLyps.calcLyapunov(tmS, tstep);
		//lyo<<myLyps;
		//tmS+=tstep;
		//pt++;
	}
	odes.print(cout, "\n");*/
	oo.close();
	printTime();
}






