

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

//the shaped pulse reader and functor class
class ShapePulse{
	private:
		int curpos;
	public:
		double amplitude;
		double dt;
		Vector<double> phases;
		ShapePulse(){}
		ShapePulse(const Vector<string> &input)
		{
			Vector<string> tmp;
			for(int i=0;i<input.size();i++)
			{
				tmp=parse_param(input[i]);
				if(tmp.size()>0){
					if(tmp[0]=="amp"){	amplitude= atof(tmp[1].c_str());	}
					if(tmp[0]=="dt"){	dt= atof(tmp[1].c_str());	}

					if(isdigit(tmp[0][0])){	phases.push_back(atof(tmp[0].c_str())*PI/180.);	}
				}
			}
			curpos=0;
		}

		void operator++(){	curpos++;	}
		void operator++(int){	curpos++;	}
		void reset(){	curpos=0;	}
		int on(){	return curpos;	}
		operator bool()
		{
			return curpos<phases.size();
		}

		hmatrix Pulse(int on, SolidSys &A)
		{
			return amplitude*(A.Ix(on)*cos(phases(curpos)) +
					A.Iy(on)*sin(phases(curpos)));
		}

		hmatrix Pulse(SolidSys &A)
		{
			return amplitude*(A.Fx()*cos(phases(curpos)) +
					A.Fy()*sin(phases(curpos)));
		}
};

//standard pulse on a 'spin'
matrix Pulse(double theta, double phase, int on, SolidSys &A)
{
	return theta*PI/180.*(A.Ix(on)*cos(phase*PI/180.) + A.Iy(on)*sin(phase*PI/180.));
}

//standard pulse on ALL the spins
matrix Pulse(double theta, double phase, SolidSys &A)
{
	return theta*PI/180.*(A.Fx()*cos(phase*PI/180.) + A.Fy()*sin(phase*PI/180.));
}

int main(int argc,char* argv[]){

	std::string infile="";
	int q=1;

//get the input parameter file name
	query_parameter(argc,argv,q++, "input file: ", infile);

//the parameter section
	Parameters pset(infile);

//parsse up the parameter file according to each subsection
	pset.addSection("spins");
	pset.addSection("shape");
	pset.addSection("params");


//start grabbing the various items from the file
	int npts=pset.getParamI("npts", "params");
	double sw=pset.getParamD("sw","params");
	string outf=pset.getParamS("outfile", "params");
	double pamp=pset.getParamD("pangle", "params");
	double pphase=pset.getParamD("pphase", "params");
	int step=pset.getParamI("isostep", "params");
	double stsize=pset.getParamD("stepsize", "params");

//SolidSys uses can read the parsed 'spins' chunk
	SolidSys sys(pset.Section("spins"));

//print out to console so we can see all is well
	cout<<sys<<endl;

//The shaped pulse class reads in the sections
	ShapePulse sp(pset.Section("shape"));

//print out some info about the shape
	cout<<"Amp::"<<sp.amplitude<<endl;
	cout<<"dt::"<<sp.dt<<endl;

//the QM parts for the simulation
	hmatrix ro=sys.Fz(),H;
	matrix det=sys.Fmi();
	string fbase=outf;
	matrix U;
	double holdiso=sys.csa[1].iso();

//going to loop over several steps in the isotropic shift
	for(int i=0;i<step;++i){
		sys.csa[1].iso(holdiso+double(i)*stsize); //set the isotropic shift for the second spin
		ro=sys.Fz(); //reset the intial density matrix
		U=Mexp(Pulse(pamp, pphase, sys), -complexi); //the inital 'prep' pulse
		sp.reset(); //restart out shaped pulse iterator at the begining
		while(sp)
		{
			H=(sp.Pulse(sys)+sys.H(0,0,0,0)); //the total hamiltonian the shaped pulse piece and the interactions
			U=Mexp(H,-complexi*sp.dt*PI2)*U; //evolve
			++sp; //advance the shape pulse iterator
			cout<<sp.on()<<"\r"; //out put some info
			cout.flush();
		}

		U=Mexp(Pulse(pamp, pphase+180.0, sys),-complexi)*U; //the final 'unprep' pulse
		ro=prop(U,ro); //propogate ro

		oneFID<SolidSys> myfid(sys, npts, sw); //set up an FID object
		Vector<complex> fid=myfid.FID(ro, det); //get the FID

		string fna=fbase+itost(i); //advance the output file iterator
		plotterFID(fid, fna, 1.0/myfid.sw()); //print this fid data to a file
	}
}


