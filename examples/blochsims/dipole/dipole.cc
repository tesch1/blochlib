

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*

loops through various pulse angles
on a dipole-dipole coupled system using the DimLessDipole object
offsets, and relaxation parmeters

*/

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}

void Info(std::string mess)
{
	cout<<mess;
	cout.flush();
}

int main(int argc,char* argv[]){




	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	Parameters pset(fn);

	int nsteps=pset.getParamI("npts");
	double tf=pset.getParamD("tf");
	double offset=pset.getParamD("offset")*PI2;
	double inBo=pset.getParamD("Bo");
	double inTemp=pset.getParamD("temperature");
	string spintype=pset.getParamS("spintype");
	string detsp=spintype;
	double t2s=pset.getParamD("T2");
	double t1s=pset.getParamD("T1");
	double moles=pset.getParamD("moles");
	double dipstr=pset.getParamD("dipole_str")*PI2;

	std::string fout=pset.getParamS("fidout");

	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("min"));
	coord<> maxs(pset.getParamCoordD("max"));

	int cv=pset.getParamI("lyps");


// Bloch set up testing


	typedef XYZfull TheShape;
	typedef XYZshape<TheShape> TheGrid;

	Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n");
	TheShape tester;
	Info("Creating total shape-grid....\n");
	TheGrid jj( gg, tester);

	typedef ListBlochParams< TheGrid, BPoptions::Particle | BPoptions::HighField, double > MyPars;
	int nsp=jj.size();
	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n");
	MyPars mypars(nsp, "1H", jj);
	nsp=mypars.size();

	Info("Setting spin parameter offsets....\n");
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles);
		mypars(j).Bo(inBo);
		mypars.temperature(inTemp);
	}

	mypars.calcTotalMo();
	mypars.print(cout);

//time train testing


//Extra ineractions
	typedef TanhScale Scaler;
	typedef Interactions<Offset<>, Relax<>,  DimLessDipole<TheGrid> > MyInteractions;
	Info("Setting Interactions....\n");
	Offset<> myOffs(mypars, offset);
	Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	DimLessDipole<TheGrid> DipDip(jj, dipstr);


	MyInteractions MyInts(myOffs, myRels,  DipDip);


//the pulse object
//The pulse list for a real pulse on protons..
	double pang=pset.getParamD("pulseangle");
	double pstep=pset.getParamD("pulsestepsize");
	int numP=pset.getParamI("pulsesteps");
	double amp=pset.getParamD("pulseamp");
	double phase=pset.getParamD("pulsephase");

	Info("Creating real pulse lists...\n");

	Pulse PP1(spintype, amp*PI2, phase*DEG2RAD); // (spin, amplitude, phase, offset)

	PP1.print(cout);

//data FID
	matrix FIDs(numP, nsteps);

	for(int kk=0;kk<numP;++kk)
	{

		double tpulse=PP1.timeForAngle((pang+double(kk)*pstep)*Pi/180., spintype);
		std::cout<<std::endl<<"On Pulse Angle: "<<(pang+double(kk)*pstep)<<" degrees "<<std::endl;
		Info("Initializing Time train for first Pulse....\n");
		TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

		Info("Initializing Time train for FID....\n");
		TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,20));

	//typedefs for Bloch parameter sets
		typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
		typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

	//THis is the BLoch solve to perform a pulse
		Info("Initializing total parameter list with a pulse....\n");
		PulseBloch myparspulse(mypars, PP1, MyInts);
		if(cv) myparspulse.calcVariational();

		if(dipstr==0) DipDip.Off();

	//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
		Info("Initializing total parameter list for FID collection....\n");
		NoPulseBloch me;
		me=(myparspulse);
		Info("Integrating first Pulse....\n");
		Vector<coord<> > tm=me.currentMag();

		stopwatch.reset();
		BlochSolver<PulseBloch > drivP(myparspulse, tm);

		drivP.setCollectionPolicy(SolverOps::FinalPoint);

		drivP.setWritePolicy(SolverOps::Hold);
		if(!drivP.solve(P1)){
			Info(" ERROR!!..could not integrate pulse P1....\n");
			return -1;
		}

		BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());

		Info("\nIntegrating FID ....\n");

		std::string lypname="lyps";
		drivP.setCollectionPolicy(SolverOps::MagAndFID);

		driv.setWritePolicy(SolverOps::Hold);
		driv.setDetect(detsp);
		if(cv) {
			driv.setLyapunovPolicy(SolverOps::LypContinous);
			driv.setLypDataFile(lypname);
		}

		if(driv.solve(F1)){
			FIDs.putRow(kk, driv.FID());
		}

	}
	matstream matout(fout);
	matout.put("vdat", FIDs);
	Vector<double> pangs(Spread<double>(pang, pang+(numP*pstep), pstep));
	matout.put("pangs", pangs);
	matout.close();
	printTime();

}



