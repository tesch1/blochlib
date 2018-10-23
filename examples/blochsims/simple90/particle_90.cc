

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;


/*
A simple 90 pulse on a PARTICLES of spins (rectangular grids)
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

	double pang=pset.getParamD("pulseangle");
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
	double tr=pset.getParamD("raddamp");
	double dipstr=pset.getParamD("dipole_str")*PI2;

	std::string fout=pset.getParamS("fidout");
	std::string magout=pset.getParamS("magout");

	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("min"));
	coord<> maxs(pset.getParamCoordD("max"));

	double amp=pset.getParamD("pulseamp");
	int cv=pset.getParamI("lyps");

	std::string dataou=pset.getParamS("trajectories", "", false);



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

//The pulse list for a real pulse on protons..
	Info("Creating real pulse lists...\n");

	Pulse PP1(spintype, amp*PI2, Pi/2.); // (spin, amplitude, phase, offset)


	Info("Setting spin parameter offsets....\n");
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles);
		mypars(j).Bo(inBo);
		mypars.temperature(inTemp);
	}

	mypars.calcTotalMo();
	mypars.print(cout);
	PP1.print(cout);
	double tpulse=PP1.timeForAngle(pang*Pi/180., spintype);

//time train testing

	Info("Initializing Time train for first Pulse....\n");
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

	Info("Initializing Time train for FID....\n");
	TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,20));

//Extra ineractions
	typedef TanhScale Scaler;
	typedef Interactions<Offset<>, Relax<>,  RadDamp, DimLessDipole<TheGrid> > MyInteractions;
	Info("Setting Interactions....\n");
	Offset<> myOffs(mypars, offset);
	Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	RadDamp RdRun(tr);
	DimLessDipole<TheGrid> DipDip(jj, dipstr);


	MyInteractions MyInts(myOffs, myRels, RdRun, DipDip);

//DipDip.Off();

//typedefs for Bloch parameter sets
	typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
	typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//THis is the BLoch solve to perform a pulse
	Info("Initializing total parameter list with a pulse....\n");
	PulseBloch myparspulse(mypars, PP1, MyInts);
	if(cv) myparspulse.calcVariational();

	if(dipstr==0)	DipDip.Off();

//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n");
	NoPulseBloch me;
	me=(myparspulse);
	Info("Integrating first Pulse....\n");
	Vector<coord<> > tm=me.currentMag();

	stopwatch.reset();
	BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");

//output trajectory data if wanted
	std::ofstream trajout;
	if(dataou!=""){
		trajout.open(dataou.c_str());
		drivP.setCollectionPolicy(SolverOps::All);
	}else{
		drivP.setCollectionPolicy(SolverOps::FinalPoint);
	}

	drivP.setWritePolicy(SolverOps::Hold);
	if(!drivP.solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n");
		return -1;
	}else if(dataou!=""){
		drivP.writeData(trajout);
	}

	BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());

	Info("\nIntegrating FID ....\n");

	std::string lypname="lyps";
	if(dataou!=""){
		driv.setCollectionPolicy(SolverOps::All);
	}else{
		drivP.setCollectionPolicy(SolverOps::MagAndFID);
	}
	driv.setWritePolicy(SolverOps::Hold);
	driv.setDetect(detsp);
	if(cv) {
		driv.setLyapunovPolicy(SolverOps::LypContinous);
		driv.setLypDataFile(lypname);
	}

	if(driv.solve(F1)){
		driv.writeSpectrum(fout);
		driv.writeMag(magout);
		if(dataou!="") driv.writeData(trajout);
	}


	printTime();

}



