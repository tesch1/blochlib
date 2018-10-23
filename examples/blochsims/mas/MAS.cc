

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*
A simple 90 pulse on a PARTICLES of spins (rectangular grids)
*/

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout << "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}

void Info(std::string mess)
{
	cout<<mess;
	cout.flush();
}


//SolverOps::All the TypeDefs To make life simple
/*** GRIDS ***/
typedef XYZfull TheShape;
typedef XYZshape<TheShape> TheGrida;	//our shape
typedef RotatingGrid<TheGrida> TheGrid; //spinning grid

/*** parameters ***/
typedef ListBlochParams< TheGrid, BPoptions::Particle | BPoptions::HighField, double > MyPars;

/*** Interactions ***/
typedef TanhScale Scaler;
typedef Interactions<Offset<>, Relax<>, BulkSus, RadDamp, DimLessDipole<TheGrid> > MyInteractions;

/*** Solver Sets ***/
typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

int main(int argc,char* argv[]){

	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	Parameters pset(fn);

	int nsteps=pset.getParamI("npts");
	double tf=pset.getParamD("tf");
	double D=pset.getParamD("bulksus");
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

	std::string fout=pset.getParamS("fidout","",false,"fid");
	std::string magout=pset.getParamS("magout","",false,"mag");
	bool savesubs=pset.getParamI("saveSubs","",false,1);
	int showProg=pset.getParamI("showProgress","",false,1);
	std::string lypname=pset.getParamS("lypout","",false,"lyps");

	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("min"));
	coord<> maxs(pset.getParamCoordD("max"));
	coord<> myAxis(pset.getParamCoordD("rotaxis"));
	double rate=pset.getParamD("spinrate")*PI2;

	double pang=pset.getParamD("pulseangle");
	double amp=pset.getParamD("pulseamp");
	double phase=pset.getParamD("pulsephase")*PI/180.0;
	int cv=pset.getParamI("lyps");

	std::string dataou=pset.getParamS("trajectories", "", false);

/*** GRIDS ***/

	if(showProg>0) Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims); //the basic grid

	if(showProg>0) Info("Creating inital shape....\n");
	TheShape tester;
	if(showProg>0) Info("Creating total shape-grid....\n");
	TheGrida jja( gg, tester);
	TheGrid jj(jja, myAxis, rate);

//this prionts out a series of grids at times so you can see the
// grid rotation
	TheGrid::iterator myit(jj);
	string fbase="grid";
	int myct=0;
	double ttf=1.0/rate/2.0;
	for(double ct=0.0;ct<ttf && rate>0.0; ct+=ttf/10.0)
	{
		myit.reset();
		string fname=fbase+itost(myct); ++myct;
		ofstream oo(fname.c_str());
		while(myit){
			oo<<myit.Point(ct)<<endl;
			++myit;
		}
	}

/****LIST BLOCH PARAMETERS *****/
	//get initial condition
	std::string initc=pset.getParamS("initcond", "", false, "AllUp");
	int IC;
	int iters=1;

//parameters
	IC=InitialCondition::getFlag(initc);
	if(IC & InitialCondition::Random)
		iters=abs(pset.getParamI("numRand", "", false, 1));

	if(showProg!=0) Info("Creating entire spin parameter list for "+itost(jj.size())+" spins....\n");
	MyPars mypars(jj.size(), "1H", jj, IC);

	if(showProg!=0) Info("Setting spin parameter offsets....\n");
	for(int j=0;j<mypars.size();j++){
		mypars(j)=spintype;
		mypars(j).moles(moles);
		mypars(j).Bo(inBo);
		mypars.temperature(inTemp);
	}

	mypars.calcTotalMo();
	if(showProg>1) mypars.print(cout);

//The pulse list for a real pulse on protons..
	if(showProg>0) Info("Creating real pulse lists...\n");
	Pulse PP1(spintype, amp*PI2, phase); // (spin, amplitude, phase, offset)
	if(showProg>0)  PP1.print(cout);
	double tpulse=PP1.timeForAngle(pang*Pi/180., spintype);

//Extra ineractions
	if(showProg>0)  Info("Setting Interactions....\n");
	Offset<> myOffs(mypars, offset);
	Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	DimLessDipole<TheGrid> DipDip(jj, dipstr);
	if(abs(rate)>0.0) DipDip.Dynamic=true; //turn on the dyanimc flag
	if(dipstr==0)	DipDip.Off();
	MyInteractions MyInts(myOffs, myRels,BsRun, RdRun, DipDip);


//The poinciar section parameters
	//to do it or not
	std::string doPoinc=pset.getParamS("doPoinc", "",false, "no");
	//the 'phase' on the poiniar torus
	double poincPhase=pset.getParamD("poincPhase", "", false, 0)*PI/180.0;
	int poincPts=0;
	if(doPoinc=="yes"  && rate>0.0){ poincPts=pset.getParamI("poincPts");	}

//Set up the Time trains
	if(showProg>0)  Info("Initializing Time train for first Pulse....\n");
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

	//altered for poincar or not...
	if(showProg>0) Info("Initializing Time train for FID....\n");
	if(doPoinc=="yes" && rate>0.0){
		nsteps=poincPts;
		tf=(PI2/rate)*double(nsteps); //collect a point every wr cycle
		if(showProg>0)  Info("Final Time for Poinciar Sections...."+dbtost(tf)+"\n");

	}
	TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,20));

//total fid and mag sums
	Vector<complex>  totFid(nsteps, 0.0);
	Vector<coord<> > totMag(nsteps, 0.0);

	for(int RunCt=0;RunCt<iters;++RunCt)
	{
		mypars.calcInitialCondition();
		if(showProg>1) mypars.print(cout);

	//THis is the BLoch solve to perform a pulse
		if(showProg>0) Info("Initializing total parameter list with a pulse....\n");
		PulseBloch myparspulse(mypars, PP1, MyInts);
		if(cv) myparspulse.calcVariational();

	//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
		if(showProg>0) Info("Initializing total parameter list for FID collection....\n");
		NoPulseBloch me;
		me=(myparspulse);
		if(showProg>0) Info("Integrating first Pulse....\n");
		Vector<coord<> > tm=me.currentMag();

		stopwatch.reset();
		BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");
		drivP.setProgressBar(SolverOps::Off);
		if(showProg>1) drivP.setProgressBar(SolverOps::On);


	//output trajectory data if wanted
		std::ofstream trajout;
		std::string curTraj=dataou+itost(RunCt);
		if(dataou!="" && savesubs){
			if(!trajout.is_open()) trajout.open(curTraj.c_str());
			drivP.setCollectionPolicy(SolverOps::All);
		}else{
			drivP.setCollectionPolicy(SolverOps::FinalPoint);
		}

		drivP.setWritePolicy(SolverOps::Hold);
		if(pang!=0.0 && amp!=0.0){
			if(!drivP.solve(P1)){
				Info(" ERROR!!..could not integrate pulse P1....\n");
				return -1;
			}else if(dataou!=""){
				drivP.writeData(trajout);
			}
		}

		BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());
		driv.setProgressBar(SolverOps::Off);
		if(showProg>1) driv.setProgressBar(SolverOps::On);

		if(showProg>0) Info("Integrating FID ....\n");

		if(dataou!=""){		driv.setCollectionPolicy(SolverOps::All);	}
		else{		drivP.setCollectionPolicy(SolverOps::MagAndFID);		}

		driv.setWritePolicy(SolverOps::Hold);
		driv.setDetect(detsp);
		std::string lypCur=lypname+itost(RunCt);
		if(cv && savesubs) {
			driv.setLyapunovPolicy(SolverOps::LypContinous);
			driv.setLypDataFile(lypCur);
		}

		std::string foutCur=fout+itost(RunCt);
		std::string magCur=magout+itost(RunCt);

		if(driv.solve(F1)){
			totFid+=driv.FID();
			totMag+=driv.M();
			if(savesubs){
				driv.writeSpectrum(foutCur);
				driv.writeMag(magCur);
				if(dataou!="") driv.writeData(trajout);
			}
		}
		trajout.close();
		if(showProg>0){
			Info("On FID "+itost(RunCt)+" ");
			printTime();
		}
	}
	std::string foutCur=fout+"TOT";
	plotterFID(totFid, foutCur, tf/double(nsteps));
	foutCur=magout+"TOT";
	std::ofstream totMagOut(foutCur.c_str());
	totMag.print(totMagOut, "\n");

}



