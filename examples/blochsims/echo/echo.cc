

#include "blochlib.h"


/* A bloch simulation of the Hahn echo.... */


//the required 2 namespaces
using namespace BlochLib;
using namespace std;

timer stopwatch;
void printTime(ofstream &oo){
	oo <<std::endl<< "Time taken: " << (stopwatch()) << " seconds\n";
}

void Info(std::string mess, ofstream &out)
{
	out<<mess;
	out.flush();
}

int main(int argc,char* argv[])
{
	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	Parameters pset(fn);

	std::string fout=pset.getParamS("fidout", "", false, "fidout");
	std::string lypname=pset.getParamS("lypname", "", false, "lyps");
	std::string logout=pset.getParamS("logout", "", false, "echo.log");
	std::string magout=pset.getParamS("magout", "", false, "magout");
	int cv=pset.getParamI("lyps");
	std::string dataou=pset.getParamS("trajectories", "", false);

	std::ofstream logf(logout.c_str());


/***************** grid setup **************/
	typedef XYZfull TheShape;
	typedef XYZshape<TheShape> TheGridG;
	typedef GradientGrid<TheGridG> TheGrid;

	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("min"));
	coord<> maxs(pset.getParamCoordD("max"));

	Info("Creating grid....\n", logf);
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n", logf);
	TheShape tester;
	Info("Creating total shape-grid....\n", logf);
	TheGridG qq( gg, tester);
	TheGrid jj(qq);


/***************** gradient setup **************/
	jj.G()=pset.getParamCoordD("grad1");

/***************** parameter setup **************/
	typedef ListBlochParams< TheGrid, BPoptions::Density | BPoptions::HighField, double > MyPars;

	double moles=pset.getParamD("moles");
	double inBo=pset.getParamD("Bo");
	double inTemp=pset.getParamD("temperature");
	std::string spintype=pset.getParamS("spintype");
	std::string detsp=spintype;

//get initial condition
	std::string initc=pset.getParamS("initcond", "", false, "AllUp");
	int IC;

	if(initc=="RandomUpDown"){
		IC=InitialCondition::RandomUpDown;
	}else if(initc=="Random" || initc=="RandomDistribution"){
		IC=InitialCondition::RandomDistribution;
	}else if(initc=="AllDown"){
		IC=InitialCondition::AllDown;
	}else{
		IC=InitialCondition::AllUp;
	}

	int nsp=jj.size();
	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n", logf);
	MyPars mypars(nsp, "1H", jj, IC);
	nsp=mypars.size();
	Info("Setting spin parameters....\n", logf);
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles/double(nsp));
		mypars(j).Bo(inBo);
		mypars(j).temperature(inTemp);
	}

	mypars.calcTotalMo();
	mypars.print(logf);


/***************** Pulse setup **************/
//The pulse list for a real pulse on protons..
	Info("Creating real pulse lists...\n", logf);

	coord<> pang1=pset.getParamCoordD("pulse1");
	coord<> pang2=pset.getParamCoordD("pulse2");

	Pulse PP1(spintype, pang1[2]*PI2, pang1[1]*Pi/180.0); // (spin, amplitude, phase, offset)
	Pulse PP2(spintype, pang2[2]*PI2, pang2[1]*Pi/180.0); // (spin, amplitude, phase, offset)


	PP1.print(cout);
	PP2.print(cout);
	double tpulse1=PP1.timeForAngle(pang1[0]*Pi/180., spintype);
	double tpulse2=PP2.timeForAngle(pang2[0]*Pi/180., spintype);


/***************** Time Train setup **************/

	double delay1=pset.getParamD("delay1");
	double tf=pset.getParamD("tf");
	int nsteps=pset.getParamI("npts");

	double tt=0.0;
	Info("Initializing Time train for first Pulse....\n", logf);
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(tt, tt+tpulse1, 10,100));

	tt+=tpulse1;
	Info("Initializing Time train for first Delay....\n", logf);
	TimeTrain<UniformTimeEngine > D1(UniformTimeEngine(tt, tt+delay1, 10,100));

	tt+=delay1;
	Info("Initializing Time train for first Pulse....\n", logf);
	TimeTrain<UniformTimeEngine > P2(UniformTimeEngine(tt, tt+tpulse2, 10,100));

	tt+=tpulse2;
	Info("Initializing Time train for FID....\n", logf);
	TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tt, tt+tf, nsteps,20));

/***************** Interactions setup **************/
	typedef TanhScale Scaler;
	typedef Interactions<Offset<MyPars>, Relax<>,
						BulkSus, RadDamp,
						DipoleDipole<TheGrid, Scaler>,
						DemagField<TheGrid, Scaler> > MyInteractions;

	double t2s=pset.getParamD("T2");
	double t1s=pset.getParamD("T1");
	double tr=pset.getParamD("raddamp");
	double scalef=pset.getParamD("scalef");
	double D=pset.getParamD("bulksus");
	double offset=pset.getParamD("offset")*PI2;
	double demag=pset.getParamD("demag", "", false, 0);


	Info("Setting Interactions....\n", logf);
	Offset<MyPars> myOffs(mypars, offset);
	Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	Scaler scal(scalef);
	DipoleDipole<TheGrid, Scaler> DipDip(jj, scal);
	DemagField<TheGrid, Scaler> DipDipDM(jj, scal);

	if(scalef==0) DipDip.off();
	if(demag==0 || scalef==0) DipDipDM.off();

	MyInteractions MyInts(myOffs, myRels,BsRun, RdRun, DipDip, DipDipDM);

//typedefs for Bloch parameter sets
	typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
	typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//THis is the BLoch solve to perform a pulse
	Info("Initializing total parameter list with a pulse....\n", logf);
	PulseBloch myparspulse(mypars, PP1, MyInts);
	if(cv) myparspulse.calcVariational();


//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n", logf);
	NoPulseBloch me;
	me=(myparspulse);
	Vector<coord<> > tm=me.currentMag();
	Info("\nInitial Condition (Mo(0))...\n", logf);
	for(int i=0;i<tm.size();++i)	logf<<"[ "<<tm[i]<<"]"<<std::endl;

	Info("Integrating first Pulse....\n", logf);

	stopwatch.reset();

//turn off grad if there is none
	if(jj.Gx()==0.0 && jj.Gy()==0.0 && jj.Gz()==0.0)
	{	myOffs.off();	}
	Info("\nGradient Strength (x,y,z) Gauss/cm ... ", logf);
	logf<<"[ "<<jj.G()<<"]"<<std::endl;

/************** First Pulse *************/
	BlochSolver<PulseBloch > drivP(myparspulse, tm);

//output trajectory data if wanted
	if(dataou!=""){
		drivP.setWritePolicy(SolverOps::Continous);
		drivP.setRawOut(dataou);
	}else{
		drivP.setWritePolicy(SolverOps::Hold);
	}
	drivP.setCollectionPolicy(SolverOps::FinalPoint);

	if(!drivP.solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n", logf);
		return -1;
	}

/************** First Delay *************/
	BlochSolver<NoPulseBloch > drivD(me, drivP.lastPoint());

	Info("\nIntegrating First Delay ....\n", logf);
	if(dataou!=""){
		drivD.setWritePolicy(SolverOps::Continous);
		drivD.setRawOut(dataou, ios::app|ios::out);
	}else{
		drivD.setWritePolicy(SolverOps::Hold);
	}

	drivD.setCollectionPolicy(SolverOps::FinalPoint);
	if(!drivD.solve(D1)){
		Info(" ERROR!!..could not integrate delay D1....\n", logf);
		return -1;
	}

/************** Second Pulse *************/
	if(dataou!=""){
		drivP.setWritePolicy(SolverOps::Continous);
		drivP.setRawOut(dataou, ios::app|ios::out);
	}

	drivP.setInitialCondition(drivD.lastPoint());
	Info("\nIntegrating Second Pulse ....\n", logf);
	myparspulse.setPulses(PP2);
	if(!drivP.solve(P2)){
		Info(" ERROR!!..could not integrate pulse P2....\n", logf);
		return -1;
	}

/************** FID *************/
	jj.G()=pset.getParamCoordD("grad2");
	Info("\nIntegrating FID ....\n", logf);
	drivD.setCollectionPolicy(SolverOps::MagAndFID);
	drivD.setDetect(detsp);
	drivD.setInitialCondition(drivP.lastPoint());

	if(cv) {
		drivD.setLyapunovPolicy(SolverOps::LypContinous);
		drivD.setLypDataFile(lypname);
	}

//output trajectory data if wanted

	if(drivD.solve(F1)){
		drivD.writeSpectrum(fout);
		drivD.writeMag(magout);
	}
	printTime(logf);
	std::cout<<"\a"<<std::endl;
	return 1;
}



