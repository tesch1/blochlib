

#include "blochlib.h"


/* Homonuclear The CRAZED pulse sequeque */

/*

RF  ---90--T1---90----T2-FID
Grad --------Gzt--2TGz------

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



typedef XYZcylinder TheShape;
typedef XYZshape<TheShape> TheGridS;
typedef GradientGrid<TheGridS > TheGrid;
typedef ListBlochParams< TheGrid > MyPars;

//Extra ineractions
typedef TanhScale Scaler;
typedef Interactions<Offset<TheGrid>, Relax<>, BulkSus, RadDamp, DemagField<TheGrid, Scaler> > MyInteractions;

//typedefs for Bloch parameter sets
typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;


int main(int argc,char* argv[]){
	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	Parameters pset(fn);

	double pang1=pset.getParamD("pulseangle1");
	double pang2=pset.getParamD("pulseangle2");
	double amp=pset.getParamD("pulseamp");
	double delay=pset.getParamD("delay");

	int nsteps=pset.getParamI("npts");
	double tf=pset.getParamD("tf");

	std::string fout=pset.GetParamS("fidout");
	std::string magout=pset.GetParamS("magout");
	int contfid=pset.getParamI("allfid","",false);

	int cv=pset.getParamI("lyps", "", false);
	std::string lypfile=pset.GetParamS("lypout", "", false, "lyps");

	std::string dataou=pset.GetParamS("trajectories", "", false);

//gradient pars
	double gradtime1=pset.getParamD("gradtime1");	//first grad pulse time
	double gradtime2=pset.getParamD("gradtime2"); //second grad pulse time (typically 2*gradetime1)

/******************/
//Grids
	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("gmin"));
	coord<> maxs(pset.getParamCoordD("gmax"));

	coord<> smins(pset.getParamCoordD("smin"));
	coord<> smaxs(pset.getParamCoordD("smax"));


	Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims);
	Info("Creating inital shape....\n");
	TheShape tester(smins, smaxs);
	Info("Creating total shape-grid....\n");
	TheGridS grids( gg, tester);
	ofstream goo("grid");
	goo<<grids<<endl;
//create the gradient grids..
	double zgrad=pset.getParamD("zgrad");

	Info("Creating Gradient map grids....\n");
	TheGrid jj(grids);

	jj.G(0.0,0.0,zgrad);
/******************/

/******************/
//set up Parameter lists
	int nsp=jj.size();
	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n");
	MyPars mypars(jj.size(), "1H", jj);
	nsp=mypars.size();

	double inBo=pset.getParamD("Bo");
	double inTemp=pset.getParamD("temperature");
	string spintype=pset.GetParamS("spintype");
	double moles=pset.getParamD("moles");
	string detsp=spintype;


	Info("setting spin parameter offsets....\n");
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles);
		mypars(j).Bo(inBo);
		mypars(j).temperature(inTemp);
	}

	mypars.calcTotalMo();
	mypars.print(cout);
/****************/
//The pulse list for a real pulse on protons..
	Info("Creating real pulse lists...\n");

	Pulse PP1(spintype, amp, 0.); // (spin, amplitude, phase, offset)
	Pulse PP2(spintype, amp,0.); // (spin, amplitude, phase, offset)

	PP1.print(cout);
	PP2.print(cout);
	double tpulse=PP1.timeForAngle(pang1*Pi/180., spintype);
	double tpulse2=PP2.timeForAngle(pang2*Pi/180., spintype);

/********************************/
//time train testing
	double tct=0;
	Info("Initializing Time train for first Pulse....\n");
	TimeTrain<UniformTimeEngine > P1(0., tpulse, 10,100);

	tct+=tpulse;

	Info("Initializing Time train for Delay....\n");
	TimeTrain<UniformTimeEngine > D1(tct, tct+delay, 128,100);
	tct+=delay;
	Info("Initializing Time train for First Gradient Pulse....\n");
	TimeTrain<UniformTimeEngine > G1(tct, tct+gradtime1, 10,100);
	tct+=gradtime1;
	Info("Initializing Time train for second Pulse....\n");
	TimeTrain<UniformTimeEngine > P2(tct, tct+tpulse2, 10,100);
	tct+=tpulse2;
	Info("Initializing Time train for Second Gradient Pulse....\n");
	TimeTrain<UniformTimeEngine > G2(tct, tct+gradtime2, 10,100);
	tct+=gradtime2;
	Info("Initializing Time train for FID....\n");
	TimeTrain<UniformTimeEngine > F1(tct, tf+tct, nsteps,100);

/*******************************/
/******************************/
//interactions
	double D=pset.getParamD("bulksus");
	double t2s=pset.getParamD("T2");
	double t1s=pset.getParamD("T1");
	double offset=pset.getParamD("offset")*PI2;

	double tr=pset.getParamD("raddamp");
	double scalef=pset.getParamD("scalef");

	Info("setting Interactions....\n");

	Offset<TheGrid> myOffs(jj, offset);
	Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	Scaler myscal(scalef);
	DemagField<TheGrid, Scaler> DipDip(jj, myscal);
	MyInteractions MyInts(myOffs, myRels,BsRun, RdRun, DipDip);

	if(scalef<=0) DipDip.Off();
/********************************/

//THis is the BLoch solve to perform a pulse

	Info("Initializing total parameter list with a pulse....\n");
	PulseBloch myparspulse(mypars, PP1, MyInts);
	if(cv) myparspulse.calcVariational();

//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n");
	NoPulseBloch me;
	me=myparspulse;

	Vector<coord<> > tm=me.currentMag();

	stopwatch.reset();

	BlochSolver<PulseBloch > drivP(myparspulse, tm);

//output trajectory data if wanted
	std::ofstream trajout;
	if(dataou!=""){
		trajout.open(dataou.c_str());
		drivP.setCollectionPolicy(SolverOps::All);
	}else{
		drivP.setCollectionPolicy(SolverOps::FinalPoint);
	}

	if(cv) drivP.setLyapunovPolicy(SolverOps::LypContinous);
	if(cv) drivP.setLypDataFile(lypfile);


//integrate the first pulse
	myOffs.off();	//turn off gradient

	drivP.setWritePolicy(SolverOps::Hold);

	Info("Integrating first Pulse ....\n");
	if(!drivP.solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n");
		return -1;
	}else{
		if(dataou!="") drivP.writeData(trajout);
	}

//setup and integrate the first delay up to the first gradient pulse
	Info("\nIntegrating for First Delay....\n");
	BlochSolver<NoPulseBloch > drivD(me, drivP.lastPoint());
	if(dataou!=""){
		drivD.setCollectionPolicy(SolverOps::All);
	}else{
		drivD.setCollectionPolicy(SolverOps::FinalPoint);
	}

	if(cv) drivD.setLyapunovPolicy(SolverOps::LypContinous);
	if(cv) drivD.setLypDataFile(lypfile);

	if(!drivD.solve(D1)){
		Info(" ERROR!!..could not integrate Delay....\n");
		return -1;
	}else{
		if(dataou!="")	drivD.writeData(trajout);
	}

//integrate the gradient pulse
	Info("\nIntegrating the First Grad Pulse....\n");
	if(gradtime1>0){
		myOffs.on();	//turn on gradient
		if(!drivD.solve(G1)){
			Info(" ERROR!!..could not integrate G1....\n");
			return -1;
		}else{
			if(dataou!="")	drivD.writeData(trajout);
		}
	}

//iniegrate the second pulse
	myOffs.off();	//turn off grad
	Info("\nIntegrating for Second Pulse....\n");
	myparspulse.setPulses(PP2);
	drivP.setInitialCondition(drivD.lastPoint());
	if(!drivP.solve(P2)){
		Info(" ERROR!!..could not pulse P2....\n");
		return -1;
	}else{
		if(dataou!="") drivP.writeData(trajout);
	}

//integrate the second gradient
 	myOffs.on();
 	std::ofstream fiout(fout.c_str());
 	std::ofstream maout(magout.c_str());

	drivD.setInitialCondition(drivP.lastPoint());
	if(gradtime2>0){
		Info("\nIntegrating the second Grad Pulse ....\n");

		if(dataou!=""){
			drivD.setCollectionPolicy(SolverOps::All);
		}else{
			drivD.setCollectionPolicy(SolverOps::MagAndFID);
		}

		if(drivD.solve(G2)){
			if(contfid){
				drivD.writeSpectrum(fiout);
				drivD.writeMag(maout);
			}
			if(dataou!="") drivD.writeData(trajout);
		}
	}

//integrate FID
	myOffs.off();
	Info("\nIntegrating for FID ....\n");
	if(dataou!=""){
		drivD.setCollectionPolicy(SolverOps::All);
	}else{
		drivD.setCollectionPolicy(SolverOps::MagAndFID);
	}

	if(drivD.solve(F1)){
		drivD.writeSpectrum(fiout);
		drivD.writeMag(maout);
		if(dataou!="") drivD.writeData(trajout);
	}


	printTime();

	cout<<"\a"<<endl;



}


