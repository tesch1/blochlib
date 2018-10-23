

#include "blochlib.h"

/* This is a  simple EPI (echo-planare-imaging)
	experiment...
	this exp is designed to calculate an entire
	directional enocding in a 1D experiment
	(not the typical '2D' one)...

	the experiement is
	(pulse)-----(FID)------(FID)...
	-------(-Gx)( Gx)(-2Gx)( Gx)...
	-------0-------------(Gy)---...
	-------( t) ( 2t)( t  )( 2t)...


*/
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

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
	int q=1;
	std::string psetfname, fname, matfname;
	query_parameter(argc,argv,q++, "Parameter file?: ", psetfname);
	Parameters pset(psetfname);
	pset.addSection("grid");

	int asORmat=pset.getParamI("asciiORmat", "", false, 0);
	if(asORmat==0){
		fname=pset.getParamS("fbase", "", false, "oepi");
	}else{
		matfname=pset.getParamS("matname", "", false, "oepi.mat");
	}

	double ang=pset.getParamD("pangle","");
	double amp=pset.getParamD("pamp", "")*PI2;
	int nsteps=pset.getParamI("npts","");
	int N=pset.getParamI("npts2d","");

	double delay=pset.getParamD("delay","");

	string spintype=pset.getParamS("spintype", "", false, "1H");
	double D=pset.getParamD("bulksus","");
	double tr=pset.getParamD("raddamp", "");;
	double scalef=pset.getParamD("scalef","");
	double offset=pset.getParamD("offset","")*PI2;
	double t2s=pset.getParamD("T2", "");
	double t1s=pset.getParamD("T1", "");

	string detsp=spintype;

	int nsp;

// Bloch set up


	typedef XYZfull TheShape;
	typedef XYZshape<TheShape> TheGrid;

	coord<int> dims(pset.getParamCoordI("dim", "grid"));
	coord<> mins(pset.getParamCoordD("min", "grid"));
	coord<> maxs(pset.getParamCoordD("max", "grid"));


	Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n");
//	TheShape tester(0,.01,0,PI2, -.01,.01);
	TheShape tester;
	Info("Creating total shape-grid....\n");
	TheGrid jj(gg,tester);

	double xgrad=pset.getParamD("xgrad","");
	double ygrad=pset.getParamD("ygrad","");
	double gradtime=pset.getParamD("gradtime","");

	std::cout<<" Gradient Pulse will be applied for "<<gradtime<<" s"<<endl;
	std::cout<<" Gx="<<xgrad<<" G/cm Gy="<<ygrad<<" G/cm"<<endl;

	Info("Creating Gradient map grids....\n");
	GradientGrid<TheGrid > mygrid(jj);
	mygrid.G(xgrad,0.,0.);
	nsp=mygrid.size();

	typedef ListBlochParams< GradientGrid<TheGrid >,BPoptions::HighField | BPoptions::Particle, double> MyPars;

	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n");
	MyPars mypars(nsp, spintype, mygrid);

//The pulse list for a real pulse on protons..
	Info("Creating real pulse lists...\n");

	Pulse PP1(spintype, amp, Pi/2.); // (spin, amplitude, phase, offset)


	Info("setting spin parameter offsets....\n");

	mypars.calcTotalMo();
	mypars.print(cout);
	PP1.print(cout);

	double tpulse=PP1.timeForAngle(ang*Pi/180., spintype);

	ang*=Pi/180.;

//time train

	Info("Initializing Time train for first Pulse....\n");
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

	double tct=tpulse;
	Info("Initializing Time train for first Gradient Pulse....\n");
	TimeTrain<UniformTimeEngine > G1(UniformTimeEngine(tct, tct+delay, 10,100));
	tct+=delay;



//Extra ineractions
	typedef TanhScale Scaler;
	typedef Interactions<Offset<MyPars>, Relax<>, BulkSus, RadDamp, DipoleDipole<TheGrid, Scaler> > MyInteractions;
	//typedef Interactions<BulkSus, RadDamp> MyInteractions;
	Info("setting Interactions....\n");
	Offset<MyPars> offs(mypars,offset);
	Relax<> rels(mypars, t2s, t1s);
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	Scaler myscal(scalef);
	DipoleDipole<TheGrid, Scaler> DipDip(mygrid, myscal);
	MyInteractions MyInts(offs, rels,BsRun, RdRun, DipDip);

	if(scalef<=0) DipDip.Off();

//typedefs for Bloch parameter sets
	typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
	typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//THis is the BLoch solve to perform a pulse
	Info("Initializing total parameter list with a pulse....\n");
	PulseBloch myparspulse(mypars, PP1, MyInts);
//	myparspulse.calcVariational();

//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n");
	NoPulseBloch me;
	me=(myparspulse);
	Info("Integrating first Pulse (Gradient Off)....\n");
	Vector<coord<> > tm=me.currentMag();

	stopwatch.reset();
	mypars.GradOff();
	BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");
	drivP.setWritePolicy(SolverOps::Hold);
	drivP.setCollectionPolicy(SolverOps::FinalPoint);
	if(!drivP.solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n");
		return -1;
	}

	Info("\nIntegrating First Gradient Pulse....\n");
	mypars.GradOn();
	mygrid.Gx(-abs(mygrid.Gx()));
	BlochSolver<NoPulseBloch> drivD(me, drivP.lastPoint(), "out");
	drivD.setWritePolicy(SolverOps::Hold);
	drivD.setCollectionPolicy(SolverOps::FinalPoint);
	if(!drivD.solve(G1)){
		Info(" ERROR!!..could not Gradient pulse P2....\n");
		return -1;
	}

	matrix FID(N, nsteps);

	for(int i=0;i<N;i++){

		ofstream fidd;

		Info("Initializing Time train for FID....\n");
		TimeTrain<UniformTimeEngine > F1(tct, tct+2.0*delay, nsteps,10);
		tct+=2.0*delay;

		if(asORmat==0){
			std::string sfid=fname+itost(i);
			ofstream fidd(sfid.c_str());
		}

		mygrid.Gx(abs(xgrad));
		drivD.setCollectionPolicy(SolverOps::FIDonly);
		drivD.solve(F1);
		if(asORmat==0){
			drivD.writeSpectrum(fidd);
		}else{
			FID.putCol(i, drivD.FID());
		}

		mygrid.Gy(ygrad);
		mygrid.Gx(-2.0*abs(xgrad));
		drivD.setCollectionPolicy(SolverOps::FinalPoint);

	//little y grad  bit
		Info("Initializing Time train Gx(refocus)--Gy gradient Pulse....\n");
		TimeTrain<UniformTimeEngine > G2a(tct, tct+gradtime, 50,10);
		tct+=gradtime;
		drivD.solve(G2a);

	//refocus x grad
		Info("Initializing Time train Gx(refocus) gradient Pulse....\n");
		TimeTrain<UniformTimeEngine > G2b(tct, tct+delay-gradtime, 50,10);
		tct+=delay-gradtime;
		mygrid.Gy(0.0);
		drivD.solve(G2b);
	}
	if(asORmat!=0){
		matstream matout(matfname, ios::out | ios::binary);
		matout.put("epi", FID);
	}
	printTime();
}



