

#include "blochlib.h"

#ifdef WITH_MPI
	#include "mpi.h"
#endif

void WriteGnu(){
	ofstream gn("plotter");
	std::string mess="set data style lines\n";
	mess+="set origin 0,0\n";
	mess+="set size 1,1\n";
	mess+="set title \"Bloch Eq Sim\"\n";
	mess+="set nokey\n";
mess+="\n";
	mess+="set multiplot\n";
	mess+="set border 31\n";
mess+="\n";
	mess+="set ylabel \"\"\n";
	mess+="set xlabel \"t\", 1\n";
mess+="\n";
	mess+="set origin 0,0\n";
	mess+="set size .3,.48\n";
	mess+="set title \"imag FID\" ,-1\n";
	mess+="plot \"fid\" u 1:6\n";
mess+="\n";
	mess+="set origin 0,.52\n";
	mess+="set title \"real FID\" ,-1\n";
	mess+="plot \"fid\" u 1:5\n";
mess+="\n";
	mess+="set xlabel \"Hz\"	,1\n";
mess+="\n";
	mess+="set origin .3,.52\n";
	mess+="set title \"real fft\" ,-1\n";
	mess+="plot \"fid\" u 2:3\n";

	mess+="set origin .3,.0\n";
	mess+="set title \"imag fft\" ,-1\n";
	mess+="plot \"fid\" u 2:4\n";

	mess+="set size .3,1\n";
	mess+="set origin .6,0\n";
	mess+="set title \"power fft\" ,-1\n";
	mess+="plot \"fid\" u 2:7\n";

	mess+="set nomultiplot\n";
	mess+="set term postscript\n";
	mess+="set output 'printme.ps'\n";
	mess+="set multiplot\n";
	mess+="set origin 0,0\n";
	mess+="set size 1,1\n";
	mess+="set title \"Bloch Eq Sim\"\n";
	mess+="set nokey\n";
mess+="\n";
	mess+="set multiplot\n";
	mess+="set border 31\n";
mess+="\n";
	mess+="set ylabel \"\"\n";
	mess+="set xlabel \"t\", 1\n";
mess+="\n";
	mess+="set origin 0,0\n";
	mess+="set size .3,.48\n";
	mess+="set title \"imag FID\" ,-1\n";
	mess+="plot \"fid\" u 1:6\n";
mess+="\n";
	mess+="set origin 0,.52\n";
	mess+="set title \"real FID\" ,-1\n";
	mess+="plot \"fid\" u 1:5\n";
mess+="\n";
	mess+="set xlabel \"Hz\"	,1\n";
mess+="\n";
	mess+="set origin .3,.52\n";
	mess+="set title \"real fft\" ,-1\n";
	mess+="plot \"fid\" u 2:3\n";

	mess+="set origin .3,.0\n";
	mess+="set title \"imag fft\" ,-1\n";
	mess+="plot \"fid\" u 2:4\n";

	mess+="set size .3,1\n";
	mess+="set origin .6,0\n";
	mess+="set title \"power fft\" ,-1\n";
	mess+="plot \"fid\" u 2:7\n";

	mess+="set nomultiplot\n";
	mess+="set border 31\n";
	mess+="set xlabel \"\"\n";
	mess+="set zlabel \"\"\n";
	mess+="set ylabel \"\"\n";
	mess+="set title \"\"\n";
	mess+="set size 1,1\n";
	mess+="set origin 0,0\n";
	mess+="set nolabel\n";
#ifdef __CYGWIN__
	mess+="set term windows\n";
#endif

#ifdef ON_WINDOWS
	mess+="set term windows\n";
#else
	mess+="set term x11\n";
#endif
	gn<<mess<<endl;
	gn.close();
}

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

int size=0;
int rank=0;
#ifdef WITH_MPI
	MPI::Init(argc, argv);
	size=MPI::COMM_WORLD.Get_size();
	rank=MPI::COMM_WORLD.Get_rank();
#endif
	

	int q=1;
	double ang=90;
	//int nsGp=2;
	int nsp=0;
	int nsteps=512;
	double tf=1;
	//double zgrad=0;
	double delay=0.003;
	double D=1;
	Vector<double> offset;
	Vector<string> spintype;
	string detsp="31P";
	Vector<double> t2s;
	Vector<double> moles;
	double tr=0;
	int sptypes=3;

	//query_parameter(argc,argv,q++, "Enter Grid size (number of spins is this number cubed): ", nsGp);
	//query_parameter(argc,argv,q++, "Enter Zgradient strength (G/cm) ", zgrad);
	//query_parameter(argc,argv,q++, "Enter Pulse Angle: ", ang);
	query_parameter(argc,argv,q++, "Enter number of spins: ", nsp);
	query_parameter(argc,argv,q++, "Enter Collection points: ", nsteps);

	sptypes=nsp;
	offset.resize(sptypes);
	for( int i=0;i<sptypes;i++){
		query_parameter(argc,argv,q++, "Enter offset for spin "+itost(i)+":", offset[i]);
	}

	t2s.resize(sptypes);
	for( int i=0;i<sptypes;i++){
		query_parameter(argc,argv,q++, "Enter T2 for spin "+itost(i)+":", t2s[i]);
	}
	spintype.resize(sptypes);
	for( int i=0;i<sptypes;i++){
		query_parameter(argc,argv,q++, "Enter spin type (i.e. 1H) for spin "+itost(i)+":", spintype[i]);
	}
	moles.resize(sptypes);
	for( int i=0;i<sptypes;i++){
		query_parameter(argc,argv,q++, "Enter moles for spin "+itost(i)+":", moles[i]);
	}
	query_parameter(argc,argv,q++, "Enter detection spin: ", detsp);
	query_parameter(argc,argv,q++, "Enter delay time: ", delay);
	query_parameter(argc,argv,q++, "Enter Shape Factor: ", D);
	query_parameter(argc,argv,q++, "Enter raditation damping term: ", tr);
	query_parameter(argc,argv,q++, "Enter Final Time: ", tf);
	double pang=90;
	query_parameter(argc,argv,q++, "Pulse Angle: ", pang);

// Bloch set up testing

/*	complex *data, *dataT;
	//data=new complex[16];
	dataT=new complex[32];
	//Vector<complex> dataV(16);
	for(int i=0;i<16;i++){
		//data[i]=2;
		dataT[i]=i;
		dataT[i+16]=2;
		//dataV[i]=2;
	}
	//dataV=FFT(dataV);
	//dataV=IFFT(dataV);
	//FFT1D_(&data[0], 16);
	//IFFT1D_(&data[0],16);
	int dis[2]={16,2};
	FFTN_(dataT, dis,2);
	for(int i=0;i<32;i++){
		cout<<"  "<<dataT[i]<<endl;;
	}

	//delete [] data;
	delete [] dataT;
	*/

	typedef XYZcylinder TheShape;
	typedef XYZshape<TheShape> TheGrid;
	int nsGp=100;
	coord<> mins(-0.01,-0.01,-0.01), maxs(0.01,0.01,0.01);
	coord<int> dims(1,1,nsGp);

	Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims);
	
	Info("Creating inital shape....\n");
	TheShape tester(0,.01,0,PI2, -.01,.01);

	Info("Creating total shape-grid....\n");
	TheGrid jj(tester, gg);
/*	ofstream ou("fidss");
	ou<<jj<<endl;
	//Vector<coord<complex> > outty;
	//Vector<coord<> > lll(jj.size());
	//for(int i=0;i<lll.size();i++){
	//	lll(i)=i+1;
	//}
	//jj.GridFFT(outty, lll );
//	ofstream kkj("opo1");
	//jj.FillVecFromVec(outty, lll);
	//jj.PrintGrid(kkj);
 //	kkj.close();
	//ofstream outgrid("opo");
	//for(int i=0;i<outty.size();i++){
	//	outgrid<<chop(Re(outty(i).x()))<<" "<<chop(Re(outty(i).y()))<<" "<<chop(Re(outty(i).z()))<<endl;
	//}
		//jj.PrintGrid(outgrid);
 	//outgrid.close();

*/
	double zgrad=double(5);
	query_parameter(argc,argv,q++, "Garient Strength: ", zgrad);

	Info("Creating Gradient map grids....\n");
	GradientGrid<TheGrid > mygrid(jj);
	mygrid.G(0,zgrad/2.,zgrad);
	nsp=mygrid.size();


	typedef ListBlochParams< GradientGrid<TheGrid > > MyPars;
	//int nsp=nsGp;
	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n");
	MyPars mypars(nsp, "1H", mygrid);
	

//The pulse list for a real pulse on protons..
	Info("Creating real pulse lists...\n");
	double amp=80000;

	Pulse PP1(spintype[0], amp*PI2, Pi/2.); // (spin, amplitude, phase, offset)
	Pulse PP2(spintype[0], amp*PI2, Pi/2.); // (spin, amplitude, phase, offset)


	if(spintype[0]!=detsp){
		SPulse tmm(detsp, amp*PI2, Pi/2);
		PP2+=tmm;
	}

	Info("Setting spin parameter offsets....\n");
	/*for(int j=0;j<sptypes;j++){
		mypars(j)=spintype[j];
		mypars(j).offset(offset[j]*2.*Pi);
		mypars(j).T1(0);
		mypars(j).T2(t2s[j]);
		mypars(j).Moles(moles[j]);
	}*/
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype[0];
		mypars(j).offset(offset[0]*2.*Pi);
		mypars(j).T1(0);
		mypars(j).T2(t2s[0]);
		mypars(j).Moles(moles[0]);
	}
	mypars.calcTotalMo();
	mypars.SetSolveUnits(Dimensionless);
	mypars.print(cout);
	PP1.print(cout);
	PP2.print(cout);
	double tpulse=PP1.TimeForAngle(pang*Pi/180., spintype[0]);

	ang*=Pi/180.;

//time train testing

	Info("Initializing Time train for first Pulse....\n");
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));


	Info("Initializing Time train for dealy....\n");
	TimeTrain<UniformTimeEngine > D1(UniformTimeEngine(tpulse, tpulse+delay, 10,100));

	Info("Initializing Time train for second Pulse....\n");
	TimeTrain<UniformTimeEngine > P2(UniformTimeEngine(tpulse+delay, tpulse+delay+tpulse, 10,100));

	//int nnss=int(PI2*max(offset, offset2))/100;
	Info("Initializing Time train for FID....\n");
	TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse+delay+tpulse, tf+tpulse+delay+tpulse, nsteps,100));
//Integration testing
	//Info("Pulsing spins....\n");
	//me.DeltaPulse(ang, 0, 0,"1H");


//Extra ineractions
	typedef Interactions<BulkSus, RadDamp> MyInteractions;
	Info("Setting Interactions....\n");
	BulkSus BsRun(D);
	RadDamp RdRun(tr);

	MyInteractions MyInts(BsRun, RdRun);


//typedefs for Bloch parameter sets
	typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
	typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//THis is the BLoch solve to perform a pulse
	Info("Initializing total parameter list with a pulse....\n");
	PulseBloch myparspulse(mypars, PP1, MyInts);

//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n");
	NoPulseBloch me;
	me=(myparspulse);

	Info("Integrating first Pulse (Gradient Off)....\n");
	Vector<coord<> > tm=me.CurrentMag();

	//mypars.GradOff();
	BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");
	drivP.SetWritePolicy(Hold);
	drivP.SetCollectionPolicy(FinalPoint);
	if(!drivP.Solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n");
		return -1;
	}
//cout<<drivP.LastPoint()<<endl;



	Info("\n\nIntegrating delay (gradient Off)....\n");
	//mypars.GradOn();
	BlochSolver<NoPulseBloch > driv(me, drivP.LastPoint(), "out");
	driv.SetWritePolicy(Hold);
	driv.SetCollectionPolicy(FinalPoint);

	/*if(!driv.Solve(D1)){
		Info(" ERROR!!..could not integrate Delay D1....\n");
		return -1;
	}

	Info("\n\nIntegrating second Pulse (gradient Off)....\n");
	myparspulse.SetPulses(PP2);
	drivP.SetInitialCondition(driv.LastPoint());
	if(!drivP.Solve(P2)){
		Info(" ERROR!!..could not integrate pulse P2....\n");
		return -1;
	}
	*/


	Info("\n\nIntegrating for FID (gradient Off)....\n");


	driv.SetWritePolicy(Hold);
	driv.SetCollectionPolicy(FIDonly);
	driv.SetInitialCondition(drivP.LastPoint());
	driv.SetDetect(detsp);

	if(driv.Solve(F1)){
		string fname="fid";
		driv.WriteSpectrum(fname);
		WriteGnu();

	}


	printTime();



#ifdef WITH_MPI
	MPI::Finalize();
#endif


 // Grid testing stuff

	//gg.Translate(10,4,9);
	//gg.Center(50,50,50);
	//o2<<jj;


/*// Gradient Grid testing
	GradientGrid<FullCubeGrid> gg(0, 0.5, 10);
	gg.Gx(60);
	ofstream oo("gout");
	gg.printData(oo, Gradx);
*/


}



