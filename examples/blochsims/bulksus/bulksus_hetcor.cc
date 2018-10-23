

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*
THis simulates the effect of the Bulk Suseptibility on
a HETCOR experiement...hopefully we shall see several echos
in the indirect dimension

a HETOCR is a 2D experiement

spin1:: 90--t--90-----
spin2:: -------90-FID

*/
timer stopwatch;
void printTime(int nrounds=1){
     std::cout <<std::endl<< "Time taken: "
     << (stopwatch()/nrounds) << " seconds";
}

void Info(std::string mess)
{
    cout<<mess<<endl;
    cout.flush();
}

int main(int argc,char* argv[])
{
    std::string fn;

//the parameter file
    query_parameter(argc,argv,1, "Enter file to parse: ", fn);
    Parameters pset(fn);

//get the basic parameters
    int nsteps=pset.getParamI("npts");
    double tf=pset.getParamD("tf");
    double inTemp=pset.getParamD("temperature");
    string spintype1=pset.getParamS("spintype1");
    string spintype2=pset.getParamS("spintype2");
    string detsp=pset.getParamS("detect");;

    double moles=pset.getParamD("moles");

    std::string fout=pset.getParamS("fidout");

    coord<int> dims(pset.getParamCoordI("dim"));
    coord<> mins(pset.getParamCoordD("min"));
    coord<> maxs(pset.getParamCoordD("max"));

    std::string dataou=pset.getParamS("trajectories", "", false);

// Grid Set up
    typedef XYZfull TheShape;
    typedef XYZshape<TheShape> TheGrid;

    Info("Creating grid....");
    Grid<UniformGrid> gg(mins, maxs, dims);

    Info("Creating inital shape....");
    TheShape tester;
    Info("Creating total shape-grid....");
    TheGrid jj( gg, tester);

//List BlochParameters
    typedef ListBlochParams<
               TheGrid,
               BPoptions::Density | BPoptions::HighField,
               double > MyPars;
    int nsp=jj.size();
    Info("Creating entire spin parameter list for "+itost(nsp)+" spins....");
    MyPars mypars(nsp, "1H", jj);
    nsp=mypars.size();

//The pulse list for a real pulse on protons..
    Info("Creating real pulse lists...");

//get the info from the pset
    coord<> pang1=pset.getParamCoordD("pulse1");
    coord<> pang2=pset.getParamCoordD("pulse2");
    double delaystep=pset.getParamD("delay");

// (spin, amplitude, phase, offset)
    Pulse PP1(spintype1, pang1[2]*PI2, pang1[1]*DEG2RAD);
	Pulse PP2(spintype1, pang2[2]*PI2, pang2[1]*DEG2RAD);
	PP2+=Pulse(spintype2, pang2[2]*PI2, pang2[1]*DEG2RAD);

//get the Bo
    double inBo=pset.getParamD("Bo");

    Info("Setting spin parameter offsets....");
    for(int j=0;j<nsp;j++){
        if(j%2==0){ mypars(j)=spintype1; }
        else{  mypars(j)=spintype2;}

        mypars(j).moles(moles/nsp);
    	mypars(j).Bo(inBo);
        mypars.temperature(inTemp);
    }

    mypars.calcTotalMo();
    mypars.print(cout);
    PP1.print(cout);
    PP2.print(cout);


//Extra interactions
    typedef Interactions<
                 Offset<>,
                 Relax<>,
                 BulkSus > MyInteractions;
    Info("Setting Interactions....");

//the offsets
//get the first offset
    double offset1=pset.getParamD("offset1")*PI2;
    double offset2=pset.getParamD("offset2")*PI2;
    Offset<> myOffs(mypars, offset1);

//Relaxation
    double t2s1=pset.getParamD("T2_1");
    double t1s1=pset.getParamD("T1_1");
    double t2s2=pset.getParamD("T2_2");
    double t1s2=pset.getParamD("T1_2");
    Relax<> myRels(mypars, (!t2s1)?0.0:1.0/t2s1, (!t1s1)?0.0:1.0/t1s1);

    for(int i=0;i<nsp;++i){
		//set the offsets and relaxtion vals
		if(i%2==0){
			myOffs.offset(i)=offset1;
			myRels.T1(i)=(!t1s1)?0.0:1.0/t1s1;
			myRels.T2(i)=(!t2s1)?0.0:1.0/t2s1;
		}else{
			myOffs.offset(i)=offset2;
			myRels.T1(i)=(!t1s2)?0.0:1.0/t1s2;
			myRels.T2(i)=(!t2s2)?0.0:1.0/t2s2;
		}
	}


//Bulk suseptibility
    double D=pset.getParamD("D");

    BulkSus myBs(D);

//total interaction obect
    MyInteractions MyInts(myOffs, myRels, myBs);

//typedefs for Bloch parameter sets
    typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
    typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//second dimension points
	int npts2D=pset.getParamI("npts2D");

//our data matrix
    matrix FIDs(npts2D, nsteps);

//get the time for the 2 90 pulse
	double tpulse1=PP1.timeForAngle(pang1[0]*Pi/180., spintype1);
	double tpulse2=PP2.timeForAngle(pang2[0]*Pi/180., spintype1);

//the time trains this one will always be the same
	Info("Initializing Time train for first Pulse....");
	TimeTrain<UniformTimeEngine >
	   P1(UniformTimeEngine(0., tpulse1, 10,10));

//loop over all our D values
    for(int kk=0;kk<npts2D;++kk){

		double curDelay=double(kk)*delaystep;
		cout<<"On delay: "<<curDelay<<" "<<kk<<"/"<<npts2D<<"    \r"; cout.flush();
	//the time trains for the dealy
		TimeTrain<UniformTimeEngine >
		   D1(UniformTimeEngine(tpulse1, tpulse1+curDelay, 10,5));

	//the time trains for the dealy
		TimeTrain<UniformTimeEngine >
		   P2(UniformTimeEngine(tpulse1+curDelay, tpulse2+tpulse1+curDelay, 10,10));

		TimeTrain<UniformTimeEngine >
		   F1(UniformTimeEngine(tpulse2+tpulse1+curDelay, tpulse2+tpulse1+curDelay+tf, nsteps,5));


    //THis is the BLoch solve to perform a pulse
        PulseBloch myparspulse(mypars, PP1, MyInts);

    //This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
        NoPulseBloch me;
        me=(myparspulse);

    //out initial condition
        Vector<coord<> > tm=me.currentMag();

        stopwatch.reset();
        BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");
		drivP.setProgressBar(SolverOps::Off);

    //integrate the Pulse
        drivP.setWritePolicy(SolverOps::Hold);
        if(!drivP.solve(P1)){
            Info(" ERROR!!..could not integrate pulse P1....");
            return -1;
        }

    //the fids initial condition is just the previous
    // integrations last point
        BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());
		driv.setProgressBar(SolverOps::Off);

    //integrate the Delay
        driv.setWritePolicy(SolverOps::Hold);
        if(!driv.solve(D1)){
            Info(" ERROR!!..could not integrate delay D1....");
            return -1;
        }

    //integrate second the Pulse
       	drivP.setWritePolicy(SolverOps::Hold);
	//set the new pulse set
		myparspulse.setPulses(PP2);
       	drivP.setInitialCondition(driv.lastPoint());
        if(!drivP.solve(P2)){
            Info(" ERROR!!..could not integrate pulse P2....");
            return -1;
        }


    //set the detection spin
        driv.setDetect(detsp);
    //set various data collection policies
       	driv.setInitialCondition(drivP.lastPoint());
        driv.setCollectionPolicy(SolverOps::MagAndFID);
        driv.setWritePolicy(SolverOps::Hold);

	 //integrate the FID
        if(driv.solve(F1)){
            FIDs.putRow(kk, driv.FID());
        }
    }

    matstream matout(fout, ios::binary | ios::out);
    matout.put("vdat", FIDs);
    matout.close();
    printTime();

}

