

#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*
A simple 90 pulse demonstrating the
effect of radiation dmaping on a few spins
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
    string spintype=pset.getParamS("spintype");
    string detsp=spintype;
    double t2s=pset.getParamD("T2");
    double t1s=pset.getParamD("T1");

    double moles=pset.getParamD("moles");

    std::string fout=pset.getParamS("fidout");
    std::string magout=pset.getParamS("magout");

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
               BPoptions::Particle | BPoptions::HighField,
               double > MyPars;
    int nsp=jj.size();
    Info("Creating entire spin parameter list for "+itost(nsp)+" spins....");
    MyPars mypars(nsp, "1H", jj);
    nsp=mypars.size();

//The pulse list for a real pulse on protons..
    Info("Creating real pulse lists...");

//get the info from the pset
    double pang=pset.getParamD("pulseangle");
    double amp=pset.getParamD("pulseamp");
    double phase=pset.getParamD("pulsephase");

// (spin, amplitude, phase, offset)
    Pulse PP1(spintype, amp*PI2, phase*DEG2RAD);

//get the Bo
    double inBo=pset.getParamD("Bo");

    Info("Setting spin parameter offsets....");
    for(int j=0;j<nsp;j++){
        mypars(j)=spintype;
        mypars(j).moles(moles);
        mypars(j).Bo(inBo);
        mypars.temperature(inTemp);
    }

    mypars.calcTotalMo();
    mypars.print(cout);
    PP1.print(cout);

//get the time for the 90 pulse
    double tpulse=PP1.timeForAngle(pang*Pi/180., spintype);

//the time trains
    Info("Initializing Time train for first Pulse....");
    TimeTrain<UniformTimeEngine >
       P1(UniformTimeEngine(0., tpulse, 10,100));

    Info("Initializing Time train for FID....");
    TimeTrain<UniformTimeEngine >
       F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,1));

//Extra interactions
    typedef Interactions<
                 Offset<>,
                 Relax<>,
                 RadDamp > MyInteractions;
    Info("Setting Interactions....");

//the offsets
//get the first offset
    double offset=pset.getParamD("offset")*PI2;
    Offset<> myOffs(mypars, offset);

//Relaxation
    Relax<> myRels(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);

//Bulk suseptibility
    double tr=pset.getParamD("tr");

    RadDamp myRD(tr);

//total interaction obect
    MyInteractions MyInts(myOffs, myRels, myRD);

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
    Info("Integrating first Pulse....\n");

//our initial condition
    Vector<coord<> > tm=me.currentMag();
	tm.print(cout, "\n");
    stopwatch.reset();
    BlochSolver<PulseBloch > drivP(myparspulse, tm, "out");
//cout<<drivP.lastPoint()<<endl;

//output trajectory data if wanted
    std::ofstream trajout;
    if(dataou!=""){
        trajout.open(dataou.c_str());
        drivP.setCollectionPolicy(SolverOps::All);
    }else{
        drivP.setCollectionPolicy(SolverOps::FinalPoint);
    }

//integrate the Pulse
    drivP.setWritePolicy(SolverOps::Hold);
    if(!drivP.solve(P1)){
        Info(" ERROR!!..could not integrate pulse P1....\n");
        return -1;
    }else if(dataou!=""){
        drivP.writeData(trajout);
    }

//the fids initial condition is just the previous
// integrations last point
    BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());

    Info("\nIntegrating FID ....\n");

//set various data collection policies
    std::string lypname="lyps";
    if(dataou!=""){
        driv.setCollectionPolicy(SolverOps::All);
    }else{
        drivP.setCollectionPolicy(SolverOps::MagAndFID);
    }
    driv.setWritePolicy(SolverOps::Hold);

//set the detection spin
    driv.setDetect(detsp);

//integrate the FID
    if(driv.solve(F1)){

    //dump out the data
        driv.writeSpectrum(fout);
        driv.writeMag(magout);
        if(dataou!="") driv.writeData(trajout);
    }


    printTime();

}

