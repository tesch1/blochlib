


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 11-27-01
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
A 'chunk' of spins in a Rotating (time dependant) nonuniform magnetic field
uses the classes found in 'rotatingfield.h"

Uses the 'stepping' fields i.e. fields that are switched NOT rotated..
this is more realistic as rotating large fields are next to impossible ot
achieve experimentally (too much current to shunt too quickly...
and of course inductance problems)...so here we use the
'StepTimeBfield' class...

How the initial bfield should be calculated..
1) it assumed that you have 3 orthogonal coils (bx, by,bz) in the
   Parameter file IN THAT ORDER...
2) the field should be tweeked BEFORE full simulation such that
   all 3 fields have the same average strength across the grid

  use 'match_field.cc' to tune you fields

3) if no directions are specified onle BZ will be usecd
4) times and directional switches should go like so...
	 t=0---------t=t1---------t=t2------...
	 |----dir0-----|---dir1----|----dir2...

 	time always starting at 0...


*/

#include "blochlib.h"
#include "rotatingfield.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

timer stopwatch;
void printTime(int nrounds=1){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}

void Info(std::string mess, ostream &oo)
{
	oo<<mess;
	oo.flush();
}


int main(int argc,char* argv[])
{


/********** Parameter File ******/
	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	Parameters pset(fn);

//need the log file name first
	std::string logf=pset.getParamS("logfile","",false, "info.log");
	ofstream log(logf.c_str());

/***** GRIDS ****/
	pset.addSection("grid");

	std::string gtype=pset.getParamS("shape", "grid", false, "full");

	coord<int> dims(pset.getParamCoordI("dim", "grid"));
	coord<> maxs(pset.getParamCoordD("max", "grid"));
	coord<> mins(pset.getParamCoordD("min", "grid"));

	/*

//clyidrical grid if wanted
	Info("Creating grid....\n",log);
	typedef XYZshape<XYZcylinder> TheGrid;
	Grid<UniformGrid> gg(mins, maxs, dims);

	coord<> shamaxs(pset.getParamCoordD("shapemax", "grid", ','));
	coord<> shamins(pset.getParamCoordD("shapemin", "grid", ','));

	Info("Creating inital shape....\n",log);
	XYZcylinder tester(shamins, shamaxs);
	Info("Creating total shape-grid....\n",log);
	XYZshape<XYZcylinder> jj( gg, tester);

	*/


	typedef XYZshape<XYZfull> TheGrid;
	Info("Creating grid....\n",log);
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n",log);
	XYZfull tester;
	Info("Creating total shape-grid....\n",log);
	XYZshape<XYZfull> jj( gg, tester);


/**** initial condition ****/
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

	coord<int> rotframe=pset.getParamCoordI("rotframe", "", ',',false);

/**** out file names and log files ****/

	std::string fout=pset.getParamS("fileout");
	std::string magout=pset.getParamS("magout");
	std::string fielddata=pset.getParamS("fieldfile","",false, "field.mat");

	double amp=pset.getParamD("amp");

	int cv=pset.getParamI("lyps");


/**** Bloch Parameters ***/
	string spintype=pset.getParamS("spintype");
	string detsp=spintype;
	double moles=pset.getParamD("moles");

	int nsp=jj.size();

//typdefs for parameters
	typedef ListBlochParams<TheGrid,BPoptions::Particle|BPoptions::HighField, coord<> > MyPars;

	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n",log);
	MyPars mypars(nsp, "1H", jj, IC);
	nsp=mypars.size();

/*** Pulse ***/
	Info("Creating real pulse lists...\n",log);

	double pang=pset.getParamD("pang");

	Pulse PP1(spintype, amp*PI2, Pi/2.); // (spin, amplitude, phase, offset)


	Info("setting spin parameter moles and Spin Types....\n",log);
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles/double(nsp));
	}

	mypars.calcTotalMo();
	PP1.print(log);
	double tpulse=PP1.timeForAngle(pang*Pi/180., spintype);

	pang*=Pi/180.;

/****** Time Train for 90 deg pulse *****/
	int nsteps=pset.getParamI("npts");
	double tf=pset.getParamD("tf");

	Info("Initializing Time train for first Pulse....\n",log);
	TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

	Info("Initializing Time train for FID....\n",log);
	TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,20));

/****** Interactions *****/
	Info("setting Interactions....\n",log);
	typedef StepTimeBfield<TheGrid> MField ;
	typedef Offset<MField> Offsets;
	typedef Relax<MyPars::Offset_T> Relaxs;
	typedef TanhScale Scaler;
	typedef DimLessDipole<TheGrid> DipDipole;
	typedef Interactions<Relaxs, Offsets, DipDipole, BulkSus, RadDamp> MyInteractions;

	double D=pset.getParamD("bulksus");
	double offset=0; offset=pset.getParamD("offset");
	double t2s=pset.getParamD("T2");
	double t1s=pset.getParamD("T2");
	double tr=pset.getParamD("raddamp");
	double dipstr=pset.getParamD("dipole_str")*PI2;

//make sure we do no tput in infinity for the relaxation
	Relaxs RelRun(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);

//field coils
	std::string magsection=pset.getParamS("section","");
	MField myB(jj,pset, magsection);

	myB.calculateField();

//find the max B (in either direction)
// such that we can sit in its rotating frame
	coord<> Boo=myB.AverageField(rotframe);
	log<<"Max Field: [ "<<myB.MaxField(rotframe)<<"]"<<endl;
	log<<"Average Field: [ "<<Boo<<"]"<<endl;

//must 'reset' the magentic field to the 'rotating frame'
	Range All(Range::Start,Range::End);
	coord<> avbx=sum(myB.BxCoil)/myB.BxCoil.size();
	coord<> avby=sum(myB.ByCoil)/myB.ByCoil.size();
	coord<> avbz=sum(myB.BzCoil)/myB.BzCoil.size();
	log<<"Averge Mag from X-Coil: [ "<<avbx<<"]"<<endl;
	log<<"Averge Mag from Y-Coil: [ "<<avby<<"]"<<endl;
	log<<"Averge Mag from Z-Coil: [ "<<avbz<<"]"<<endl;
	//myB.BxCoil(SolverOps::All)-=avbx;
	//myB.ByCoil(SolverOps::All)-=avby;
	//myB.BzCoil(SolverOps::All)-=avbz;

	myB.writeMatlab(fielddata);


//get the pulse field parameters chunk
	pset.addSection("pulsefield");
	myB.readPulseParams(pset);

//We take the 'norm' of the field at t=0.0, and scale
// the field by this amount as numerical stability
// at 1Mhz freq is not good, so we scale time
	double scaleB=myB.norm(0.0);
	myB.BxCoil/=scaleB;
	myB.ByCoil/=scaleB;
	myB.BzCoil/=scaleB;
	log<<"Time and Bfield Scaled By..."<<scaleB<<endl;

	myB.writeMatlab("scaledf.mat");

//pulsefield list (as we want initial condition to be Bfield(t=0) )
	Offsets OffRun(myB, mypars);

//display the pulse set
	for( int i=0;i<myB.tsplit.size();++i){
		log<<myB.tsplit[i]<<"\t"<<myB.direction[i]<<endl;
	}

//remaining interactions
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	DipDipole DipDip(jj,dipstr);
	DipDip.Dynamic=true;
	if(dipstr==0){ DipDip.Dynamic=false; DipDip.off();	}

	Info("Grid Point.....................Magnetic Field.................Offset.........................\n",log);
	TheGrid::iterator myIt(jj);
	while(myIt){
		log<<myIt.Point()<<" ["<<myB.Bfield(myIt.curpos())<<"] ["<<OffRun.offset(myIt.curpos())/PI2<<"] "<<endl;
		++myIt;
	}


	MyInteractions MyInts(RelRun,OffRun, DipDip, BsRun, RdRun);

/**** The Massive Bloch containers ***/
//typedefs for Bloch parameter sets
	typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
	typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

//THis is the BLoch solve to perform a pulse
	Info("Initializing total parameter list with a pulse....\n",log);
	PulseBloch myparspulse(mypars, PP1, MyInts);
	if(cv) myparspulse.calcVariational();

//The offset class sets the corect Bo to each of the points
//The 'Bloch' sets up the correct Initial Condition
// so now we print the initial starter condition
	mypars.print(log);

//This is the Bloch solver to Collect the FID (i.e. has no pusles...FASTER)
	Info("Initializing total parameter list for FID collection....\n",log);
	NoPulseBloch me;
	me=(myparspulse);
	Info("Integrating first Pulse....\n",log);
	Vector<coord<> > tm=me.currentMag();

	stopwatch.reset();
	BlochSolver<PulseBloch > drivP(myparspulse, tm);
	drivP.setWritePolicy(SolverOps::Hold);
	drivP.setCollectionPolicy(SolverOps::FinalPoint);
	if(!drivP.solve(P1)){
		Info(" ERROR!!..could not integrate pulse P1....\n",log);
		return -1;
	}

	BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint(), "out");

	Info("\nIntegrating FID ....\n",log);

	std::string lypname="lyps";

	driv.setWritePolicy(SolverOps::Hold);
	driv.setCollectionPolicy(SolverOps::MagAndFID);
	driv.setDetect(detsp);
	if(cv) {
		driv.setLyapunovPolicy(SolverOps::LypContinous);
		driv.setLypDataFile(lypname);
	}

	if(driv.solve(F1)){
		driv.writeSpectrum(fout);
		driv.writeMag(magout);
		if(cv) WriteGnuplotLyp(lypname, nsp*3);
	}


	printTime();
	return 0;

}



