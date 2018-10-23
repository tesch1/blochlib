


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 12-30-01
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
	if(MPIworld.master())
	{	oo<<mess;	oo.flush();}
}




//why this function?  becuase it is the easiest way to
// allow the usr to choose multiple grid types in the config file
// (i.e. Cylinders vs Rectangles, etc)
template<class TheGrid>
class Dooer{
	public:
		Dooer(){};

		int run(TheGrid &jj, Parameters &pset, ostream &log)
		{
			int nsteps=pset.getParamI("npts");
			double tf=pset.getParamD("tf");
			double D=pset.getParamD("bulksus");
			double offset=0; offset=pset.getParamD("offset");
			string spintype=pset.getParamS("spintype");
			string detsp=spintype;
			double t2s=pset.getParamD("T2");
			double t1s=pset.getParamD("T1");
			double moles=pset.getParamD("moles");
			double tr=pset.getParamD("raddamp");
			double dipstr=pset.getParamD("dipole_str")*PI2;

			double pang=pset.getParamD("pang");

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

			coord<int> rotframe=pset.getParamCoordI("rotframe", "", ',',false);

			std::string fout=pset.getParamS("fileout");
			std::string magout=pset.getParamS("magout");
			std::string fielddata=pset.getParamS("fieldfile","",false, "field.mat");

			double amp=pset.getParamD("amp");

			int cv=pset.getParamI("lyps");
		// Bloch set up testing

			int nsp=jj.size();
			typedef ListBlochParams<TheGrid,BPoptions::Particle|BPoptions::HighField, coord<> > MyPars;
			//int nsp=nsGp;
			Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n",log);
			MyPars mypars(nsp, "1H", jj, IC);
			nsp=mypars.size();

		//The pulse list for a real pulse on protons..
			Info("Creating real pulse lists...\n",log);

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

		//time train testing

			Info("Initializing Time train for first Pulse....\n",log);
			TimeTrain<UniformTimeEngine > P1(UniformTimeEngine(0., tpulse, 10,100));

			Info("Initializing Time train for FID....\n",log);
			TimeTrain<UniformTimeEngine > F1(UniformTimeEngine(tpulse, tpulse+tf, nsteps,20));

		//Extra ineractions
			typedef WrBfield<TheGrid> MField ;
			typedef Offset<MField> Offsets;
			typedef Relax<typename MyPars::Offset_T> Relaxs;
			typedef DimLessDipole<TheGrid> DipDipole;

			typedef Interactions<Relaxs, Offsets, DipDipole, BulkSus, RadDamp> MyInteractions;
			Info("setting Interactions....\n",log);

			Relaxs RelRun(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);

		//calclate B fields
			std::string magsection=pset.getParamS("section","");
			MField myB(jj,pset, magsection);
			myB.calculateField();

		//find the max B (in either direction)
		// such that we can sit in its rotating frame
			coord<> Boo=myB.AverageField(rotframe);
			log<<"Max Field: [ "<<myB.MaxField(rotframe)<<"]"<<endl;
			log<<"Average Field: [ "<<Boo<<"]"<<endl;

		//must 'reset' the magentic field to the 'rotating frame'
			Range all(Range::Start,Range::End);
			coord<> avbx=sum(myB.BxCoil)/myB.BxCoil.size();
			coord<> avby=sum(myB.ByCoil)/myB.ByCoil.size();
			coord<> avbz=sum(myB.BzCoil)/myB.BzCoil.size();
			myB.BxCoil(all)-=avbx;
			myB.ByCoil(all)-=avby;
			myB.BzCoil(all)-=avbz;

		//dump the field info to a matlab file
			 if( MPIworld.master()) myB.writeMatlab(fielddata);

		//get the rotation speed of the fields..
			myB.wr=pset.getParamCoordD("fieldwr", "", false);

			Offsets OffRun(myB, mypars);


			BulkSus BsRun(D);
			RadDamp RdRun(tr);
			DipDipole DipDip(jj,dipstr);

			Info("Grid Point.....................Magnetic Field.................Offset.........................\n",log);
			if( MPIworld.master()) {
				typename TheGrid::iterator myIt(jj);
				while(myIt){
					log<<myIt.Point()<<" ["<<myB.Bfield(myIt.curpos())<<"] ["<<OffRun.offset(myIt.curpos())/PI2<<"] "<<endl;
					++myIt;
				}
			}


			MyInteractions MyInts(RelRun,OffRun, DipDip, BsRun, RdRun);


		//DipDip.Off();

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


			if(dipstr==0){	DipDip.Off(); DipDip.Dynamic=false;	}


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

		//solve for the 90 pulse
			if(!drivP.solve(P1)){
				Info(" ERROR!!..could not integrate pulse P1....\n",log);
				return -1;
			}


			BlochSolver<NoPulseBloch > driv(me, drivP.lastPoint());

			Info("\nIntegrating FID ....\n",log);

			std::string lypname="lyps";

			driv.setWritePolicy(SolverOps::Hold);

		//see if we want to write the trajectories or not...
			std::string dataou=pset.getParamS("trajectories", "", false);
			std::ofstream trajout;
			if(dataou!="" && MPIworld.master() ){
				trajout.open(dataou.c_str());
				driv.setCollectionPolicy(SolverOps::All);
			}else{
				driv.setCollectionPolicy(SolverOps::MagAndFID);
			}

			driv.setDetect(detsp );
			if(cv && MPIworld.master()) {
				driv.setLyapunovPolicy(SolverOps::LypContinous);
				driv.setLypDataFile(lypname);
			}

			if(driv.solve(F1)){
				if(MPIworld.master())
				{
					driv.writeSpectrum(fout);
					driv.writeMag(magout);
					if(dataou!="")	driv.writeData(trajout);
					if(cv) WriteGnuplotLyp(lypname, nsp*3);
				}
			}


			printTime();
			return 0;
	}
};


int main(int argc,char* argv[]){

	MPIworld.start(argc, argv);

	std::string fn;
	if(MPIworld.master())
		query_parameter(argc,argv,1, "Enter file to parse: ", fn);

	MPIworld.scatter(fn);
	Parameters pset(fn);

	pset.addSection("grid");

	std::string gtype=pset.getParamS("shape", "grid", false, "full");
	std::string logf=pset.getParamS("logfile","",false, "info.log");
	ofstream log(logf.c_str());

	coord<int> dims(pset.getParamCoordI("dim", "grid"));
	coord<> maxs(pset.getParamCoordD("max", "grid"));
	coord<> mins(pset.getParamCoordD("min", "grid"));

//oi i say.....memory problem in compilation of more then one grid type..
//alas,


	/*if(gtype=="cylinder"){
		Info("Creating grid....\n",log);
		Grid<UniformGrid> gg(mins, maxs, dims);

		coord<> shamaxs(pset.getParamCoordD("shapemax", "grid", ','));
		coord<> shamins(pset.getParamCoordD("shapemin", "grid", ','));

		Info("Creating inital shape....\n",log);
		XYZcylinder tester(shamins, shamaxs);
		Info("Creating total shape-grid....\n",log);
		XYZshape<XYZcylinder> jj( gg, tester);
		Dooer<XYZshape<XYZcylinder> > doo;
		return doo.run(jj, pset,log);

	}else if(gtype=="rect"){
		Info("Creating grid....\n",log);
		Grid<UniformGrid> gg(mins, maxs, dims);

		coord<> shamaxs(pset.getParamCoordD("shapemax", "grid", ','));
		coord<> shamins(pset.getParamCoordD("shapemin", "grid", ','));

		Info("Creating inital shape....\n",log);
		XYZrect tester(shamins, shamaxs);
		Info("Creating total shape-grid....\n",log);
		XYZshape<XYZrect> jj( gg, tester);
		Dooer<XYZshape<XYZrect> > doo;
		return doo.run(jj, pset,log);
	}else{*/
		Info("Creating grid....\n",log);
		Grid<UniformGrid> gg(mins, maxs, dims);

		Info("Creating inital shape....\n",log);
		XYZfull tester;
		Info("Creating total shape-grid....\n",log);
		XYZshape<XYZfull> jj( gg, tester);
		Dooer<XYZshape<XYZfull> > doo;

		int mplo=doo.run(jj, pset,log);
		MPIworld.end();
		return mplo;
}



