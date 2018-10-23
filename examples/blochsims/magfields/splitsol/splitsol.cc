


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-30-02
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
splitsol.cc---
	a program that calculates the fiel inside
	a split solonoid, then attempts to calculate
	the 'time' for a 90 pulse (the inhomegnatity profile)
*/

#include "blochlib.h"

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

class myBFcalc: public StaticField
{

	public:
		typedef  double Offset_T;
		Vector<coord<> > Bf;

		Offset_T Bfield(int i)
		{ return sqrt(Bf[i].x()*Bf[i].x()+ Bf[i].y()*Bf[i].y());	}

		Offset_T Bfield( double t, int i)
		{ return Bfield(i);	}
};

//The Grids
typedef XYZfull TheShape ;
typedef XYZshape<TheShape> TheGrid;

//the magnetic field caluclator
typedef MultiBiot<TheGrid > MField;

typedef ListBlochParams< TheGrid, BPoptions::Particle|BPoptions::HighField, double > MyPars;

//Extra ineractions
typedef Offset<myBFcalc> Offsets;
typedef Interactions<Offsets > MyInteractions;

//typedefs for Bloch parameter sets
typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;



int main(int argc, char *argv[])
{

	MPIworld.start(argc, argv);
	std::string parse="";
	int q=1;
	if(MPIworld.master())
		query_parameter(argc, argv, q++, "\n\tEnter Parmeter set file name:", parse);

	MPIworld.scatter(parse);

	Parameters pset(parse);
	pset.addSection("parameters");
	std::string choose=pset.getParamS("section","parameters");
	typedef XYZshape<XYZfull> TheGrid;
	pset.addSection("grid");
	Grid<UniformGrid> g1(pset.getParamCoordD("min", "grid"),
	pset.getParamCoordD("max", "grid"),
	pset.getParamCoordI("dim", "grid"));
	TheGrid g2(g1, XYZfull());
	MField mycoil(g2,pset, choose);
	mycoil.Controller=MPIworld;
	if(mycoil.Bfield.size()==0) mycoil.calculateField();

	if(MPIworld.master()){
		mycoil.writeMatlab(pset.getParamS("matout", "parameters", false, "field.mat"));
		mycoil.write(pset.getParamS("textout", "parameters", false, "shape.boit"));
	}


	//now attempt to simulate the 'offsets' in the xyplane
	if(MPIworld.master()){
		int nsteps=pset.getParamI("npts", "parameters");
		double tf=pset.getParamD("tf", "parameters");
		double offset=0; offset=pset.getParamD("offset", "parameters");
		string spintype=pset.getParamS("spintype", "parameters", false, "1H");
		string fout=pset.getParamS("fout", "parameters", false, "data");
		string mout=pset.getParamS("mout", "parameters", false, "mag");
		string detsp=spintype;

	//parameters
		MyPars mypars(g2.size(), spintype, g2);
		mypars.calcTotalMo();
		mypars.print(std::cout);

	//the offset interaction with the field
		coord<> maxF=mycoil.maxField(), aveF=mycoil.averageField();
		cout<<"Max Field: [ "<<maxF<<"]"<<endl;
		cout<<"Average Field: [ "<<aveF<<"]"<<endl;
		myBFcalc bfc;
		bfc.Bf=mycoil.Bfield;
		Range All(Range::Start, Range::End);
		bfc.Bf(SolverOps::All)-=aveF;

		Offsets OffRun(bfc, mypars);
		MyInteractions inters(OffRun);

		TimeTrain<UniformTimeEngine > F1(0, tf, nsteps,1);

		NoPulseBloch me(mypars, inters);
		Vector<double> ppars=pset.getParamVectorD("pulse", "parameters")*DEG2RAD;
		if(ppars.size()<2){
			cerr<<endl<<"pulse must have 2 inputs {algne},{phase}"<<endl;
			return 0;
		}
		me.deltaPulse(ppars[0], ppars[1]);

		BlochSolver<NoPulseBloch > drivD(me, me.currentMag());

	//output trajectory data if wanted
		drivD.setCollectionPolicy(SolverOps::MagAndFID);
		drivD.setWritePolicy(SolverOps::Hold);

	//solve the FID and write it to a file
		if(drivD.solve(F1)){
			drivD.writeSpectrum(fout);
			drivD.writeMag(mout);
		}

		printTime();

	//ring a bell when we are done
		std::cout<<"\a"<<std::endl;

	}


	MPIworld.end();
}



