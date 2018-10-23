


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
match_bfield.cc---
	a litte progam to match 'average' Bfields
	across the grid to some value

	Bz/By=tan(acos(1/sqrt(3)))
	Bz/Bx=tan(acos(1/sqrt(3)))
	(want them to form the magic angle)

	and we want the zfield to give us an Larmor frequencey
	for protons of 1 MHz
*/

#include "blochlib.h"
#include "rotatingfield.h"

//MUST include this to use the minimzer
#include "minuit/minuit.h"
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

//here is our MINUIT function wrapper for the interaface
double MagCalc(double *amps, int npar, int iflag);
void fcn (int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil)
{
  *fcnval=MagCalc(xval, npar, iflag);
}


//need to make our 'field' a global var as well as the parameter vars
//want Bz to give 1 MHz field
double ztarget;

//want z-x and z-y to form the amgic angle
double yztarget;
double xztarget;

typedef StepTimeBfield<XYZshape<XYZfull> > MField ;
MField myB;

double MagCalc(double *amps, int npar, int iflag)
{

	char nm[10];
	int  dumm;
	double error, dummy, val;
	myB.Coils[0].amps=amps[0];
	myB.Coils[1].amps=amps[1];
	myB.Coils[2].amps=amps[2];
	myB.calculateField();
	coord<> Boo,avbx, avbz, avby;
//dump some progress info
	if(MPIworld.master()){
		Boo=myB.AverageField(coord<int>(1,1,1) );
		cout<<"Max Field: [ "<<myB.MaxField(coord<int>(1,1,1))<<"]"<<endl;
		cout<<"Average Field: [ "<<Boo<<"]"<<endl;
		cout<<"Main Frequency: "<<GAMMA1H/PI2*Boo.z()*1.e-4<<endl;

		avbx=sum(myB.BxCoil)/myB.BxCoil.size();
		cout<<"Averge Mag from X-Coil: [ "<<avbx<<"]"<<endl;
		avby=sum(myB.ByCoil)/myB.ByCoil.size();
		cout<<"Averge Mag from Y-Coil: [ "<<avby<<"]"<<endl;
		avbz=sum(myB.BzCoil)/myB.BzCoil.size();
		cout<<"Averge Mag from Z-Coil: [ "<<avbz<<"]"<<endl;
		cout<<"currents: ["<<myB.Coils[0].amps<<" "<<myB.Coils[1].amps<<" "<<myB.Coils[2].amps<<"]"<<endl;

		cout<<"xz Angle: "<<Boo.z()/Boo.x()*180.0/PI<<endl;
		cout<<"yz Angle: "<<Boo.z()/Boo.y()*180.0/PI<<endl<<endl;


	}
	MPIworld.scatter(Boo);
	MPIworld.scatter(avbx);
	MPIworld.scatter(avbz);
	MPIworld.scatter(avby);

	double chiz=0, chiy=0, chix=0;
	MNPOUT(1, nm, val, error, dummy, dummy, dumm);
	if(error!=0) chix=tan(Boo.z()/Boo.x())-tan(xztarget);
	MNPOUT(2, nm, val, error, dummy, dummy, dumm);
	if(error !=0) chiy=tan(Boo.z()/Boo.y())-tan(yztarget);
	MNPOUT(3, nm, val, error, dummy, dummy, dumm);
	if(error!=0) chiz=(GAMMA1H/PI2*Boo.z()*1.e-4 -ztarget);


	return sqrt(chiz*chiz+chix*chix+chiy*chiy);
}


int main(int argc, char *argv[])
{

/*** start MPI if we can ***/
	MPIworld.start(argc, argv);

	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	MPIworld.scatter(fn);

	Parameters pset(fn);

	std::string logfn=pset.getParamS("logfile", "", false, "info.log");
	std::ofstream log(logfn.c_str());

/*** Get The grid ***/
	Info("Creating grid....\n",log);

	pset.addSection("grid");
	coord<> mins=pset.getParamCoordD("min", "grid");
	coord<> maxs=pset.getParamCoordD("max", "grid");
	coord<int> dims=pset.getParamCoordI("dim", "grid");
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n",log);
	XYZfull tester;
	Info("Creating total shape-grid....\n",log);
	XYZshape<XYZfull> jj( gg, tester);

/*** Our field type ***/
	coord<int> rotframe=pset.getParamCoordI("rotframe", "", ',',false);

	std::string fielddata=pset.getParamS("fieldfile","",false, "field.mat");

	std::string magsection=pset.getParamS("section","");
	myB.read(jj,pset, magsection);


/*** the fitting input parameters FOR RANK 1 ONLY***/

	pset.addSection("fitting");
	double toler=pset.getParamD("tolerance", "fitting", false, 1e-4);

//want Bz to give 1 MHz field
	ztarget=pset.getParamD("ztarget", "fitting", false, 1.0e6);

//want z-x and z-y to form the amgic angle
	yztarget=pset.getParamD("yztarget", "fitting", false,acos(1/sqrt(3.0))*180.0/PI)*PI/180.0;
	xztarget=pset.getParamD("xztarget", "fitting", false,acos(1/sqrt(3.0))*180.0/PI)*PI/180.0;

//these are the valid 'MINUIT parameter strings
//they should look like ' 15  ''Lambda Mass''  1.2, 0.1'
	std::string zampR=pset.getParamS("zbounds", "fitting");
	std::string xampR=pset.getParamS("xbounds", "fitting");
	std::string yampR=pset.getParamS("ybounds", "fitting");


	int err;
	char *moo=NULL;
	MNINIT(5,6,7);
	moo=const_cast<char *>(xampR.c_str());
	MNPARS(moo, err);
	moo=const_cast<char *>(yampR.c_str());
	MNPARS(moo, err);
	moo=const_cast<char *>(zampR.c_str());
	MNPARS(moo, err);

	int MAXLINE=256;
	char *command; command=new char[MAXLINE];
	snprintf (command, MAXLINE, "MIGRAD");
	MNCOMD (minuitfcn, command, err, NULL);
	delete command;

/*** dump out the good data ***/
	if(MPIworld.master())	myB.writeMatlab(fielddata);

/*** Stop MPI ***/
	MPIworld.end();
}



