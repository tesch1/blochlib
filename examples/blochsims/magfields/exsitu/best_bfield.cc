


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
best_bfield.cc---
	a litte progam to find the best geometry/currents
	for a potential design of an exsitu magnet
	the requirements are a large region of linear
	inhomogenaity to within
*/

#include "blochlib.h"

//MUST include this to use the minimzer
#include "minuit/minuit.h"

//the mag fitter parameters container to
// make setting these params much easier
#include "magfitter.h"

//The Dcoil Shape
#include "Dcoil.h"

//include the 'exsitu Sim'
#include "exsituSim.h"
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


/***Global Variables and Definitions ***/
std::ofstream logFile;	//Log File out stream

//need to make our 'field' a global var as well as the parameter vars
//want Bz to give 1 MHz field
double ztarget;

//want z-x and z-y to form the amgic angle
double yztarget;
double xztarget;

//The data container for the Magnetic field fitter data
MagFitter myFitter;

//The Grids
typedef XYZrect TheShape ;
typedef XYZshape<XYZrect> TheGrid;
Vector<TheGrid> myGrids;

//the magnetic field caluclator
typedef MultiBiot<TheGrid > MField;
Vector<MField> myFields;
MField *curB;

/*** Master 'Fitter' function ***/
double MagCalc(double *amps, int npar, int iflag);
double MagCalc(double *amps, int npar, int iflag)
{

	double  val;
	myFitter.setCoilParams(amps); // set the new values in the parameters
	curB->read(myFitter.coilParams(),"", false); //reread the coils
	coord<> Boo(0.0,0.0,0.0);
	if(MPIworld.master()) logFile<<myFitter<<endl;
	coord<> totGsize=0;
	for(int i=0;i<myGrids.size();++i){
		curB->grid=myGrids[i]; //set the grid of the Bfield Calculator
		totGsize+=myGrids[i].dim();
		curB->calculateField();	//calculate the new field
		if(MPIworld.master())
			Boo+=curB->AverageField()*myGrids[i].dim();
	}

	MPIworld.scatter(Boo);

	double chiz=0, chiy=0, chix=0;
	Boo/=totGsize;
	if(MPIworld.master()){
		//logFile<<"Max Field: [ "<<myB.MaxField(coord<int>(1,1,1))<<"]"<<endl;
		logFile<<"Average Field: [ "<<Boo<<"]"<<endl;
		logFile<<"Main Frequency: "<<GAMMA1H/PI2*Boo.z()*1.e-4<<endl;
		logFile<<"xz Angle: "<<Boo.z()/Boo.x()*180.0/PI<<endl;
		logFile<<"yz Angle: "<<Boo.z()/Boo.y()*180.0/PI<<endl<<endl;
	}
	if(myFitter.used(0)) chix=tan(Boo.z()/Boo.x())-tan(xztarget);
	if(myFitter.used(1)) chiy=tan(Boo.z()/Boo.y())-tan(yztarget);
	if(myFitter.used(2)) chiz=(GAMMA1H/PI2*Boo.z()*1.e-4 -ztarget);


	return sqrt(chiz*chiz+chix*chix+chiy*chiy);
}

//here is our MINUIT function wrapper for the interaface
void fcn (int npar, double* grad, double* fcnval, double* xval, int iflag, void* futil)
{
  *fcnval=MagCalc(xval, npar, iflag);
}


int main(int argc, char *argv[])
{

/** Insert the Dcoils into the Coil list ***/
	BiotFunctions.insert("Dhelix", Biot_Dcoil);
	BiotFunctions.insert("Dcircle", Biot_Dcircle);

/*** start MPI if we can ***/
	MPIworld.start(argc, argv);

	std::string fn;
	query_parameter(argc,argv,1, "Enter file to parse: ", fn);
	MPIworld.scatter(fn);

//master parameter set
	Parameters pset(fn);

	pset.addSection("params"); //auxillary parameters

/** log file name ***/
	std::string logfn=pset.getParamS("logfile", "params", false, "output.log");
	if(logfn!="") logFile.open(logfn.c_str());

/*** Get The grid ***/
	Info("Creating grids....\n",logFile);
//here we can have a couple of grids for slice planes
//planes of out 3D world
	pset.addSection("grids");
	Parameters subp(pset.section("grids"));
	std::string subbase=subp.getParamS("base","",false, "grid");
	int maxFit=subp.getParamI("numgrids","",false, 1000000);
	int numGrids=0;
//count the number of params present
	while( (numGrids+1)<=maxFit && subp.addSection(subbase+itost(numGrids+1)))
	{	numGrids++;		}

//add a parameter to our master list
	if(MPIworld.master()) cout<<"Number of Grids: "<<numGrids<<endl;
	if(numGrids==0){
		std::cerr<<"ERROR: Number of grids must be >=1"<<std::endl;
		exit(1);
	}
	myGrids.resize(numGrids);
	numGrids=0;
	while(numGrids<myGrids.size())
	{
		std::string secName=subbase+itost(numGrids+1);
		coord<> mins=subp.getParamCoordD("min", secName);
		coord<> maxs=subp.getParamCoordD("max", secName);
		coord<int> dims=subp.getParamCoordI("dim", secName);
		Grid<UniformGrid> gg(mins, maxs, dims);
		myGrids[numGrids]=TheGrid(gg, TheShape(mins, maxs));
		++numGrids;
	}

/*** Our field type ***/
	pset.addSection("fitter");
	std::string CoilSec=pset.getParamS("section", "params");
	Vector<std::string> CoilList=parse_param(CoilSec, ',');
	//can be 'fit' or 'field' to get just fields or to fit the data
	std::string toDo=pset.getParamS("toDo", "params", false, "fit");
	if(toDo!="fit" && CoilList.size()>0){
		myFields.resize(CoilList.size());
		for(int i=0;i<CoilList.size();++i){
			pset.addSection(CoilList[i]);
			myFields[i].read(pset,CoilList[i], false); //reread the coils
			curB=&myFields[0];
		}
	}else if(CoilList.size()>0){
		myFields.resize(1);
		myFitter.read(pset.section("fitter"), pset.section(CoilList[0]));
		logFile<<myFitter<<endl;
		myFields[0].read(myFitter.coilParams(),"", false); //reread the coils
		myFields[0].grid=myGrids[0]; //set the grid
		curB=&myFields[0];
	}


/*** the fitting input parameters ***/
//want Bz to give 1 MHz field
	ztarget=pset.getParamD("ztarget", "fitter", false, 1.0e5);

//want z-x and z-y to form the amgic angle
	yztarget=pset.getParamD("yztarget", "fitter", false,acos(1/sqrt(3.0))*180.0/PI)*PI/180.0;
	xztarget=pset.getParamD("xztarget", "fitter", false,acos(1/sqrt(3.0))*180.0/PI)*PI/180.0;

	if(toDo=="fit")
	{
		int err;
		MNINIT(5,6,7);

	//set the fitting parameters
		for(int i=0;i<myFitter.size();++i){
			char *moo=const_cast<char *>(myFitter.minuitString(i).c_str());
			MNPARS(moo, err);
		}

	/*** The MINUIT run command ***/
		int MAXLINE=256;
		char *command; command=new char[MAXLINE];
		snprintf (command, MAXLINE, "MIGRAD");
		MNCOMD (minuitfcn, command, err, NULL);
		delete command;
	}

/*** dump out the good data ***/

//get the file names
	std::string matBase=pset.getParamS("matout", "params", false, "field");
	if(matBase.find(".mat")<matBase.size())
		matBase=matBase.substr(0, matBase.size()-4);
	std::string textBase=pset.getParamS("textout", "params", false, "field");
	if(textBase.find(".biot")<textBase.size())
		textBase=textBase.substr(0, textBase.size()-5);

//just dump the data
	if(toDo!="fit"){
	//we dump out a histogram of B1 v B0 if we have 2 coils and the B1 and B0 flags
		std::string B1coil=pset.getParamS("B1", "params", false, "");
		std::string B0coil=pset.getParamS("B0", "params", false, "");
		int totSize=0, curOff=0;
		for(int i=0;i<myGrids.size();++i) totSize+=myGrids[i].size();
		Vector<coord<> > B1(totSize,0.0), B0(totSize,0.0);

	//pointers to the B0 and B1 coils data sets
		MField *B0field=NULL;
		MField *B1field=NULL;

	//The fields at each Grid...(so we only need to calcualate them once)
		Vector<Vector<coord<> > > B0Fields(myGrids.size());
		Vector<Vector<coord<> > > B1Fields(myGrids.size());

		for(int j=0;j<myFields.size();++j){
			curB=&myFields[j];
			curOff=0;
			for(int i=0;i<myGrids.size();++i){
				curB->grid=myGrids[i];
				curB->calculateField();

			//get the histogram profiles for B1
				if(B1coil==CoilList[j]){
					B1field=curB;
					B1Fields[i]=curB->Bfield;
					double tmpSum=0.0;
					for(int k=curOff, l=0;k<(myGrids[i].size()+curOff);++k, ++l)
					{	B1[k]=curB->Bfield[l];	}
				}
			//get the histogram profiles for B0
				if(B0coil==CoilList[j]){
					B0field=curB;
					B0Fields[i]=curB->Bfield;
					double tmpSum=0.0;
					for(int k=curOff, l=0;k<(myGrids[i].size()+curOff);++k, ++l)
					{	B0[k]=curB->Bfield[l];		}
				}

				curOff+=myGrids[i].size();

				if(MPIworld.master()){
					curB->writeMatlab(matBase+itost(j)+"_"+itost(i)+".mat");
					curB->write(textBase+itost(j)+"_"+itost(i)+".biot");
				}
			}
		}
		if(B0coil != "" && B1coil!=""){
			std::string homoProf=pset.getParamS("histogram", "params", false, "fieldHist.mat");
			matstream profs(homoProf);
			profs.put("B0",B0);
			profs.put("B1",B1);
		}

		cout<<"time for Field Calc:"<<stopwatch()<<endl;
		stopwatch.reset();

	//Simulate a spin system wth offsets given by each of the
	// B0 grid points
		if(toDo=="sim"){
			if(B0coil=="" || !B0field){
				std::cerr<<std::endl<<"Error: to perform spin sim, you must tell "<<std::endl;
				std::cerr<<" which coil is the 'B0' field "<<std::endl;
				exit(1);
			}
			int npts=pset.getParamI("npts", "params", false, 1024);
			double sw=pset.getParamD("sw", "params", false, 20000.0);
			Parameters spinSec(pset.section("params"));
			spinSec.addSection("spins");
			SolidSys sys(spinSec.section("spins"));
			if(MPIworld.master()) cout<<sys<<endl;

			HamiltonianGen myHams;
			matrix roeq=myHams.Hamiltonian(sys,pset.getParamS("roeq", "params", false, "Ix"));
			matrix detect=myHams.Hamiltonian(sys,pset.getParamS("detect", "params", false, "Ip"));
			exsituSim myfid(sys, npts, sw);
			myfid.B0fields=B0Fields;
			myfid.B1fields=B1Fields;
			Vector<complex> fid=myfid.FID(roeq, detect);
			if(MPIworld.master()) plotterFID(fid, "fid", 1.0/sw);
			cout<<"time for FID Calc:"<<stopwatch()<<endl;
		}

//if we just did a fit
	}else{
		if(MPIworld.master()){
			curB->writeMatlab(matBase+itost(myGrids.size()-1)+".mat");
			curB->write(textBase+itost(myGrids.size()-1)+".biot");
		}
		for(int i=0;i<myGrids.size()-1;++i){
			curB->grid=myGrids[i];
			curB->calculateField();
			if(MPIworld.master()){
				curB->writeMatlab(matBase+itost(i)+".mat");
				curB->write(textBase+itost(i)+".biot");
			}
		}

	}

/*** Stop MPI ***/
	MPIworld.end();
}



