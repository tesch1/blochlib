

#include "blochlib.h"
#include "crazed.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

// PART 2 NEEDED BECUASE THE COMPILSE FARTS OUT WHEN COMPILE THE WHIO PROGRAM
// AS ONE FILE...sheesh...this make functions out make a grid


/* Homonuclear The CRAZED pulse sequeque */

/*

RF  ---90--T1---90----T2-FID
Grad --------Gzt--2TGz------

*/

timer stopwatch;
void printTime(int nrounds){
	 std::cout <<std::endl<< "Time taken: " << (stopwatch()/nrounds) << " seconds\n";
}


void Info(std::string mess)
{
	cout<<mess;
	cout.flush();
}




TheGrid SimHolder::GenGrid()
{
	coord<int> dims(pset.getParamCoordI("dim"));
	coord<> mins(pset.getParamCoordD("min"));
	coord<> maxs(pset.getParamCoordD("max"));


	Info("Creating grid....\n");
	Grid<UniformGrid> gg(mins, maxs, dims);

	Info("Creating inital shape....\n");
	TheShape tester;
	Info("Creating total shape-grid....\n");
	TheGridS grids( gg, tester);

//create the gradient grids..
	double zgrad=pset.getParamD("zgrad");

	Info("Creating Gradient map grids....\n");
	TheGrid out(grids);
	out.G(0.0,0.0,zgrad);
	return jj;
}

void SimHolder::GenParams(){

//Set up Parameter lists
	int nsp=jj.size();
	Info("Creating entire spin parameter list for "+itost(nsp)+" spins....\n");
	MyPars mypars(jj.size(), "1H", jj);
	nsp=mypars.size();

	double inBo=pset.getParamD("Bo");
	double inTemp=pset.getParamD("temperature");
	string spintype=pset.GetParamS("spintype");
	double moles=pset.getParamD("moles");


	Info("Setting spin parameter offsets....\n");
	for(int j=0;j<nsp;j++){
		mypars(j)=spintype;
		mypars(j).moles(moles);
		mypars(j).Bo(inBo);
		mypars(j).temperature(inTemp);
	}

	mypars.calcTotalMo();
	mypars.print(cout);
}



void SimHolder::GenInteracts()
{
	double D=pset.getParamD("bulksus");
	double t2s=pset.getParamD("T2");
	double t1s=pset.getParamD("T1");
	double offset=pset.getParamD("offset")*PI2;

	double tr=pset.getParamD("raddamp");
	double scalef=pset.getParamD("scalef");
//typedef Interactions<BulkSus, RadDamp> MyInteractions;
	Info("Setting Interactions....\n");

	myOffs=Offset<TheGrid>(jj, offset);
	myRels=Relax<>(mypars, (!t2s)?0.0:1.0/t2s, (!t1s)?0.0:1.0/t1s);
	BulkSus BsRun(D);
	RadDamp RdRun(tr);
	Scaler myscal(scalef);
	DipDip=DemagField<TheGrid, Scaler>(jj, myscal);
	MyInts=MyInteractions(myOffs, myRels,BsRun, RdRun, DipDip);

	if(scalef<=0) DipDip.Off();
}



