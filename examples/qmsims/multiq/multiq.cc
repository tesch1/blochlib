#include "blochlib.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;


/* Performs a simulation of the
 MultiQuantum Excitation scheme using a
 phase modified postC7 to pump and convert the
 multi-quantum excitation...

 to excite 'N' orders
 phi=0, 2Pi/2N, 4Pi/N...2Pi-->i*Pi/N for the i'th point in the FID

|excitation | conversion | point
 i*(pC7)       i*(pC7)pi/2    FID

 for i=0;

 | (pC7)_0 | (pC7)_pi/2 |

 for i=1

 | (pC7)_0(pC7)_Pi/N |(pC7)_pi/2 (pC7)_pi/2 |

 etc...

 in order for this simulation to immitate
 reality, MANY spins must be simulated, thus 'brute-force'
 calculation of the pC7 sequences would be too lengthy

 so we use both the periodic property of the C7 and
 the fact that a global phase shift to the sequence is
 simply a Z rotation of the original C7...thus we
 need to calculate the first pC7, and simply phase shift
 the rest

 */

int main(int argc,char* argv[]){

/*** grab the config file ***/
	std::string fname;
	int q=1;
	query_parameter(argc,argv,q++, "Enter File To Parse: ", fname);
	Parameters pset(fname);

/*** add the usable sections ***/
	pset.addSection("spins");
	pset.addSection("parameters");

/*** Grid Set up ***/
	coord<> mins=pset.getParamCoordD("min","spins");
	coord<> maxs=pset.getParamCoordD("max","spins");
	coord<int> dims=pset.getParamCoordI("dim","spins");
	Grid<UniformGrid> BasicG(mins,maxs,dims);
	XYZshape<> MasterGrid(BasicG, XYZfull());

/*** spin system set up ***/
	SolidSys mysys(BasicG.size()); //have as many spins as grid points

	//now we must add all the dipole-dipole couplings
	double maxD=pset.getParamD("maxD","spins");
	std::string neighbor=pset.getParamS("neighbor","spins", false, "all");
	double Vol, difflen;

	XYZshape<>::iterator Grid1(MasterGrid);
	while(Grid1){
		XYZshape<>::iterator Grid2(MasterGrid);
		while(Grid2){
			if(Grid1.curpos()!=Grid2.curpos() && Grid2.curpos()>Grid1.curpos()){
				Vol=cube(min(Grid2.dr()));
				difflen=cube(norm(Grid2.Point()-Grid1.Point()));
				if(neighbor=="nearest" && (abs((Vol/difflen)-1.0)) < 1.0e-6){
					mysys.dip.push_back(Dip(maxD*Vol/difflen, Grid1.curpos(),Grid2.curpos()));
				}else if(neighbor!="nearest"){
					mysys.dip.push_back(Dip(maxD*Vol/difflen, Grid1.curpos(),Grid2.curpos()));
				}

			}
			++Grid2;
		}
		++Grid1;
	}
cout<<mysys<<endl;
	mysys.setSpinMats(); //set all the basic spin matricies
	mysys.setCrystalAs(); //set the cyrstal orientation axis
	mysys.setMats(mysys); //set the 'Zero' and 'One' matrix

/*** get pulse sequence parameters ***/
	int N=pset.getParamI("N","parameters");
	double cycleTime=pset.getParamD("cycleTime","parameters");
	double wr=1/cycleTime;
	double maxtstep=pset.getParamD("maxtstep","parameters");

	std::string initS=pset.getParamS("initro", "parameters");

	HamiltonianGen hamgen;
	matrix roeq=hamgen.Hamiltonian(mysys, initS);


}
