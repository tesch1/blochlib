

#include "blochlib.h"

//This extends the SolidSys class
// and overwrites the two
// Hamiltonian Functions required by
// the \'oneFID\' object..
//it simply returns the matrix generated from
// the input hamiltonian
using namespace BlochLib;
using namespace std;

class HamilSys: public SolidSys
{
    public:
        std::string hamil;
        HamiltonianGen myGen;

        HamilSys(int nspins, std::string Hamil=""):
            SolidSys(nspins),
            hamil(Hamil),
            myGen(Hamil)
        {}

        hmatrix Hamiltonian(double t1, double t2,double wr=0.0)
        {
            SpinSys & tm=(SpinSys &)(*this); //cast a ptr back to SPinSys
            return myGen.Hamiltonian(tm, hamil, theRotations.theta, theRotations.phi);
        }

        hmatrix Hamiltonian(double wr, double rot, double alpha, double beta, double t1, double t2)
        {
            SpinSys & tm=(SpinSys &)(*this); //cast a ptr back to SPinSys
            return myGen.Hamiltonian(tm, hamil, alpha, beta);
        }
};


int main(int argc,char* argv[]){

	//start MPI
    MPIworld.start(argc, argv);

    std::string inf="";
    if(MPIworld.master()) query_parameter(argc,argv,1, "InputFile:: ", inf);

    MPIworld.scatter(inf);

  //make a parameter set
    Parameters pset(inf);
  //decalare our \'HamilSys\'
    HamilSys mysys(
        pset.getParamI("numspins"),
        pset.getParamS("hamil")
        );

  //nnum points in fid
    int npts=pset.getParamI("npts");
  //the sweep width
    double sw=pset.getParamD("sw");
  //our detection matrix
    std::string detectST=pset.getParamS("detect");
  //our starting matrix
    std::string roeqST=pset.getParamS("roeq");
  //powder average type
    std::string aveType=pset.getParamS("aveType");
    int thetaStep=pset.getParamI("thetaSteps", "", false);
    int phiStep=pset.getParamI("phiSteps", "", false);

  //decalre our powder
    powder zcw(aveType, thetaStep, phiStep);

  //using this hamiltonianGen to create the
  //detect and roeq matrices
    HamiltonianGen mygen;
    matrix detect=mygen.Hamiltonian(mysys, detectST);
    matrix roeq=mygen.Hamiltonian(mysys, roeqST);

  //our fid vector
    Vector<complex> fid(npts, 0);
  //decalre an \'oneFID\' object
    oneFID<HamilSys> myfid(mysys, npts, sw);

  //set the MPI controller of the oneFID
    myfid.Controller=MPIworld;

  //collect the FID
    fid=myfid.FID(zcw, roeq, detect);
    complex sm=0;
  //dc offset correct it
    sm=sum(fid(Range(npts/2, npts)))/double(npts/2.0);
    fid-=sm;
  //print out the data
    if(MPIworld.master()) plotterFID(fid, "fid", 1.0/myfid.sw());

  //end MPI
    MPIworld.end();
}

