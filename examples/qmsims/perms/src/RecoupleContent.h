
#ifndef _RecoupleContent_h_
#define _RecoupleContent_h_

//#include "fullbinarytree.h"
#include <string>
#include "blochlib.h"
//#include "RecoupleSequence.h"
#include "RecoupleSubUnits.h"
#include "spindex.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;



//this class is the work horse for determination of
// tensorial/spin tensor content of a given recoupling sequence
class RecoupleContent{

	private:
		static const int Save=0x0001;
		static const int Run=0x0002;
		bool killFile;
	public:

		std::string outFileName; //the output file name if "GenerateAndQuit"
		std::fstream *outFile; //the out file handle


		SolidSys sys; //the spin system
		powder pows; //the powder angles

		Vector<Vector<int> > trains; //generated lists of trains

		Vector<TensorGen> tensors; //the lists of our tensors

		double maxtstep; //the max time step for propogation
		double wr; //spinning speed
		double rotorang; //rotor angle
		std::string pulseON; //spin to pulse

		matrix traces; //holds the traces for each tensor and train

		RecoupleSubUnits SubUnits; //recoupling subunits

		int progress;

/*** Permutation Calculating ***/


	//this calulates the 'length' of the permutation
	// RecoupleSubUnits generate 'n' Propogators, but we
	// typically only what N combinations of them
	// so the value returned by this function is almost always
	// N, or the permutation length
		int permutationLength(int order);

	//generates the Next K Subset of a generic list of
	// N numbers (0, 1,2,3...) where we only want
	// a sub set of n numbers...of course there are many posibilities
	//if 'more'  is false, the list initializes, if it is true it
	// calcs the next subset
		void nextKSubset(int n, int k, int *subset,bool &more);

	//this fills up the 'index' list out to the appropriate length
	// as given by 'permutationLength(order)'  it simply
	// returns a list with int Mod SubUnits.size()
		void generateLists(int permLen, int *propList);

	//master function that generates allthe valid reocuple pulse trains
	// fill up the 'trains' vector with indexes that are to match
	// the 'props' in 'SubUnits'
		void generateTrains(int order);
	//generates trains and dumps them to a file...
		void generateTrains(int order, std::fstream &outF);
	//generates trains and dumps them to a file...
		void generateTrains(int order, std::string outF);
	//reads trains from a file...
		void generateTrains(std::fstream &outF);
	//reads trains from a file...
		void generateTrains(std::string outF);

	//generates the 'master' trains and the SubUnit trains
		void generateTrains(int order,Parameters &in);
	//generates trains and dumps them to a file...
		void generateTrains(int order,Parameters &in, std::string outF);
	//generates trains and dumps them to a file...
		void generateTrains(int order,Parameters &in, std::fstream &outF);

/*** Constructors ***/
		RecoupleContent():
			killFile(false), outFileName(""), outFile(NULL),
			maxtstep(1.0e-6)
		{}

		RecoupleContent(std::string outname):
			killFile(false),outFileName(outname), outFile(NULL),
			maxtstep(1.0e-6)
		{}

		RecoupleContent(int len, double mindt=1.0e-6):
			killFile(false), outFileName(""), outFile(NULL),
			maxtstep(mindt),progress(false)
		{	generateTrains(len); 	}//generateTrains(len, type, baseSym, factorSym, baseamp); 	}

		RecoupleContent(int len, double wrin, double rotor, std::string pulse, double mindt=1.0e-6):
			killFile(false),outFileName(""), outFile(NULL),
			maxtstep(mindt),wr(wrin), rotorang(rotor), pulseON(pulse),progress(false)
		{	generateTrains(len); }//generateTrains(len, type, baseSym, factorSym, baseamp); 	}


		~RecoupleContent()
		{	if(killFile){ outFile->close(); delete outFile;	}	}

/*** The work horse Propogator and Trace calculators ***/

	//generate a propogator from a train 'which' and dpowder angle
		matrix generateProp(int which, powder::iterator &powit);

	//generate an FID from the train 'which'
	// the master runner for MPI master
		Vector<complex> FID(int which, matrix &roeq, matrix &detect,int npts=256, int napps=1);

	//generate an FID from the train 'which'
	// the slave runner for MPI slave
		Vector<complex> FID(int which, matrix &roeq, matrix &detect,int npts, int napps, int powpt);

	//generate a 2Q transfer profile
	// the master runner for MPI master
		complex transfer(matrix &roeq, matrix &detect, int &napps);

	//generate an 2Q from the train 'which'
	// the slave runner for MPI slave
		complex transfer(matrix &roeq, matrix &detect, int &napps, int powpt);

	//calculate the tensorial components for the train 'which'
	//napps is the number of sucsessive applications of the train (i.e. as in a 2D experiment)
		Vector<complex> calcTrace(int which, int napps);
		void calcTrace(int napps);
		void calcTrace(){	calcTrace(1);	}

	//total number of tensors inside the tensors vector
		int tensorSize() const;

	//calculates the sequence time for the trina 'which'
		double sequenceTime(int which);


	//gets the name of the squence 'which'as generated
	// from the units of 'SubUnits'
		std::string sequenceName(int which);

	// writes ONE train to a file (binary)
		void writeOne(std::ostream &out);

	// reads ONE train to a file (binary)
		void readOne(std::istream &out);

	// reads in the Trains to a file (binary)
		void read(std::string out);
		void read(std::fstream &out);

	// writes out the Trains to a file (binary)
		void write(std::string out);
		void write(std::fstream &out);

};

std::fstream &operator<<(std::fstream &out, RecoupleContent &oo);
std::fstream &operator>>(std::fstream &out, RecoupleContent &oo);


#endif
