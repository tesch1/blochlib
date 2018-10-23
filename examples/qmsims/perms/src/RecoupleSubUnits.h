

#ifndef _RecoupleSubUnits_h_
#define _RecoupleSubUnits_h_

//#include "fullbinarytree.h"
#include <string>
#include "blochlib.h"
#include "pulsedata.h"
#include "spindex.h"
#include "sequenceparse.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;




/*
RecoupleSubUNit reads in two strings or parameters
	1) mode--> "CpCpp", or "180" for either permutation C or baring C mode
	2) subunit --> a string composed of characters like Types "1" or "2"...
	  such that the small "1" represents a 'normal' sequence
	  and a '2" represents a barred sequence...so this subunit is
	  acctually N strings (a Vector<std::string>)
	  where the first string is the first subunit, and the secoond string is the
	  next subunit...

	this class generats the two 'RecoupleTrains' in acordance to these subunits
	and are subsequently used in the 'RecoupleContent to generate te permuataion list

*/

class RecoupleSubUnits
{
	private:

/*** parse and choppig the 'Units' string ***/
	//private function that parses the input
	// like "oowwee, ooOOEE, " into 'oowwee' 'ooOOEE'
		Vector<std::string> parseUnits(std::string in);

	//private function that parses the input
	// like "oowwee" intop 'o' 'o' 'w' 'w' 'e' 'e'
		Vector<char> chopUnits(std::string in);

/*** sub units ***/
	//holds the valid chars that can be in the 'parseUnits' string
	// like 'o', 'O','w','W' etc
		Vector<char> subs_;

	//will take in a string like "o, O, w, W" (char units) and
	//parses them into individual pieces
		Vector<char> parseSubUnits(std::string);

/*** types ***/
	//holds 'what' type each char in 'subs' is supposed to be
	// these should be defined in a Parameter section
		Vector<std::string> types_;

	//checks that the types are allowed (R, Rp, Cp...)
		bool checkType(std::string);

	//will take in a string like "R, Rp, Cpp" (type units) and
	//parses them into individual pieces
		Vector<std::string> parseTypes(std::string);

//		RecoupleTrain generateTrains(Vector<char> in, std::string intype, int baseSym, int factorSym, double baseamp);

	//calculates the number of Dt time stesp for a specific
	//value of Dt
		int calcPropSteps(double maxtstep, double bt, double et);

	public:


/*** Constructors ***/
		RecoupleSubUnits();
		RecoupleSubUnits(std::string intype, std::string inunits, std::string insubs);
		RecoupleSubUnits(Vector<std::string> intype,  Vector<std::string> inunits,  Vector<std::string> insubs);

		void setParams(std::string insubs,std::string intype,  std::string inunits);
		void setParams(Vector<std::string> insubs,Vector<std::string> intype,  Vector<std::string> inunits);

/*** The Units and their propogators ***/
	//the 'Units' list
	//like "ooOO" "OOoo", etc
		Vector<std::string> units;

	//Since we assume that all each char in 'subs' is
	// a periodic propogator, we need to calculate them ONCE
	// for each powder angle...this is the list of each char propogator
	// and its generation function
		Vector<matrix> subProps;

	//the container for the sequneces....
	// each character in 'subUnits'
		Vector<Vector<PulseData> > trains;

	//this is the sequence time for each 'subunit'
		Vector<double> subPropTime;

	//grabs and generates all the pulse trains
	// it returns the rotor speed which should be defined
	// in the parameters section
		double generateTrains(Parameters &in);

	//you need to sys.setPowderAngles before you call this
		void generateSubProps(
				SolidSys &sys,	// sys--> the Hamiltonian generator
	 	        double wr,
				double dtmax	// dtmax --> maximal time step allowed
			);

	//typically the 'units' string is going to be a number of 'subs'
	// and these too will be periodic and thus only need to be calculated
	// once for each powder angle and we can use the 'subProps' to calc
	// them as well SOO YOU MUST call 'generateSubProps' then 'generateProps'
	// in that order...
	// NOTE: if 'units' has a length=1 then is simply uses the subProps
		Vector<matrix> Props;
		Vector<double> propTime; //time for each propogator

	//you need to sys.setPowderAngles before you call this
	// this function will call 'generateSubProps'
		void generateProps(
				SolidSys &sys,	// sys--> the Hamiltonian generator
	 	        double wr,
				double dtmax	// dtmax --> maximal time step allowed
			);

	//you need to sys.setPowderAngles before you call this
	// and call generateSubProps
	// this function will NOT call 'generateSubProps'
		void generateProps();

/*** Auxilliary THings  ***/

	//Debug Flag
	//if 0 no debug
	//if 1 print out the "Text" representation as to what is going on
	//if 2 print out even more text info...
		int debugflag;

	//prop size...the number of 'units'
		int size() const {	return units.size();	}
		int propSize() const {	return Props.size();	}
		int subPropSize() const { return subProps.size();	}

};




#endif
