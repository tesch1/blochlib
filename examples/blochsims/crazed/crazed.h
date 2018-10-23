
#ifndef _crazed_h__
#define _crazed__h__

#include "blochlib.h"
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

extern timer stopwatch;
void printTime(int nrounds=1);
void Info(std::string mess);


typedef XYZfull TheShape;
typedef XYZshape<TheShape> TheGridS;
typedef GradientGrid<TheGridS > TheGrid;
typedef ListBlochParams< TheGrid > MyPars;

//Extra ineractions
typedef TanhScale Scaler;
typedef Interactions<Offset<TheGrid>, Relax<>, BulkSus, RadDamp, DemagField<TheGrid, Scaler> > MyInteractions;

//typedefs for Bloch parameter sets
typedef Bloch< MyPars, Pulse, MyInteractions> PulseBloch;
typedef Bloch< MyPars, NoPulse, MyInteractions> NoPulseBloch;

class SimHolder{
	private:
	public:
		parameters pset;
		TheGrid jj;
		MyPars mypars;
		MyInteractions MyInts;

		Offset<TheGrid> myOffs;
		Relax<> myRels;
		BulkSus BsRun;
		RadDamp RdRun;
		Scaler myscal;
		DemagField<TheGrid, Scaler> DipDip;

		PulseBloch myparspulse;
		NoPulseBloch me;

		TheGrid GenGrid();
		void GenParams();
		void GenInteracts();



		SimHolder(parameters &inp):
			pset(inp), jj(GenGrid())
		{
			GenParams();
			GenInteracts();
		}
};




#endif

