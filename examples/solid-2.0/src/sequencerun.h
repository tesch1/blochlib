

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-27-02
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
 	sequencerun.h -->
 	this maintains the pulse sections and performs the
 	propogation and fid collection (the main runner class)

 	it is public to sequence parse so it can also parse
 	'non-section sequences'

 	it adds  these function for scripting
 	(it also has the sequenceparse ones well)

 	use(seqname)
 	  --> use this subsection
 	use(seqname, repeat)
 	  --> use the subsection reapeating it n times
 	use(seqname, repeat, hold)
 	  --> use the subsection reapeating it n times
 	  --> and 'held' so that it is only appiled once in 2D cycles

 	fid()
 	  --> collect an fid at this poin in the prop run
 	  --> if the sequence is ptop it will collect one point

 	fid(index)
 	  --> collect an fid at this poin in the prop run
 	  --> if the sequence is ptop it will collect one point
 	  --> and ADDS it to an fid at index

 	savefid()
 	  --> saves an fid into a file...using the default out name

 	savefid(name)
 	  --> saves an fid into a file of name 'name'

 	ptop()
 	  --> spcifies that the sequence is point to point

	pulse, delay, rotor, cycler, (from the 'sequenceParse') are also available

	atlerSys(what, number)
	  --> changes a coupling or some internal param to the current spin system

NOTE:: any variables decalred outside the 'subsections'
in the pulses parameter set will be global to all the
subsections....

it also hold the main ParamSet

*/

#ifndef __sequencerun_h__
#define __sequencerun_h__ 1

#include "sequenceparse.h"
#include "paramset.h"
#include "multipowder.h"
#include "multisystem.h"
#include "propogation.h"
#include "pulsedata.h"

/*** This is a helper class that hold the subsections
It allows the 'use' command below to save several CPU clicks

if the sequence is ptop, then the sections need to be
calculated only once for each powder angle

it also maintains that sections propogator matrixs
and the PulseData list.  It will update the
Pulse data lists time when the main sequence runs
as its start time will probably change
***/
using namespace BlochLib;
using namespace std;



class subSequence{
	friend class SequenceRun;
private:


	//the pulse data Vector
		Vector<PulseData> myPulses;

	//the input string vector of the paraser
		Vector<std::string> data;

		bool haveParsed;
	public:

	//a ptr to the main propogator object
	// to be set by the SequenceRun class
		Propogation *props;

	//the hold and repeat flags
		bool hold;
		int repeat;

		int lastPow;
		double lastT;

		std::string cycler;
		bool docycler;

		subSequence();
		subSequence(const subSequence &cp);
		void operator=(const subSequence &cp);

		subSequence(const Vector<std::string> &indata, Propogation &inp);
		~subSequence();
	//performs the propogator..simply calculates the thing if
	// it needs to...

	//the propogator
		matrixs U;
		matrixs holdU;

		void propogator(int powpt, int pt=-1, double startT=0);
		void propogate(matrixs &ro, int powpt, int fidpt);

		Vector<PulseData> &parse();

};

std::ostream &operator<<(std::ostream &otr, subSequence oo);

/***************** The main Sequence class **********/
class SequenceRun:
	public SequenceParse
{

	private:

	//these are the current running params for the main sequence
		double amp, phase, offset, currentTime;
		std::string on, usesys, usepow, cycler;

	//this adds the global parameters to the
	// top of any 'sub' section so that the parameters
	// carry to them...if they
	// are set to different things in the  sections,
	// then they will be overrdden
	//	void addGlobalPars(Vector<std::string>

	//these advance the propogator 'U' for the MAIN sequence
	// not for the Sub Sequences...the 'subSequence' class handles those
		//void doPulse(std::string inSqe);
		//void doDelay(std::string inSqe);

		void doPulseData(std::string inSqe);

		void setAmplitude(std::string inSqe);
		void setOn(std::string inSqe);
		void setRotor(std::string inSqe);
		void setOffset(std::string inSqe);
		void setSystem(std::string inSqe);
		void setPowder(std::string inSqe);

		void doCycler(std::string inSqe);


	//this should take in a string like "ptop()"
	// sets the point to point flag
		void addPtop(std::string inSqe);

	//this should take in a string like "A=B"
	// and adds the var to the GLOBAL vars
		void addGVar(std::string inSqe);

	//Shows the current pulse sequence
		void Shows(std::string inSqe);

	//sets the sim as a 2D
		void to2D(std::string inSqe);

	//Creates the Sequence for a 'use' command
	// SHOWS it (it does not do it).
		void showUsed(std::string inSqe);

	//The FID dooer....
		void doFID(std::string inSqe);

	//save the fid at its current count....
		void saveFID(std::string inSqe);

	//this resets the FID(s) to be 0 everywhere
	// obeys the command //reset()
		void ResetPars(std::string inSqe);

	//sets a new detection input string...
	// the detection is set in the 'doFID' function
	// so all we need to do here is reset the input string
	// this could be done with a 'detect=SPinOp" command as well
	// but ust to make the sysntax similar to 'ro(SPinOp)'
		void setDetect(std::string inSqe);

	//sets ro to the roeq if nothing is in the argument,
	//is something is pressent is sets ro to that value...
		void setRo(std::string inSqe);

	//sets/changes a spin system coupling param
		void doAlterSys(std::string inSqe);

	//dump the propogator info are the time
		void showProps(std::string inSqe);
	
	//if somethign happens, like a [ctlr]-C or an error
	// setting this flag will dump the current state
	// of the system for looking at
		void setDumpOnDie(std::string inSqe);
	
	//this is the 'global' pulse program
	// (everything NOT in the subsections)
		Vector<std::string> mainProc;

	//these are the sub pulse procs
		typedef std::map<std::string, subSequence > SubPulseMap;
		typedef std::map<std::string, subSequence >::iterator SubPulseMapIter;
		typedef std::map<std::string, subSequence >::const_iterator SubPulseMapConstIter;
		SubPulseMap subProc;

	//for each 'used' subSequence' it gets put into the ptr
	//list for later propogation...
		Vector<subSequence> PropList;
		int numProp;

	//this holds the fids...for 1Ds
		Vector<scomplex> fid1D;

	//this holds the 2D fids
		matrixs fid2D;

	//the pointopoint flag
		bool isPtoP;

	//is the sequence 2D?
		bool is2D;

	//is the powder done so we can save?
		bool canSave;

	//the inital density matrixs
		matrixs roeq, ro;

	//the detection matrixs
		matrixs detect;

	//the current propogator
		matrixs U;

	//the Propogrator class
		Propogation myProp;

	//the current powder point index
		int powpt;
	//the current 2D point index
		int pt2D;
	//has the data been saved yet?
		bool notSaved;
	
	//the dump on die flag
		bool dumpOnDie;
	public:
		
	//the main parameter set...
		ParamSet Params;

	//the main spinSys object
		MultiSystem Systems;

	//the main Power list
		MultiPowder Powders;

		Vector<std::string> totalSequence; //the total puls sequence at any time

	//set this flag to 'true' if we want to perform
	//the sequence and to 'false' if we just wish to
	//show the parser 'working'
		bool toRun;

		SequenceRun();
		SequenceRun(Parameters &pset);
		SequenceRun(Vector<std::string> &pset);

		//Paraser the pulses ection into the 'global' items and the
		// subPulse sections
		void parsePulse(Parameters &pset);

		void setParams(Parameters &pset);
		void setParams(const Vector<std::string> &pset);

		void setSystems(Parameters &pset);
		void setSystems(const Vector<std::string> &pset);

		void setPowders(Parameters &pset);
		void setPowders(const Vector<std::string> &pset);


	//this is the master decider function
	// and snaps it to the proper above function
		bool decide(std::string inSqe);

		void printSeq(std::ostream &oo) ;
		friend std::ostream &operator<<(std::ostream &oo, SequenceRun &out);

	//the master runner....
		int run();
	
	//dumps the current state of the system...
		void dumpState();
};


#endif
