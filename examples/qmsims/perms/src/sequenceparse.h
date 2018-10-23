

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-25-02
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
 	sequenceparse.h --> this takes in a 'vector string' (usually from
 	a Parameter set, and will create a 'PulseData' vector  from the input


 	to set a Pulse the function 'pulse' is used and can be

 	pulse(time) used if amplitude is set already via the 'amplitude' function
 	pulse(time, phase)
 	pulse(time, phase, amplitude)
 	pulse(time, phase, amplitude, offset)

 	the amplitude, is in Hz, phase in degreed, and offset in Hz

 	delays are set via

 	delay(time) --> NONE incremented delay
 	delay(time, dt) --> incremented delay (a 2D) type step

 	rotor anlges are set via

 	rotor(angle) --> angles in degrees

 	cycler traces are set via

 	cycler(SpinOp)

	To change the spin system use

	spinsys(sysname)

	To change the powder system use

	powder(powsec)
*/

#ifndef __sequenceparse_h__
#define __sequenceparse_h__ 1

#include "blochlib.h"
#include "pulsedata.h"

using namespace BlochLib;
using namespace std;



/***************** The main Sequence class **********/
class SequenceParse :
	public ScriptParse
{
	private:


	//creates a valid line from for the PulseData to read
	//proper line to the sequence vector
		void addLine(std::string inSqe);
		bool addGvars;

	//adds some global vars to the parser list
		void addGs();

	//adds things to our pulse list
		void addPulses(std::string curLin);

	//several 'continous' variables
		double curRot, curAmp, curOff;
		std::string curSys, curPow, curOn, curCycler,
				rWhite, inside;
		double currentTime;

	public:
		Vector<std::string> sequence; //the 'output' vector
		Vector<PulseData> pulses;

		SequenceParse();
		SequenceParse(bool gv);

		SequenceParse(Parameters &pset, bool addG=false);
		SequenceParse(const Vector<std::string> &pset, bool addG=false);

		SequenceParse(const SequenceParse &pset);
		void operator=(const SequenceParse &pset);

	//this is the master decider function
	// and snaps it to the proper above function
		virtual bool decide(std::string inSqe);

	//this is the master decider function
	// and snaps it to the proper above function
		bool decideSeP(std::string inSqe);
};


#endif
