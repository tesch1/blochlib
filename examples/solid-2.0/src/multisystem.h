/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-26-02
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
	multisystem.h -->
	maintains a list to a bunch of SolidSystems

	the one requirement is that the System does NOT change
	is Hilbert Space dimensions...the things that should change are
	the specific interactions...

	it parses a parameter set like
	-----

	#these are the global (or intial parameters)
	# the numspins and the type cannot change in the sub sections
	# if no subsections are present, the the global is the default
	# spin system returned
	numspins 2
	T 1H 0
	T 1H 1
	D 2000 0 1

	# the first spin section (will overwrite any global interactions)
	spin1{
		C 5000 -2556 0 1
		C -5000 -2556 0 0
		#D 2134 0 1
	}

	# the second spin section (will overwrite any global interactions)
	spin2{
		D 2534 0 1
		J 500 0 1
	}



*/


#ifndef __multisystem_h__
#define __multisystem_h__ 1

#include "blochlib.h"
#include <map>

using namespace BlochLib;
using namespace std;


class MultiSystem
{
	private:
		typedef std::map<std::string,SolidSys> SysMap;
		typedef std::map<std::string,SolidSys>::iterator SysMapIter;

		SysMap Systems; //contains the powder object by name

		Vector<std::string> global; //the global part of the system

	//cleans the input parameter set for easy and correct reading
		Vector<std::string> preParse(const Vector<std::string> &pset);

	public:

		MultiSystem(){}

		MultiSystem(Parameters &pset);

	//this will grab all the
	//subsections named 'powder1...powderN' and add them
	// to the std::map
		void parse(Parameters &pset);

	//add a powder to the list from a parameter set list
		void addSystem(std::string name,Parameters &PowderPset);

	//add a powder to the list from a SolidSys
		void addSystem(std::string name,const SolidSys &sys);

	//returns the pointer to the powder in the list
		SolidSys *getSystem(std::string name);

	//sets the Bfield in all the systems
		void setBfield(double BF);

		friend std::ostream &operator<<(std::ostream &oo,const MultiSystem &out);
};


#endif
