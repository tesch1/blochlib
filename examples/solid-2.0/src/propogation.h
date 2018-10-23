
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
	Propogation.h -->does the propogation on a list of
	Pulse data's...it does not care how they are generated...
	(that is the job for the SequenceParser and other classes)

	important things that need to be declared before using it

	A VEctor of Pulse Datas

	A Map of SolidSYstems sets

	A Map of Powder Angles sets

	A ParamSet class that acts as a bridge between
	the ParerGlobalVars and the system
	(things like wr, maxtstep, etc are defined in that class)

*/

#ifndef __Propogation_h__
#define __Propogation_h__ 1

#include "pulsedata.h"
#include "paramset.h"
#include "blochlib.h"
#include "multipowder.h"
#include "multisystem.h"
using namespace BlochLib;
using namespace std;



class Propogation
{
	public:
		Vector<PulseData> pdata;

		MultiSystem *Systems;
		MultiPowder *Powders;
		ParamSet *Params; //ptr to paramset object

		powder *curPow;
		SolidSys *curSys;

		bool powChanged;

		Propogation();

		~Propogation();


		void update2Dtime(int pt);

		void setPowder(std::string use);

		void setSystem(std::string use);

		bool is2D();

		matrixs propogate(matrixs &ro, Vector<nonUniqueMapEle<std::string,int> > &powpt,  int pt=-1);

		matrixs propogate(matrixs &ro, int powpt,  int pt=-1);

		matrixs propogator(Vector<nonUniqueMapEle<std::string,int> > &powpt,  int pt=-1);

		matrixs propogator(int powpt,  int pt=-1);

};


std::ostream &operator<<(std::ostream &oo,const Propogation &out);


#endif

