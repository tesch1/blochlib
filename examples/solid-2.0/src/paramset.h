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
	paramset.h -->

	maintains a list of the normal NMR parameters
	it use the 'ParserGlobalVars' as its container for number types
	(so that they can be updated byt any assignments later in
	the sequence parser)

	strings are simply maintained statically...

	to get bools simply use 'getB("name")'
	to get ints simply use 'getI("name")'
	to get doubles simply use 'getD("name")'
	to get strings simply use 'getS("name")'
*/

#ifndef __paramset_h__
#define __paramset_h__ 1

#include "blochlib.h"

using namespace BlochLib;
using namespace std;



class ParamSet
{
	private:

	//this is the list of number params
	//that the constructor assignes to the ParserGlobalVars

		static std::string doubleList[10];
		static double doubleListdef[10];

		static std::map<std::string, std::string> StringList;

	//these are the strign params
		static std::string stringList[4];
		static std::string stringListdef[4];
		void init();
		void init(const Vector<std::string> &pset);

		Parser myParse;

	public:
		ParamSet();
		ParamSet(Parameters &pset);
		ParamSet(const Vector<std::string> &pset);

		void parse(Parameters &pset);
		void parse(const Vector<std::string> &pset);

		bool getB(std::string name);
		int getI(std::string name);
		double getD(std::string name);
		std::string getS(std::string name);

		void  set(std::string name, bool in);
		void  set(std::string name, int in);
		void  set(std::string name, double in);
		void  set(std::string name, std::string in);

		friend std::ostream &operator<<(std::ostream &oo,const ParamSet &out);
};





#endif
