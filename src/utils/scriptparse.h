

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
 	scriptparser.h --> this takes in a 'vector string'
 	the input allows for variable setting, and looping over parameters,
 	and simple if, elseif, else structures...

 	this is designed to be a base for specific
 	'functional' classes....i.e. classes that will add various functions
 	to the script.  This just acts as a basic loop and if control parser

 	it holds the 'Parser' object that holds the variables and does
 	the evaluation

 	subclasses should overwrite the 'decide' function that determins
 	what to do with things other then 'loop' or 'if' types of structores


 	variables are simply declared via a simple syntax

 	A=B (a gets set from B) where B can be some expression

 	loops are performed via the

 	loop i=1:40

 	...do things

 	end

 	syntax where the i=1:40 means loop from 1 to 40

	and

	if(t==9)
	 ...
	elseif(moo=34)
	 ...
	else
	 ...
	end

	and a 'print' command

	print(thing)

	that prints that item to the screen

*/

#ifndef __scriptparse_h__
#define __scriptparse_h__ 1

#include "parser.h"
#include "utils/params.h"

BEGIN_BL_NAMESPACE



/***************** The main Sequence class **********/
class ScriptParse{

	private:

	//this function parse the 'loop' keyword
		int parseLoops(Vector<std::string> &pset, int start=0);

	//this function parse the 'if' keyword
		int parseIfs(Vector<std::string> &pset, int start=0);

	public:
		Vector<std::string> data; //the input data vector

		Parser myParse; //the parser object

		ScriptParse(){}
		ScriptParse(Parameters &pset);
		ScriptParse(const Vector<std::string> &pset);

		ScriptParse(const ScriptParse &pset);
		void operator=(const ScriptParse &pset);

		virtual ~ScriptParse(){}

		bool parse();
		bool parse(Parameters &pset);
		bool parse(const Vector<std::string> &pset);

	//this should take in a string like "A=B" and add the
	//proper line to the sequence vector
		void addVar(std::string inSqe);

	//this should take in a string like "A=B" and add the
	//proper line to the sequence vector Asd it to the GLOBAL list
		void addGlobalVar(std::string inSqe);

	//prints the argument
		void doPrint(std::string inSqe);

	//this is the master decider function
	// and snaps it to the proper above function
	// this is designed to be redefined by a subclass
		virtual bool decide(std::string inSqe);
};

END_BL_NAMESPACE



#endif
