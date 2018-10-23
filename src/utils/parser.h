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
  parser.h --> this class takes in a string, parses it and
  then evalutates the expression

  it has the option to both Add functions AND to add variables
  to the master maps...

  the initially there are only 2 variables set..Pi and E

  there are a variety of class that comprise
  a parsed object, they are not teribly iimportant
  (they do the stack analysis and iterative evaluation)

  There are 2 types of variable declarations

  the Global one, which ALL parser objects
  will be able to see, and the internal parser variabls
  which only that object will be able to see.

  global variables must be declared in a 'main' or another function
  before they are used...this container is GLOBAL to everything

  internal variables must be declared using the objects insertVar function

  this class is also templated to allow
  complex, or double (or smaller type) to be used as the output

  if you are using the complex version, the 'complex(#,#)' function
  makes a complex number for you in the input expression

 */

#ifndef __Parser_Data_h__
#define __Parser_Data_h__ 1

#include "container/Vector/Vector.h"
#include "utils/utils.h"

#include <stdio.h>
#include <string>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <map>

BEGIN_BL_NAMESPACE


/***** DATA ELEMENT CLASSES *****/
//these are the Stack containers
//basic data tructures
template<class T>
class  ParserStackElement
{
	friend class Parser;
	int    ElNo;		//the element number in the stack
	std::string OpType;		//the type of object (number, variable, function, symbol)
	struct {
		int LetNum;	//an array index
		std::string VarName; //variable name
		T _Value;	//a numeric value
		std::string MFun;	//a function name
	  } OpSpecifier;
};


/********************************************************/
/*			THe Variable Maps and Typedefs 			    */
/********************************************************/

typedef std::map<std::string, double > VarMap;
typedef std::map<std::string, double >::iterator VarMapIter;
typedef std::pair<std::string, double> aVar;

extern VarMap ParserGlobalVars; //the GLOABL variables object

/******************************/
/****** Global Vars *******/
/******************************/
#define N_pi  3.141592653589793238
#define N_e   2.718281828459045235


/**** THe MAster Parser Class ****/
class Parser
{
	private:
		//the error found global var
		static char SYMBOL[15];
		static char MATHSYMBOL[13];

		Vector<double> NStack; //the working stack evalutations
		int tos; //the stack posistion

		//some numerical limits
		static double lnlim;      // Smallest number allowed
		static double loglim;    // Smallest number allowed in call to log10() *

		bool ErrorFound;		//error flag

		int n;	//holds the current parser string posiston
		std::string chunk; //the piece in the parser we are on.
		std::string SymbGroup; //the current 'type'

		Vector<ParserStackElement<double> > usingStack; //the master stack for this clas

		VarMap internalVars; //the internal variables

		std::string InExpr; //the input expression

	//Variable testing
		//is int input a function
		int isMatFunc(std::string matf);
		//is the input a variable
		int isVariable(std::string variable);
		//is the input a 'math' operation
		int isSymbol(char ch);
		//is the input a 'math' operation
		int isMathSymbol(char ch);

	//the master parseing function set
	//gose from 1..9  ina a recusive way
		void GetNext();

		//gets the first chunk and starts the parser
		void First();

		//get a + or - (operations)
		void Third();

		//gets '*' or a '/'
		void Fourth();

		//gets '^'
		void Fifth();

		//gets the + or - (sign changes)
		void Sixth();

		//mathches '(' and  ')'
		void Seventh();

		//gets numbers or variables (x,y,Pi, or e)
		void Eighth();

		//gets the math functions
		void Ninth();

	//the computation parts
		//add the input to the stack
		void Push(double j);
		//get the current input, and drease the stack size
		void Pop(double &j);
		//updates the 'Nstack' with the computed value of
		// MF, removing the 'old' values from the stack
		bool Compute(std::string MF);

	public:

		Parser();
		Parser(std::string inp);
		Parser(const Parser &in);

		void operator=(const Parser &in);

	//the master parsing function...
		bool parse(std::string expr);

	//the master parsing function...
		bool parse(const char *expr);

	//add an interanl variable
		void addVar(std::string name, double val);
		void updateVar(std::string name, double val);
		double getVar(std::string name);
		static double getGlobalVar(std::string name);
		static	void addGlobalVar(std::string name, double val);
	//returns the value for the proper input variables
	//the vars 'x', 'y', 'z', and 't'
	// are internal vars iniitally set to 0
		double evaluate();
		double evaluate(double x, double y=0.0, double z=0.0, double t=0.0);

		inline double operator()()
		{	return evaluate();	}

		inline double operator()( double x, double y=0.0, double z=0.0, double t=0.0)
		{	return evaluate(x,y,z,t);	}

		//will be false if a failure happen
		inline operator bool()	{	return ErrorFound;	}
		inline bool fail()	{	return ErrorFound;	}
};


END_BL_NAMESPACE


#endif

