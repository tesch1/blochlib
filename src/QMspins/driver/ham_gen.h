/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 08.10.02
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
  ham_gen.h --> this class takes in a string, parses it and
  then evalutates the expression into a spinoperator matrix set

  it is a mas improvement on the old parser

  it alows for complex numbers (at least the 'complexi'
  variable to be used, to make numbers complex)

  'levels' have been allowed (you can use parenthesis ())

  complex variables can be assigned to the mix....

  the global variables pi, e, and complexi are automatically added

  this one will be a bit slower then the last one, however,
  it is much more powerful....

 */

#ifndef __HamiltonianGen_Data_h__
#define __HamiltonianGen_Data_h__ 1

#include "container/Vector/Vector.h"
#include "utils/utils.h"
#include "container/complex.h"
#include "container/matrix/matrix.h"
#include "QMspins/space_ten.h"
#include "QMspins/spin_ten.h"
#include "QMspins/spinsys.h"

#include <stdio.h>
#include <string>
#include <ctype.h>
#include <math.h>
#include <iostream>
#include <map>

BEGIN_BL_NAMESPACE


/***** DATA ELEMENT CLASSES *****/
//this is the basic spinop parsed object
class HamiltonianGenSpinOp{
	friend class HamiltonianGen;
	char type; //'I' or 'T'
	int rankl; //the l rank (if 'type=='T')
	int rankm; //the m rank (if 'type=='T')

	char cart; //the carteis type (x,y,z,p,m,e) (if type='I');

	int spin1; //the first spin
	int spin2; //the second spin (if Type='T');

	public:
		inline bool isOp(){	return (type=='T' || type=='I');	}
		HamiltonianGenSpinOp():
			type('N'), rankl(0), rankm(0), cart('N'), spin1(-1), spin2(-1)
		{}
		inline HamiltonianGenSpinOp operator=(const HamiltonianGenSpinOp &cp)
		{
			type=cp.type; rankl=cp.rankl; rankm=cp.rankm;
			cart=cp.cart; spin1=cp.spin1; spin2=cp.spin2;
			return *this;
		}

};

//this is the basic Space Op parsed object
class HamiltonianGenSpaceOp{
	friend class HamiltonianGen;
	int rankl; //the l rank
	int rankm; //the m rank

	public:
		HamiltonianGenSpaceOp():
			rankl(-10), rankm(-10)
		{}

		inline bool isOp(){	return (rankl!=-10 && rankm!=-10);	}

		inline HamiltonianGenSpaceOp operator=(const HamiltonianGenSpaceOp &cp)
		{		rankl=cp.rankl; rankm=cp.rankm; 	return *this; }
};

//these are the Stack containers
//basic data tructures
class  HamiltonianGenStackElement
{
	friend class HamiltonianGen;
	int    ElNo;		//the element number in the stack
	std::string OpType;		//the type of object (number, variable, function, symbol)
	struct {
		int LetNum;	//an array index
		std::string VarName; //variable name
		complex _Value;	//a numeric value
		std::string MFun;	//a function name
		HamiltonianGenSpinOp spinOp; //the parsed spin op function
	  	HamiltonianGenSpaceOp spaceOp; //the parsed space op function
	  } OpSpecifier;
};

/// This is the stack object
// it holds either a number or a matrix, and a tag that
// tells me which one to use...
class HamiltonianGenStack{
	friend class HamiltonianGen;

	char type; //either 'C' or 'M'
	complex number;
	matrix mat;
	public:
		HamiltonianGenStack():	type('C'){}
		inline void operator=(HamiltonianGenStack &rhs)
		{	mat=rhs.mat; type=rhs.type; number=rhs.number;	}

};


/********************************************************/
/*			THe Variable Maps and Typedefs 			    */
/********************************************************/

typedef std::map<std::string, complex > CVarMap;
typedef std::map<std::string, complex >::iterator CVarMapIter;
typedef std::pair<std::string, complex> aCVar;

extern CVarMap HamiltonianGenGlobalVars; //the GLOABL variables object

/******************************/
/****** Global Vars *******/
/******************************/
#define N_pi  3.141592653589793238
#define N_e   2.718281828459045235


/**** THe MAster HamiltonianGen Class ****/
class HamiltonianGen
{
	private:
		//the error found global var
		static char SYMBOL[15];
		static char MATHSYMBOL[13];

		Vector<HamiltonianGenStack> NStack; //the working stack evalutations
		int tos; //the stack posistion

		//some numerical limits
		static double lnlim;      // Smallest number allowed
		static double loglim;    // Smallest number allowed in call to log10() *

		bool ErrorFound;		//error flag

		int n;	//holds the current parser string posiston
		std::string chunk; //the piece in the parser we are on.
		std::string SymbGroup; //the current 'type'

		Vector<HamiltonianGenStackElement> usingStack; //the master stack for this clas

		CVarMap internalVars; //the internal variables

		std::string InExpr; //the input expression

	//Variable testing
	//is int input a function
		int isMatFunc(std::string matf);
	//is the input a variable
		int isVariable(std::string variable);
	//is the input a spin operator
		int isSpinOp(std::string variable);
	//is the input a Spactial operator
		int isSpaceOp(std::string variable);
	//is the input a 'math' operation
		int isSymbol(char ch);
	//is the input a 'math' operation
		int isMathSymbol(char ch);

	//this splits a spinOp into its correct pieces
		HamiltonianGenSpinOp chopSpinOp(std::string in);

	//this returs the spaceop value
		matrix doSpinOp(SpinSys &A, HamiltonianGenSpinOp &sp);

	//this splits a Space into its correct pieces
		HamiltonianGenSpaceOp chopSpaceOp(std::string in);

	//this returs the spaceop value
		complex doSpaceOp(HamiltonianGenSpaceOp &sp, double theta, double phi, double gamma);

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
		void Push(HamiltonianGenStack j);
		//get the current input, and drease the stack size
		void Pop(HamiltonianGenStack &j);
		//updates the 'Nstack' with the computed value of
		// MF, removing the 'old' values from the stack
		bool Compute(std::string MF);

	//some math operations cannot be completed on
	// matricies (things like >, <, etc) this throws
	// an error if the user are doing soemthing like
		void matErr(std::string thing);

	public:

		HamiltonianGen();
		HamiltonianGen(std::string inp);
		HamiltonianGen(const HamiltonianGen &in);

		void operator=(const HamiltonianGen &in);

	//the master parsing function...
		bool parse(std::string expr);

	//the master parsing function...
		bool parse(const char *expr);

	//add an interanl variable
		void addVar(std::string name, complex val);
		void updateVar(std::string name, complex val);
		complex getVar(std::string name);
		void addGlobalVar(std::string name, complex val);
	//returns the value for the proper input variables
	//the vars 'x', 'y', 'z', and 't'
	// are internal vars iniitally set to 0
		matrix evaluate(SpinSys &A);
		matrix evaluate(SpinSys &A, double theta, double phi, double gamma=0.0);

		inline matrix operator()(SpinSys &A)
		{	return evaluate(A);	}

		inline matrix operator()(SpinSys &A, double theta, double phi, double gamma=0.0)
		{	return evaluate(A, theta, phi, gamma);	}

		inline matrix Hamiltonian(SpinSys &A)
		{	return evaluate(A);	}

		inline matrix Hamiltonian(SpinSys &A, double theta, double phi, double gamma=0.0)
		{	return evaluate(A, theta, phi, gamma);	}

		inline matrix Hamiltonian(SpinSys &A, std::string p)
		{
			if(p!=InExpr) parse(p);
			return evaluate(A);
		}

		inline matrix Hamiltonian(SpinSys &A, std::string p, double theta, double phi, double gamma=0.0)
		{
			if(p!=InExpr) parse(p);
			return evaluate(A, theta, phi, gamma);
		}

		//will be false if a failure happen
		inline operator bool()	{	return ErrorFound;	}
		inline bool fail()	{	return ErrorFound;	}
};

END_BL_NAMESPACE



#endif

