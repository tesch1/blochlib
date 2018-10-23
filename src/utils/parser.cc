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
  parser.cc --> this class takes in a string, parses it and
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

#ifndef __Parser_Data_cc__
#define __Parser_Data_cc__ 1

#include "parser.h"

BEGIN_BL_NAMESPACE



//the static constants
//the error found global var
char Parser::SYMBOL[15] = "+-*/%^=()<>&|!";
char Parser::MATHSYMBOL[13] = "+-*/%^=<>&|!";

//some numerical limits
double Parser::lnlim=1.0e-36;      // Smallest number allowed
double Parser::loglim=1.0e-10 ;    // Smallest number allowed in call to log10() *

VarMap ParserGlobalVars;

/*** Constrcutors ***/
Parser::Parser():
	ErrorFound(false), n(0), chunk(""),SymbGroup(""),InExpr("")
{
	VarMapIter jj=ParserGlobalVars.find("e");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("e", N_e));
	else ParserGlobalVars["e"]=N_e;

	jj=ParserGlobalVars.find("pi");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("pi", N_pi));
	else ParserGlobalVars["pi"]=N_pi;

	jj=ParserGlobalVars.find("Pi");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("Pi", N_pi));
	else ParserGlobalVars["Pi"]=N_pi;

	jj=ParserGlobalVars.find("PI");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("PI", N_pi));
	else ParserGlobalVars["PI"]=N_pi;

	jj=ParserGlobalVars.find("deg2rad");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("deg2rad", DEG2RAD));
	else ParserGlobalVars["deg2rad"]=DEG2RAD;

	jj=ParserGlobalVars.find("rad2deg");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("rad2deg", RAD2DEG));
	else ParserGlobalVars["rad2deg"]=RAD2DEG;

	internalVars.insert(aVar("x", 0));
	internalVars.insert(aVar("y", 0));
	internalVars.insert(aVar("z", 0));
	internalVars.insert(aVar("t", 0));
	NStack.resize(50,0);
}

Parser::Parser(std::string inp):
	ErrorFound(false), n(0), chunk(""),SymbGroup(""),InExpr("")
{
	VarMapIter jj=ParserGlobalVars.find("e");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("e", N_e));
	else ParserGlobalVars["e"]=N_e;

	jj=ParserGlobalVars.find("pi");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("pi", N_pi));
	else ParserGlobalVars["pi"]=N_pi;

	jj=ParserGlobalVars.find("Pi");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("Pi", N_pi));
	else ParserGlobalVars["Pi"]=N_pi;

	jj=ParserGlobalVars.find("PI");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("PI", N_pi));
	else ParserGlobalVars["PI"]=N_pi;

	jj=ParserGlobalVars.find("deg2rad");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("deg2rad", DEG2RAD));
	else ParserGlobalVars["deg2rad"]=DEG2RAD;

	jj=ParserGlobalVars.find("rad2deg");
	if(jj==ParserGlobalVars.end())	ParserGlobalVars.insert(aVar("rad2deg", RAD2DEG));
	else ParserGlobalVars["rad2deg"]=RAD2DEG;

	internalVars.insert(aVar("x", 0));
	internalVars.insert(aVar("y", 0));
	internalVars.insert(aVar("z", 0));
	internalVars.insert(aVar("t", 0));
	NStack.resize(50,0);
	parse(inp);
}

Parser::Parser(const Parser &in)
{
	*this=in;
}

void Parser::operator=(const Parser &in)
{
	internalVars=in.internalVars;
	ErrorFound=in.ErrorFound;
	n=in.n;
	chunk=in.chunk;
	SymbGroup=in.SymbGroup;
	InExpr=in.InExpr;
	tos=in.tos;
	usingStack=in.usingStack;
}

void Parser::addVar(std::string name, double val)
{
	VarMapIter i=internalVars.find(name);
	if(i==internalVars.end())
	{
		VarMapIter jj=ParserGlobalVars.find(name);
		if(jj==ParserGlobalVars.end()){
			internalVars.insert(aVar(name, val));
		}else{
			ParserGlobalVars[name]=val;
		}
	}else{
		internalVars[name]=val;
	}
}

void Parser::addGlobalVar(std::string name, double val)
{
	VarMapIter jj=ParserGlobalVars.find(name);
	if(jj==ParserGlobalVars.end()){
		ParserGlobalVars.insert(aVar(name, val));
	}else{
		ParserGlobalVars[name]=val;
	}
}

void Parser::updateVar(std::string name, double val)
{
	addVar(name, val);
}

double Parser::getVar(std::string name)
{
	VarMapIter i=internalVars.find(name);
	if(i==internalVars.end())
	{
		VarMapIter jj=ParserGlobalVars.find(name);
		if(jj==ParserGlobalVars.end()){
			std::cerr<<"Parser::Warning...variable "<<name<<"not found, setting to 0"<<std::endl;
			return 0;
		}
		return jj->second;
	}
	return i->second;
}

double Parser::getGlobalVar(std::string name)
{
	VarMapIter jj=ParserGlobalVars.find(name);
	if(jj==ParserGlobalVars.end()){
		std::cerr<<"Parser::Warning...variable "<<name<<"not found, setting to 0"<<std::endl;
		return 0;
	}
	return jj->second;
}

//find a function
int Parser::isMatFunc(std::string matf)
{

  if      ( matf=="abs"   ) return(1);
  else if ( matf=="acos"   ) return(1);
  else if ( matf=="asin"  ) return(1);
  else if ( matf=="atan"   ) return(1);
  else if ( matf=="cos"    ) return(1);
  else if ( matf=="exp"    ) return(1);
  else if ( matf=="floor"  ) return(1);
  else if ( matf=="frac"   ) return(1);
  else if ( matf=="int"    ) return(1);
  else if ( matf=="ln"     ) return(1);
  else if ( matf=="log"    ) return(1);
  else if ( matf=="pow"    ) return(1);
  else if ( matf=="round"  ) return(1);
  else if ( matf=="ceil"  ) return(1);
  else if ( matf=="sin"    ) return(1);
  else if ( matf=="sqrt"   ) return(1);
  else if ( matf=="tan"    ) return(1);
  else if ( matf=="trunc"  ) return(1);
  else if ( matf=="sinh"   ) return(1);
  else if ( matf=="cosh"   ) return(1);
  else if ( matf=="tanh"   ) return(1);
  else if ( matf=="CHSSGN" ) return(1);
  else if ( matf=="FAC"    ) return(1);
  else return(0);

}

//is the input a variable
int Parser::isVariable(std::string variable)
{
	//find the interanl vars first
	VarMapIter i=internalVars.find(variable);
	if(i==internalVars.end()) {
		//see if it is in the global one
		VarMapIter jj=ParserGlobalVars.find(variable);
		if(jj==ParserGlobalVars.end())	return 0;
		else return 1;
	}else return 1;

}

//is the input a 'math' operation
int Parser::isSymbol(char ch)
{
	int j = 0;
	while (j<=14 && ch!=SYMBOL[j]) j++;
	return (j<15)?1:0;
}

//is the input a 'math' operation
int Parser::isMathSymbol(char ch)
{
	int j = 0;
	while (j<=12 && ch!=MATHSYMBOL[j]) j++;
	return (j<13)?1:0;
}

//The main Parser
bool Parser::parse(const char *expr)
{	return parse(std::string(expr));	}

//The main Parser
bool Parser::parse(std::string expr)
{

  usingStack.resize(0);

  ErrorFound = false;
  n = 0;
  InExpr=expr;

  First();
  return(!ErrorFound);

}

//gets the 'next' thing in the expression  and attempts to figure out
// what it could be ...a function (MATFUNC), a variable (VARIABLE),
// a number (NUMBER), or a symbol (SYMBOL)
void Parser::GetNext()
{
	chunk="";

	if(n>=int(InExpr.size())) return;

	while (isspace(InExpr[n])) n++;

	if (isSymbol(InExpr[n])){
		SymbGroup="SYMBOL";
		chunk+=InExpr[n];
		n++;
		if(n<int(InExpr.size())&& isMathSymbol(InExpr[n-1])  && isMathSymbol(InExpr[n]) ){
			SymbGroup="SYMBOL";
			chunk+=InExpr[n];
			n++;
		}

	}else if (isalpha(InExpr[n])){
		while (n<int(InExpr.size()) && !isSymbol(InExpr[n]) && (InExpr[n-1]!=']' ))
		{
			if (!isspace(InExpr[n])) chunk+=(InExpr[n]);
			n++;
		}
		if (isMatFunc(chunk)) SymbGroup="MATFUNC";
		else if (isVariable(chunk)) SymbGroup="VARIABLE";
		else{
			ErrorFound = 1;
			BLEXCEPTION(std::string("Parser::Prase Error...thing \"")+chunk+"\" for \""+InExpr+"\" not found")
		}

	}else if(n<int(InExpr.size()) &&  (isdigit( InExpr[n] ) || (InExpr[n]=='.'))) {

		/*if(InExpr[n]=='.'){
			j = n+1;
			while ( InExpr[j] != '\0') { j++; }
			while (j>=n)
			{
				InExpr [j+1] = InExpr [j];
				j--;
			}
			InExpr[n]= '0' ;
		}*/
		while( n<int(InExpr.size())  )
		{
		   	if(isdigit(InExpr[n]) || InExpr[n]=='.'){
				chunk+=InExpr[n];
				SymbGroup="NUMBER";
				n++;
			}else if(InExpr[n]=='e'){
				if(n+1<int(InExpr.size())){
					if(InExpr[n+1]=='-' || InExpr[n+1]=='+' || isdigit( InExpr[n+1])){
						chunk+=InExpr[n];
						chunk+=InExpr[n+1];
						SymbGroup="NUMBER";
						n+=2;
					}
				}else{
					break;
				}
			}else{
				break;
			}
		}
	}
//cout<<endl<<"first: "<<chunk<<" | "<<InExpr<<" | "<<n<<" | "<<SymbGroup<<endl;

}

//the first round of parsing
void Parser::First()
{

  GetNext();
  if (!chunk.empty() && !ErrorFound) Third();
  else ErrorFound = 1;

}

void Parser::Third()
{
	char sign, secS='\0';

	Fourth();
	sign = chunk[0];
	if(chunk.size()>1) secS=chunk[1];

	while( sign == '+' || sign == '-'||
	   sign == '<' || sign == '>'|| sign == '%'||
	   sign == '&' || sign == '|' || sign == '!' || sign == '='
	   )
	{
		GetNext();
		Fourth();
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType="MATFUNC" ;
		if( sign == '+' ) Stack.OpSpecifier.MFun= "PLUS" ;
		else if(sign=='-') Stack.OpSpecifier.MFun= "MINUS" ;
		else if(sign=='>' && secS=='\0') Stack.OpSpecifier.MFun= "GREATER" ;
		else if(sign=='<' && secS=='\0') Stack.OpSpecifier.MFun= "LESS" ;
		else if(sign=='=' && secS=='=') Stack.OpSpecifier.MFun= "EQEQ" ;
		else if(sign=='!' && secS=='=') Stack.OpSpecifier.MFun= "NEEQ" ;
		else if(sign=='>' && secS=='=') Stack.OpSpecifier.MFun= "GREATEREQ" ;
		else if(sign=='<' && secS=='=') Stack.OpSpecifier.MFun= "LESSEQ" ;
		else if(sign=='&' && secS=='&') Stack.OpSpecifier.MFun= "ANDAND" ;
		else if(sign=='|' && secS=='|') Stack.OpSpecifier.MFun= "OROR" ;
		else if(sign=='%') Stack.OpSpecifier.MFun= "MOD" ;
		usingStack.push_back(Stack); //put the thing in the stack
		sign = chunk[0];
	}

} //End Third


void Parser::Fourth()
{
	char sign;

	Fifth();
	sign = chunk[0];
	while( sign == '*' || sign == '/' ) {
		GetNext();
		Fifth();
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		if( sign == '*' ) Stack.OpSpecifier.MFun= "MULT" ;
		else Stack.OpSpecifier.MFun= "DIV" ;
		usingStack.push_back(Stack); //put the thing in the stack
		sign = chunk[0];
	}

} // End Fourth


// This routine keeps parsing when a '^' symbol is recognized.
void Parser::Fifth()
{

	Sixth();
	if( chunk[0] == '^' ) {
		GetNext();
		Fifth();
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "pow" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

} /* End Fifth */


// This routine keeps parsing while a '+' symbol or '-' symbol is
// OPerator (not change of sign like 'Third()')
void Parser::Sixth()
{
	char sign;

	sign = ' ';
	if(SymbGroup== "SYMBOL" && chunk[0] == '+' || chunk[0] == '-'  | chunk[0] == '!') {
		sign = chunk[0];
		GetNext();
	}

	Seventh();

	if( sign == '-' ) {
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "CHSSGN" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

	if( sign == '!' ) {
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "NOT" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

} // End Sixth


//keep parseing through a '(' ')' set
void Parser::Seventh()
{

	if( chunk[0] == '(' && SymbGroup=="SYMBOL") {
		GetNext();
		Third(); //parse the insides until we get a ')'
		if (chunk[0] != ')') ErrorFound = 1;
		GetNext();
	}
	else Eighth();

} // End Seventh



// This routine keeps parsing while a number or variable is recognized.
void Parser::Eighth()
{

	if ( SymbGroup== "NUMBER" )
	{
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "NUMBER" ;
		Stack.OpSpecifier._Value =std::atof(chunk.c_str());
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else if ( SymbGroup== "VARIABLE"  ) {
		ParserStackElement<double> Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "VARIABLE" ;
		Stack.OpSpecifier.VarName =  chunk ;
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else Ninth();

} // End Eighth


// This routine keeps parsing while a mathematical function is recognized
void Parser::Ninth()
{
	if (  SymbGroup== "MATFUNC"  ) {
		std::string TempFunk =  chunk ;
		GetNext();
		if( chunk[0] == '(' ) {
			GetNext();
			Third();
			if (chunk[0] != ')') ErrorFound = 1;
			ParserStackElement<double> Stack;
			Stack.ElNo = usingStack.size()-1;
			Stack.OpType="MATFUNC" ;
			Stack.OpSpecifier.MFun= TempFunk ;
			usingStack.push_back(Stack); //put the thing in the stack
			GetNext();
		}
	}else
		ErrorFound = 1;

} /* End Ninth */


/*******************************************************************/
/********** THE COMPUTATION FUNCTIONS *****************************/
//add the input to the stack
void Parser::Push(double j)
{
  NStack[tos] = j;
  tos++;
}

//get the current input, and drease the stack size
void Parser::Pop(double &j)
{
  tos--;
  j = NStack[tos];
}


//updates the 'Nstack' with the computed value of
// MF, removing the 'old' values from the stack
bool Parser::Compute (std::string MF)
{
	static int SINGULAR=0;
	double t1, t2;

	if(MF=="PLUS") {
		Pop(t1); //get the rhs (reduce stack)
		Pop(t2); //get the lhs (reduce stack)
		Push(t2+t1); //add the new computed value to the stack
	}else if( MF=="MINUS") {
		Pop(t1);
		Pop(t2);
		Push(t2-t1);
	}else if( MF=="GREATER") {
		Pop(t1);
		Pop(t2);
		Push(double(t2>t1));
	}else if( MF=="LESS") {
		Pop(t1);
		Pop(t2);
		Push(double(t2<t1));
	}else if( MF=="EQEQ") {
		Pop(t1);
		Pop(t2);
		Push(double(t2==t1));
	}else if( MF=="NOT") {
		Pop(t1);
		Push(double(!t1));
	}else if( MF=="NEEQ") {
		Pop(t1);
		Pop(t2);
		Push(double(t2!=t1));
	}else if( MF=="GREATEREQ") {
		Pop(t1);
		Pop(t2);
		Push(double(t2>=t1));
	}else if( MF=="LESSEQ") {
		Pop(t1);
		Pop(t2);
		Push(double(t2<=t1));
	}else if( MF=="ANDAND") {
		Pop(t1);
		Pop(t2);
		Push(double(int(t2)&&int(t1)));
	}else if( MF=="OROR") {
		Pop(t1);
		Pop(t2);
		Push(double(int(t2)||int(t1)));
	}else if( MF=="MOD") {
		Pop(t1);
		Pop(t2);
		Push(double(int(t2)%int(t1)));
	}else if( MF=="MULT") {
		Pop(t1);
		Pop(t2);
		if ((fabs(t1) > 10e+10) || (fabs(t2) > 10e+10)) {
			Push (pow(sqrt(t2)*sqrt(t1),2.0));
			if (SINGULAR) return true;    //a better mult for large values
			SINGULAR = 0;
		}else{
			Push(t2*t1);
		}
	}else if( MF=="DIV") {
		Pop(t1);
		Pop(t2);
		if ((fabs(t1) < loglim) && (fabs(t2) < loglim)) { //a divide by 0 error
			Push(0);
			SINGULAR = 1;
			return true;
		}else{
			if (fabs(t1) < lnlim) {  //underflow....
				 Push(0); return true;
	 		}else Push(t2/t1);
		}
	}else if( MF=="pow") {
		Pop(t1);
		Pop(t2);
		Push(pow(t2,t1));
	}else if ( MF=="abs") {
		Pop(t1);
		Push(fabs(t1));
	}else if ( MF== "acos" ) {
		Pop(t1);
		if (fabs(t1) > 1){  Push(0); return true;}
		Push(acos(t1));
	}else if (MF=="atan") {
		Pop(t1);
		Push(atan(t1));
	}else if( MF=="asin") {
		Pop(t1);
		if ( fabs(t1) > 1 ){	Push(0); return true;	}
		Push(asin(t1));
	}else if( MF=="cos") {
		Pop(t1);
		Push(cos(t1));
	}else if( MF=="exp") {
		Pop(t1);
		Push(exp(t1));
	}else if( MF=="frac") {
		Pop(t1);
		Push(t1-(int)t1);
	}else if( MF=="int") {
		Pop(t1);
		Push((double)((int)t1));
	}else if( MF=="ln") {
		Pop(t1);
		if (t1 > lnlim) Push(log(t1));
		else return true;
	}else if( MF=="log") {
		Pop(t1);
		if (t1 > loglim) Push(log10(t1));
		else return true;
	}else if (MF=="round") {
		Pop(t1);
		if (t1 >= 0.0) Push((int)(t1+0.5));
		else Push( (int)(t1- 0.5 ));
	}else if (MF=="ceil") {
		Pop(t1);
		Push( ceil(t1));
	}else if (MF=="sin") {
		Pop(t1);
		Push(sin(t1));
	}else if (MF=="sqrt") {
		Pop(t1);
		if (t1< 0.0 ){	Push(0);  return true;}
		Push(sqrt(t1));
	}else if (MF=="tan") {
		Pop(t1);
		if ( fabs(cos(t1)) < lnlim ){ Push(0);  return true;}
		else Push(tan(t1));
	}else if (MF=="trunc") {
		Pop(t1);
		Push( (int)t1);
	}else if( MF=="floor") {
		Pop(t1);
		Push(floor(t1));
	}else if( MF=="sinh") {
		Pop(t1);
		Push(sinh(t1));
	}else if(MF=="cosh") {
		Pop(t1);
		Push(cosh(t1));
	}else if (MF=="tanh") {
		Pop(t1);
		Push(tanh(t1));
	}else if( MF=="CHSSGN") {
		Pop(t1);
		Push(-t1);
	}else if( MF=="print") {
		Pop(t1);
		std::cout<<t1<<std::endl;
	}
	/*else  if (MF=="FAC") {
		int fp;
		Pop(t1);
		fp=floor(t1+0.5);
		if (t1<0.0)
			if (fabs(t1-fp)<1e-13) *ErrorEval = 1;
			else Push(Factorial(t1));
		else
			Push(Factorial(t1));
	}*/
	return false;
}


/**** the evaluators ****/

//returns the value for the proper input variables
//the vars 'x', 'y', 'z', and 't'
// are internal vars iniitally set to 0
double Parser::evaluate()
{
	return evaluate(0.,0.,0.,0.);
}

double Parser::evaluate(double x, double y, double z, double t)
{
	int i;
	tos    = 0;
	ErrorFound=false;
	//regesiter the 4 vars
	updateVar("x", x);
	updateVar("y", y);
	updateVar("z", z);
	updateVar("t", t);


	//Evaluate the stack for input value of x
	for (i=0;(i<usingStack.size()) && !(ErrorFound);i++) {
	//cout<<"Index: "<<i<<" | "
	//<<usingStack->Stack[i].OpType<<" | "
	//<< usingStack->Stack[i].OpSpecifier.MFun<<endl;

	// If stack element is a number, then add it to the calc stack
		if( usingStack[i].OpType== "NUMBER")
		{
			Push( usingStack[i].OpSpecifier._Value );
		}
		// If stack element is a variable, then push it
		/// find it then push in in to the list
		else if( usingStack[i].OpType== "VARIABLE"  )
		{
			double curV=getVar(usingStack[i].OpSpecifier.VarName);
			Push(curV);
		}

		// If stack element is an operator or a function, then compute
		else if( usingStack[i].OpType== "MATFUNC")
		{
			ErrorFound=Compute( usingStack[i].OpSpecifier.MFun );
		}
	}

	ErrorFound = !(ErrorFound);
	return (NStack[0]); //return the last value calced
}

END_BL_NAMESPACE


#endif

