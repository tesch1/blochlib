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
  ham_gen.cc --> this class takes in a std::string, parses it and
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

#ifndef __HamiltonianGen_Data_cc__
#define __HamiltonianGen_Data_cc__ 1

#include "ham_gen.h"

using namespace std;

BEGIN_BL_NAMESPACE



//the static constants
//the error found global var
char HamiltonianGen::SYMBOL[15] = "+-*/%^=()<>&|!";
char HamiltonianGen::MATHSYMBOL[13] = "+-*/%^=<>&|!";

//some numerical limits
double HamiltonianGen::lnlim=1.0e-36;      // Smallest number allowed
double HamiltonianGen::loglim=1.0e-10 ;    // Smallest number allowed in call to log10() *

CVarMap HamiltonianGenGlobalVars;

/*** Constrcutors ***/
HamiltonianGen::HamiltonianGen():
	ErrorFound(false), n(0), chunk(""),SymbGroup(""),InExpr("")
{
	CVarMapIter jj=HamiltonianGenGlobalVars.find("e");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("e", N_e));
	else HamiltonianGenGlobalVars["e"]=N_e;

	jj=HamiltonianGenGlobalVars.find("pi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("pi", N_pi));
	else HamiltonianGenGlobalVars["pi"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("Pi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("Pi", N_pi));
	else HamiltonianGenGlobalVars["Pi"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("PI");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("PI", N_pi));
	else HamiltonianGenGlobalVars["PI"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("deg2rad");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("deg2rad", DEG2RAD));
	else HamiltonianGenGlobalVars["deg2rad"]=DEG2RAD;

	jj=HamiltonianGenGlobalVars.find("rad2deg");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("rad2deg", RAD2DEG));
	else HamiltonianGenGlobalVars["rad2deg"]=RAD2DEG;

	jj=HamiltonianGenGlobalVars.find("complexi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("complexi", complex(0.0,1.0)));
	else HamiltonianGenGlobalVars["complexi"]= complex(0.0,1.0);

	NStack.resize(50);
}

HamiltonianGen::HamiltonianGen(std::string inp):
	ErrorFound(false), n(0), chunk(""),SymbGroup(""),InExpr("")
{
	CVarMapIter jj=HamiltonianGenGlobalVars.find("e");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("e", N_e));
	else HamiltonianGenGlobalVars["e"]=N_e;

	jj=HamiltonianGenGlobalVars.find("pi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("pi", N_pi));
	else HamiltonianGenGlobalVars["pi"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("Pi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("Pi", N_pi));
	else HamiltonianGenGlobalVars["Pi"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("PI");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("PI", N_pi));
	else HamiltonianGenGlobalVars["PI"]=N_pi;

	jj=HamiltonianGenGlobalVars.find("deg2rad");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("deg2rad", DEG2RAD));
	else HamiltonianGenGlobalVars["deg2rad"]=DEG2RAD;

	jj=HamiltonianGenGlobalVars.find("rad2deg");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("rad2deg", RAD2DEG));
	else HamiltonianGenGlobalVars["rad2deg"]=RAD2DEG;

	jj=HamiltonianGenGlobalVars.find("complexi");
	if(jj==HamiltonianGenGlobalVars.end())	HamiltonianGenGlobalVars.insert(aCVar("complexi", complex(0.0,1.0)));
	else HamiltonianGenGlobalVars["complexi"]= complex(0.0,1.0);

	NStack.resize(50);
	parse(inp);
}

HamiltonianGen::HamiltonianGen(const HamiltonianGen &in)
{
	*this=in;
}

void HamiltonianGen::operator=(const HamiltonianGen &in)
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

void HamiltonianGen::addVar(std::string name, complex val)
{
	CVarMapIter i=internalVars.find(name);
	if(i==internalVars.end())
	{
		CVarMapIter jj=HamiltonianGenGlobalVars.find(name);
		if(jj==HamiltonianGenGlobalVars.end()){
			internalVars.insert(aCVar(name, val));
		}else{
			HamiltonianGenGlobalVars[name]=val;
		}
	}else{
		internalVars[name]=val;
	}
}

void HamiltonianGen::addGlobalVar(std::string name, complex val)
{
	CVarMapIter jj=HamiltonianGenGlobalVars.find(name);
	if(jj==HamiltonianGenGlobalVars.end()){
		internalVars.insert(aCVar(name, val));
	}else{
		HamiltonianGenGlobalVars[name]=val;
	}
}

void HamiltonianGen::updateVar(std::string name, complex val)
{
	addVar(name, val);
}

complex HamiltonianGen::getVar(std::string name)
{
	CVarMapIter i=internalVars.find(name);
	if(i==internalVars.end())
	{
		CVarMapIter jj=HamiltonianGenGlobalVars.find(name);
		if(jj==HamiltonianGenGlobalVars.end()){
			std::cerr<<"HamiltonianGen::Warning...variable "<<name<<"not found, setting to 0"<<std::endl;
			return 0;
		}
		return jj->second;
	}
	return i->second;
}


//find a function
int HamiltonianGen::isMatFunc(std::string matf)
{

  if      ( matf=="abs"   ) return(1);
  else if ( matf=="acos"   ) return(1);
  else if ( matf=="asin"  ) return(1);
  else if ( matf=="atan"   ) return(1);
  else if ( matf=="cos"    ) return(1);
  else if ( matf=="real"    ) return(1);
  else if ( matf=="imag"    ) return(1);
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
int HamiltonianGen::isVariable(std::string variable)
{
	//find the interanl vars first
	CVarMapIter i=internalVars.find(variable);
	if(i==internalVars.end()) {
		//see if it is in the global one
		CVarMapIter jj=HamiltonianGenGlobalVars.find(variable);
		if(jj==HamiltonianGenGlobalVars.end())	return 0;
		else return 1;
	}else return 1;

}

//is the input a 'math' operation
int HamiltonianGen::isSymbol(char ch)
{
	int j = 0;
	while (j<=14 && ch!=SYMBOL[j]) j++;
	return (j<15)?1:0;
}

//is the input a 'math' operation
int HamiltonianGen::isMathSymbol(char ch)
{
	int j = 0;
	while (j<=12 && ch!=MATHSYMBOL[j]) j++;
	return (j<13)?1:0;
}

//is the input a spin operator
// the generic form can be either
// T##_#,#
// T#m#_#,#
// I(x,y,z,m,p,e)_#
// I(x,y,z,m,p,e)
//this splits a spinOp into its correct pieces
HamiltonianGenSpinOp HamiltonianGen::chopSpinOp(std::string in)
{
	std::string tmS=removeWhite(in);
	HamiltonianGenSpinOp out;
	int gotm=0;
	char tmC[2]; tmC[0]='\0'; tmC[1]='\0';
	if(tmS.size() <2) return out; //not a spinop
	if(tmS[0]=='T'){
		out.type='T';
		if(tmS.size()<5) return HamiltonianGenSpinOp(); //not a T tensor

		if(isdigit(tmS[1])){ tmC[0]=tmS[1]; out.rankl=std::atoi(tmC);	}//the fisr number

		if(tmS[2]=='m' && isdigit(tmS[3])){ gotm=1; tmC[0]=tmS[3]; out.rankm=-1*std::atoi(tmC);	}
		else if(isdigit(tmS[2])){	gotm=0;tmC[0]=tmS[2];  out.rankm=std::atoi(tmC);	}
		else{ return HamiltonianGenSpinOp();	} //not a T type tensor

		if(tmS[3+gotm]!='_') return HamiltonianGenSpinOp(); //not a spin op

		if(isdigit(tmS[4+gotm])){	tmC[0]=tmS[4+gotm]; out.spin1=std::atoi(tmC);	}

		if(int(tmS.size())<(7+gotm)){ out.spin2=out.spin1; return out;}

		if(tmS[5+gotm]!=',') return out;

		if(isdigit(tmS[6+gotm])){ tmC[0]=tmS[gotm+6]; out.spin2=std::atoi(tmC);	}
		return out;
	}else if(tmS[0]=='I'){
		out.type='I';
		if(tmS[1]=='x' || tmS[1]=='y' ||tmS[1]=='z'|| tmS[1]=='m' || tmS[1]=='p' ||tmS[1]=='e'){
			out.cart=tmS[1];
		}else{
			return HamiltonianGenSpinOp();
		}

		if(tmS.size()>=4){
			if(tmS[2]!='_') return HamiltonianGenSpinOp();
			if(isdigit(tmS[3])){ tmC[0]=tmS[3]; out.spin1=std::atoi(tmC);	}
		}
		return out;
	}else{
		return HamiltonianGenSpinOp();
	}
}


int HamiltonianGen::isSpinOp(std::string variable)
{
	HamiltonianGenSpinOp tm=chopSpinOp(variable);
	return tm.isOp();
}

//this returs the spinop value
matrix HamiltonianGen::doSpinOp(SpinSys &A, HamiltonianGenSpinOp &sp)
{
	switch(sp.type)
	{
		case 'T':
			switch(sp.rankl)
			{
				case 0:	return T0(A, sp.spin1);
				case 1: return T1(A, sp.spin1, sp.rankm);
				case 2: return T2(A, sp.spin1, sp.spin2, sp.rankm);
				default:
					BLEXCEPTION(std::string(" Spin Tensor Rank ")+itost(sp.rankl)+ std::string(" Not Available...."))
			}
		case 'I':
			switch(sp.spin1)
			{
				case -1:
					switch(sp.cart)
					{
						case 'x': return A.Fx();
						case 'y': return A.Fy();
						case 'z': return A.Fz();
						case 'e': return A.Fe();
						case 'p': return A.Fp();
						case 'm': return A.Fmi();
						default:
							BLEXCEPTION(std::string(" Spin Tensor Rank ")+itost(sp.cart)+ std::string(" Not Available...."))


					}
				default:
					switch(sp.cart)
					{
						case 'x': return A.Ix(sp.spin1);
						case 'y': return A.Iy(sp.spin1);
						case 'z': return A.Iz(sp.spin1);
						case 'e': return A.Ie(sp.spin1);
						case 'p': return A.Ip(sp.spin1);
						case 'm': return A.Imi(sp.spin1);
						default:
							BLEXCEPTION(std::string(" Spin Tensor Rank ")+itost(sp.cart)+ std::string(" Not Available...."))


					}
			}
		default:
			BLEXCEPTION(std::string(" Spin Tensor Rank ")+itost(sp.type)+ std::string(" Not Available...."))
			return A.Fe();
	}
}

//is the input a spin operator
// the generic form can be either
// A##
// A#m#
//this splits a space Op into its correct pieces
HamiltonianGenSpaceOp HamiltonianGen::chopSpaceOp(std::string in)
{
	std::string tmS=removeWhite(in);
	HamiltonianGenSpaceOp out;
	int gotm=0;
	char tmC[2]; tmC[0]='\0'; tmC[1]='\0';
	if(tmS.size() <3) return out; //not a spinop
	if(tmS[0]=='A'){

		if(isdigit(tmS[1])){ tmC[0]=tmS[1]; out.rankl=std::atoi(tmC); }//the fisr number
		else{	return HamiltonianGenSpaceOp();	}

		if(tmS.size()>3){
			if(tmS[2]=='m' && isdigit(tmS[3])){ tmC[0]=tmS[3]; gotm=1; out.rankm=-1*std::atoi(tmC);	}
		}else if(isdigit(tmS[2])){
			 tmC[0]=tmS[2];
			 out.rankm=std::atoi(tmC);
		}else{
			return HamiltonianGenSpaceOp();
		}

		return out;
	}else{
		return HamiltonianGenSpaceOp();
	}
}
//is the input a Spactial operator
int HamiltonianGen::isSpaceOp(std::string variable)
{
	HamiltonianGenSpaceOp tm=chopSpaceOp(variable);
	return tm.isOp();
}

//this returs the spaceop value
complex HamiltonianGen::doSpaceOp(HamiltonianGenSpaceOp &sp, double theta, double phi, double gamma)
{
	switch(sp.rankl)
	{
		case 0:	return complex(1.0,0.0);
		case 1:	return A1(theta, phi, sp.rankm);
		case 2:	return A2(theta, phi, sp.rankm);
		default:
			BLEXCEPTION(std::string(" Spin Tensor Rank ")+itost(sp.rankl)+ std::string(" Not Available...."))
			return complex(0.0,0.0);
	}
}


//The main HamiltonianGen
bool HamiltonianGen::parse(const char *expr)
{	return parse(std::string(expr));	}

//The main HamiltonianGen
bool HamiltonianGen::parse(std::string expr)
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
void HamiltonianGen::GetNext()
{
	chunk="";

	if(n>=int(InExpr.size())) return;

	while (isspace(InExpr[n])) n++;

	if (isSymbol(InExpr[n])){
		SymbGroup="SYMBOL";
		chunk+=InExpr[n];
		n++;
		if(n<int(InExpr.size()) && isMathSymbol(InExpr[n-1]) && isMathSymbol(InExpr[n]) ){
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
		if (isMatFunc(chunk)){ SymbGroup="MATFUNC";	}
		else if (isVariable(chunk)){ SymbGroup="VARIABLE";	}
		else if(isSpinOp(chunk)){ SymbGroup="SPINOP";	}
		else if(isSpaceOp(chunk)){ SymbGroup="SPACEOP";	}
		else{
			ErrorFound = 1;
			BLEXCEPTION(std::string("HamiltonianGen::Prase Error...thing \"")+chunk+std::string("\" for \"")+InExpr+std::string("\" not found"))
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
//cout<<std::endl<<"first: "<<chunk<<" | "<<InExpr<<" | "<<n<<" | "<<SymbGroup<<std::endl;

}

//the first round of parsing
void HamiltonianGen::First()
{

  GetNext();
  if (!chunk.empty() && !ErrorFound) Third();
  else ErrorFound = 1;

}

void HamiltonianGen::Third()
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
		HamiltonianGenStackElement Stack;
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


void HamiltonianGen::Fourth()
{
	char sign;

	Fifth();
	sign = chunk[0];
	while( sign == '*' || sign == '/' ) {
		GetNext();
		Fifth();
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		if( sign == '*' ) Stack.OpSpecifier.MFun= "MULT" ;
		else Stack.OpSpecifier.MFun= "DIV" ;
		usingStack.push_back(Stack); //put the thing in the stack
		sign = chunk[0];
	}

} // End Fourth


// This routine keeps parsing when a '^' symbol is recognized.
void HamiltonianGen::Fifth()
{

	Sixth();
	if( chunk[0] == '^' ) {
		GetNext();
		Fifth();
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "pow" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

} /* End Fifth */


// This routine keeps parsing while a '+' symbol or '-' symbol is
// OPerator (not change of sign like 'Third()')
void HamiltonianGen::Sixth()
{
	char sign;

	sign = ' ';
	if(SymbGroup== "SYMBOL" && chunk[0] == '+' || chunk[0] == '-'  | chunk[0] == '!') {
		sign = chunk[0];
		GetNext();
	}

	Seventh();

	if( sign == '-' ) {
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "CHSSGN" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

	if( sign == '!' ) {
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "MATFUNC" ;
		Stack.OpSpecifier.MFun= "NOT" ;
		usingStack.push_back(Stack); //put the thing in the stack
	}

} // End Sixth


//keep parseing through a '(' ')' set
void HamiltonianGen::Seventh()
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
void HamiltonianGen::Eighth()
{

	if ( SymbGroup== "NUMBER" )
	{
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "NUMBER" ;
		Stack.OpSpecifier._Value =std::atof(chunk.c_str());
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else if ( SymbGroup== "VARIABLE"  ) {
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "VARIABLE" ;
		Stack.OpSpecifier.VarName =  chunk ;
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else if (SymbGroup=="SPINOP"){
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "SPINOP" ;
		Stack.OpSpecifier.spinOp = chopSpinOp(chunk) ;
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else if (SymbGroup=="SPACEOP"){
		HamiltonianGenStackElement Stack;
		Stack.ElNo = usingStack.size()-1;
		Stack.OpType= "SPACEOP" ;
		Stack.OpSpecifier.spaceOp = chopSpaceOp(chunk) ;
		usingStack.push_back(Stack); //put the thing in the stack
		GetNext();
	}
	else Ninth();

} // End Eighth


// This routine keeps parsing while a mathematical function is recognized
void HamiltonianGen::Ninth()
{
	if (  SymbGroup== "MATFUNC"  ) {
		std::string TempFunk =  chunk ;
		GetNext();
		if( chunk[0] == '(' ) {
			GetNext();
			Third();
			if (chunk[0] != ')') ErrorFound = 1;
			HamiltonianGenStackElement Stack;
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
void HamiltonianGen::Push(HamiltonianGenStack j)
{
  NStack[tos] = j;
  tos++;
}

//get the current input, and drease the stack size
void HamiltonianGen::Pop(HamiltonianGenStack &j)
{
  tos--;
  j = NStack[tos];
}

void HamiltonianGen::matErr(std::string thing)
{
	BLEXCEPTION(std::string(" the operation \"")+thing+std::string("\" is not alowed on matrices..."))
}


//updates the 'Nstack' with the computed value of
// MF, removing the 'old' values from the stack
bool HamiltonianGen::Compute (std::string MF)
{
	//static int SINGULAR=0;
	HamiltonianGenStack t1, t2, tmStack;

	if(MF=="PLUS") {
		Pop(t1); //get the rhs (reduce stack)
		Pop(t2); //get the lhs (reduce stack)
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=t2.number+t1.number;
		}else if(t1.type=='M' && t2.type=='C'){
			tmStack.type='M';
			tmStack.mat=t1.mat+t2.number;
		}else if(t1.type=='C' && t2.type=='M'){
			tmStack.type='M';
			tmStack.mat=t1.number+t2.mat;
		}else{
			tmStack.type='M';
			tmStack.mat=t1.mat+t2.mat;
		}
		Push(tmStack); //add the new computed value to the stack
	}else if( MF=="MINUS") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=t2.number-t1.number;
		}else if(t1.type=='M' && t2.type=='C'){
			tmStack.type='M';
			tmStack.mat=t2.number-t1.mat;
		}else if(t1.type=='C' && t2.type=='M'){
			tmStack.type='M';
			tmStack.mat=t2.mat-t1.number;
		}else{
			tmStack.type='M';
			tmStack.mat=t2.mat-t1.mat;
		}
		Push(tmStack); //add the new computed value to the stack
	}else if( MF=="GREATER") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number>t1.number), 0.0);
		}else{
			matErr(">");
		}
		Push(tmStack);
	}else if( MF=="LESS") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number<t1.number), 0.0);
		}else{
			matErr("<");
		}
		Push(tmStack);
	}else if( MF=="EQEQ") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number==t1.number), 0.0);
		}else{
			matErr("==");
		}
		Push(tmStack);
	}else if( MF=="NOT") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(!int(t1.number.Re()),!int(t1.number.Im()));
		}else{
			matErr("!");
		}
		Push(tmStack);
	}else if( MF=="NEEQ") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number!=t1.number), 0.0);
		}else{
			matErr("!=");
		}
		Push(tmStack);
	}else if( MF=="GREATEREQ") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number>=t1.number), 0.0);
		}else{
			matErr(">=");
		}
		Push(tmStack);
	}else if( MF=="LESSEQ") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number<=t1.number), 0.0);
		}else{
			matErr("<=");
		}
		Push(tmStack);
	}else if( MF=="ANDAND") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number.Re()&&t1.number.Re()), 0.0);
		}else{
			matErr("&&");
		}
		Push(tmStack);
	}else if( MF=="OROR") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(t2.number.Re()||t1.number.Re()), 0.0);
		}else{
			matErr("&&");
		}
		Push(tmStack);
	}else if( MF=="MOD") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=complex(double(int(t2.number.Re())%int(t1.number.Re())), 0.0);
		}else{
			matErr("&&");
		}
		Push(tmStack);
	}else if( MF=="MULT") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=t2.number*t1.number;
		}else if(t1.type=='M' && t2.type=='C'){
			tmStack.type='M';
			tmStack.mat=t1.mat*t2.number;
		}else if(t1.type=='C' && t2.type=='M'){
			tmStack.type='M';
			tmStack.mat=t1.number*t2.mat;
		}else{
			tmStack.type='M';
			tmStack.mat=t2.mat*t1.mat;
		}
		Push(tmStack); //add the new computed value to the stack
	}else if( MF=="DIV") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=t2.number/t1.number;
		}else if(t1.type=='M' && t2.type=='C'){
			tmStack.type='M';
			tmStack.mat=t2.number/t1.mat;
		}else if(t1.type=='C' && t2.type=='M'){
			tmStack.type='M';
			tmStack.mat=t2.mat/t1.number;
		}else{
			tmStack.type='M';
			tmStack.mat=t2.mat/t1.mat;
		}
		Push(tmStack); //add the new computed value to the stack
	}else if( MF=="pow") {
		Pop(t1);
		Pop(t2);
		if(t1.type=='C' && t2.type=='C'){
			tmStack.type='C';
			tmStack.number=pow(t2.number,t1.number);
		}else if(t1.type=='M' && t2.type=='C'){
			matErr("pow");
		}else if(t1.type=='C' && t2.type=='M'){
			matrix tmMat=t2.mat;
			if(t1.number.Re()>0){
				for(int i=0;i<int(t1.number.Re());++i)
					tmMat=tmMat*t2.mat;
			}else{
				for(int i=0;i<int(t1.number.Re());++i)
					tmMat=tmMat/t2.mat;
			}
			tmStack.type='M';
			tmStack.mat=tmMat;
		}else{
			matErr("pow");
		}
		Push(tmStack); //add the new computed value to the stack
	}else if ( MF=="abs") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=abs(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=abs(t1.mat);
		}
		Push(tmStack);
	}else if ( MF== "acos" ) {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=acos(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=acos(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="atan") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=atan(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=atan(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="real") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=t1.number.Re();
		}else{
			tmStack.type='M';
			tmStack.mat=Re(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="imag") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=t1.number.Im();
		}else{
			tmStack.type='M';
			tmStack.mat=Im(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="asin") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=asin(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=asin(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="cos") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=cos(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=cos(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="exp") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=exp(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=exp(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="int") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=double(int(t1.number.Re()));
		}else{
			matErr("int");
		}
		Push(tmStack);
	}else if( MF=="ln") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=log(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=log(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="log") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=log10(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=log10(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="ceil") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=ceil(t1.number.Re());
		}else{
			tmStack.type='M';
			tmStack.mat=ceil(Re(t1.mat));
		}
		Push(tmStack);
	}else if (MF=="sin") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=sin(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=sin(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="sqrt") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=sqrt(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=sqrt(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="tan") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=tan(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=tan(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="floor") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=floor(t1.number.Re());
		}else{
			tmStack.type='M';
			tmStack.mat=floor(Re(t1.mat));
		}
		Push(tmStack);
	}else if( MF=="sinh") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=sinh(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=sinh(t1.mat);
		}
		Push(tmStack);
	}else if(MF=="cosh") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=cosh(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=cosh(t1.mat);
		}
		Push(tmStack);
	}else if (MF=="tanh") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=tanh(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=tanh(t1.mat);
		}
		Push(tmStack);
	}else if( MF=="CHSSGN") {
		Pop(t1);
		if(t1.type=='C'){
			tmStack.type='C';
			tmStack.number=-(t1.number);
		}else{
			tmStack.type='M';
			tmStack.mat=-(t1.mat);
		}
		Push(tmStack);
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
matrix HamiltonianGen::evaluate(SpinSys &A )
{
	return evaluate(A,0., 0., 0.);
}

matrix HamiltonianGen::evaluate(SpinSys &A, double theta, double phi, double gamma)
{
	int i;
	tos    = 0;
	ErrorFound=false;
	HamiltonianGenStack tmStack;


	//Evaluate the stack for input value of x
	for (i=0;(i<usingStack.size()) && !(ErrorFound);i++) {
	//cout<<"Index: "<<i<<" | "
	//<<usingStack->Stack[i].OpType<<" | "
	//<< usingStack->Stack[i].OpSpecifier.MFun<<std::endl;

	// If stack element is a number, then add it to the calc stack
		if( usingStack[i].OpType== "NUMBER")
		{
			tmStack.type='C';
			tmStack.number=usingStack[i].OpSpecifier._Value;
			Push( tmStack );
		}
		// If stack element is a variable, then push it
		/// find it then push in in to the list
		else if( usingStack[i].OpType== "VARIABLE"  )
		{
			tmStack.type='C';
			tmStack.number=getVar(usingStack[i].OpSpecifier.VarName);
			Push(tmStack);
		}

		else if( usingStack[i].OpType== "SPACEOP"  )
		{
			tmStack.type='C';
			tmStack.number=doSpaceOp(usingStack[i].OpSpecifier.spaceOp, theta, phi, gamma);
			Push(tmStack);
		}

		else if( usingStack[i].OpType== "SPINOP"  )
		{
			tmStack.type='M';
			tmStack.mat=doSpinOp(A,usingStack[i].OpSpecifier.spinOp);
			Push(tmStack);
		}

		// If stack element is an operator or a function, then compute
		else if( usingStack[i].OpType== "MATFUNC")
		{
			ErrorFound=Compute( usingStack[i].OpSpecifier.MFun );
		}
	}

	ErrorFound = !(ErrorFound);
	return (NStack[0].mat); //return the last value calced
}

END_BL_NAMESPACE


#endif

