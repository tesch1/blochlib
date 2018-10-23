

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
 	scriptparser.cc --> this takes in a 'vector string'
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

*/
#ifndef __scriptparse_cc__
#define __scriptparse_cc__ 1

#include "scriptparse.h"


BEGIN_BL_NAMESPACE



/**************** Sequence Parser bits *************/
ScriptParse::ScriptParse(Parameters &pset)
{
	data=pset.section("");
	parse(data);
}

ScriptParse::ScriptParse(const Vector<std::string> &pset)
{
	data=pset;
	parse(data);
}

ScriptParse::ScriptParse(const ScriptParse &rhs)
{
	if(this==&rhs) return;
	*this=rhs;
}

void ScriptParse::operator=(const ScriptParse &rhs)
{
	if(this==&rhs) return;
	data=rhs.data;
	myParse=rhs.myParse;
}

bool ScriptParse::parse(Parameters &pset)
{
	data=pset.section("");
	return parse(data);
}



//this function parse the 'loop' keyword
int ScriptParse::parseLoops(Vector<std::string> &pset, int start)
{
	Vector<std::string> loops, varn;
	std::string loopArg, lhs="", rhs="";
	double step=1,  end=0;
	int loopEnd=0;
	bool gotloop=false;
	int i=0;
	for(i=start;i<pset.size();++i)
	{
		loopArg=removeWhite(pset[i]);
		if(loopArg.substr(0,5)=="loop("){
			loopArg=getInside(loopArg);
			loops=parse_param(loopArg, ':');
			//cout<<loops<<std::endl;
			if(loops.size()!=2 && loops.size()!=3 && loops[0]!="" && loops[1]!="")
			{
				std::string _mess = " Bad 'loop' usage for \""+pset[i]+"\"";
				_mess+=" should be \"loop(var=begin:end)\" or \"loop(var=begin:step:end)\", ";;
				BLEXCEPTION(_mess)
			}

			//get the '=' sign
			varn=parse_param(loops[0], '=');
			if(varn.size()!=2)
			{
				std::string _mess = " Bad 'loop' variable usage for \""+pset[i]+"\"";
				_mess+=" should be \"loop(var=begin:end)\" or \"loop(var=begin:step:end)\", ";;
				BLEXCEPTION(_mess)
			}
			lhs=varn[0];
			rhs=varn[1];
			if(loops.size()==2)
			{
				step=1.0;
				myParse.parse(loops[1]);
				end=myParse();
			}

			if(loops.size()==3)
			{
				myParse.parse(loops[1]);
				step=myParse();
				myParse.parse(loops[2]);
				end=myParse();
			}

			myParse.parse(rhs);
			myParse.addVar(lhs, myParse());
			//cout<<"LHS: "<<lhs<<" RHS: "<<rhs<<" "<<myParse.getVar(lhs)<<std::endl;
			//cout<<"step: "<<step<<" end: "<<end<<std::endl;
			gotloop=true;
		}

		if(gotloop && lhs!="")
		{
			//find the 'end' of the list
			Vector<std::string> LoopInside;
			int numL=0;
			//pset.print(cout, "\n");
			for(int KK=i+1;KK<pset.size();++KK)
			{
				loopArg=removeWhite(pset[KK]);
				if(pset[KK].find("{")<pset[KK].size() || pset[KK].find("}")<pset[KK].size())
				{
					BLEXCEPTION(" bad syntax (section ended ('{' '}' problem) before loop 'end' was found)")
				}

				if(loopArg.substr(0,5)=="loop(" || loopArg.substr(0,3)=="if(")
				{	numL++;		}

				if(loopArg==("end"))
				{
					if(numL<=0){ loopEnd=KK; i=KK; break; }
					numL--;
				}
				//cout<<"loop: "<<loopArg<<" | "<<numL<<" | "<<KK<<" | "<<loopEnd<<std::endl;

				LoopInside.push_back(pset[KK]);
			}
			if(loopEnd==0)
			{
				BLEXCEPTION(" bad syntax (section ended before loop 'end' was found)")
			}
			while(myParse.getVar(lhs)<=end)
			{
				//cout<<myParse.getVar(lhs)<<std::endl;
				//LoopInside.print(cout,"\n");
				parse(LoopInside); //recurse through the loop (in case we have another)
				myParse.updateVar(lhs, myParse.getVar(lhs)+step);
				//cout<<myParse.getVar(lhs)<<std::endl;
			}
			gotloop=false;
			return i;
		}
	}
	return start;
}

//this function parse the 'if' keyword
// returns the parsed index...
int ScriptParse::parseIfs(Vector<std::string> &pset, int start)
{
	Vector<std::string> Ifinside, Elinside;
	Vector<Vector<std::string> > ElIfinside;
	std::string loopArg;
	int iftest=0;
	Vector<int> elifTest;
	bool gotif=false, gotelif=false, gotel=false, dids=false;
	int i=0;
	//cout<<"*****************"<<std::endl;
	//pset.print(cout, "\n");
	for(i=start;i<pset.size();++i)
	{
		loopArg=removeWhite(pset[i]);
		gotelif=false;
		int numL=0;
		if(loopArg.substr(0,2)=="if"
		   && loopArg.substr(0,3)!="if(")
		{
			BLEXCEPTION(std::string(" Bad 'if' Syntax requires an argument")+
				"\n at line \""+pset[i]+"\"")
		}
//cout<<"Pre Anything: "<<pset[i]<<" | "<<numL<<std::endl;

		if(loopArg.substr(0,3)=="if(")
		{
			loopArg=getInside(loopArg);
			if(loopArg=="")
			{
				BLEXCEPTION(" 'if' requires an argument")
			}
			myParse.parse(loopArg);
			//std::cerr<<loopArg<<std::endl;
			iftest=int(myParse());
			//std::cerr<<iftest<<std::endl;
			gotif=true;
			//grab the stuff inside the if
			numL=0;
			for(int KK=i+1;KK<pset.size();KK++)
			{
				loopArg=removeWhite(pset[KK]);
				//dive into a new if if present
				//std::cerr<<loopArg<<" K: "<<KK<<" i:" <<i<<std::endl;

				if(loopArg.substr(0,5)=="loop(" || loopArg.substr(0,3)=="if(")
				{	numL++;		}


				//Ifinside.print(cout, "\n");
				if(loopArg==("end")){		numL--;	}

				if(numL<=0 && (loopArg.substr(0,4)=="else"))
				{ i=KK-1;	break;}
				if(numL==-1 && (loopArg.substr(0,3)=="end"))
				{  i=KK;	break;}

				Ifinside.push_back(pset[KK]);
			}
			if(loopArg.substr(0,3)=="end" && numL<=0) break;
		}

		if(i<pset.size()) loopArg=removeWhite(pset[i]);
//cout<<"After IF*********************** "<<i<<" "<<loopArg<<" numL: "<<numL<<std::endl;
		if(loopArg.substr(0,6)=="elseif")
		{
			Vector<std::string> tmL;
			if(loopArg.substr(0,7)!="elseif(")
			{
				BLEXCEPTION(std::string(" Bad 'else if' Syntax requires an argument")+
				"\n at line \""+pset[i]+"\"")
			}
			loopArg=getInside(loopArg);
			if(loopArg=="")
			{
				BLEXCEPTION(" 'else if' requires an argument")
			}
			myParse.parse(loopArg);
			elifTest.push_back(int(myParse()));
			gotelif=true;
			//grab the stuff inside the elseif
			numL=0;
			for(int KK=i+1;KK<pset.size();KK++)
			{
				loopArg=removeWhite(pset[KK]);
//cout<<"Pre ElseIF: "<<pset[KK]<<" | "<<numL<<" test: "<<test<<std::endl;
				//dive into a new if if present

				if(loopArg.substr(0,5)=="loop(" || loopArg.substr(0,3)=="if(" )
				{	numL++;		}

				if(loopArg==("end")){		numL--;		}

				if(numL<=0 && (loopArg.substr(0,4)=="else"))
				{ ElIfinside.push_back(tmL); i=KK-1; break;}
				if(numL==-1 && (loopArg.substr(0,3)=="end"))
				{ ElIfinside.push_back(tmL); i=KK; break;}

				tmL.push_back(pset[KK]);
			}
			if(loopArg.substr(0,3)=="end" && numL<=0) break;
		}

		if(i<pset.size()) loopArg=removeWhite(pset[i]);
		if(loopArg.substr(0,4)=="else" && !gotelif)
		{
			gotel=true;
			//if the if test failed we do this one
			//grab the stuff inside the elseif
			numL=0;
			for(int KK=i+1;KK<pset.size();KK++)
			{
				loopArg=removeWhite(pset[KK]);

				if(loopArg.substr(0,5)=="loop(" ||  loopArg.substr(0,3)=="if(" )
				{	numL++;		}

				if(loopArg==("end")){		numL--;		}
				if(numL==-1 && loopArg.substr(0,3)=="end")
				{	i=KK;	break;	}
				Elinside.push_back(pset[KK]);
			}
			break;
		}
	}
//now do the correct thing...
	int retend=Ifinside.size()+start;
	for(int jj=0;jj<ElIfinside.size();++jj) retend+=ElIfinside[jj].size();
	retend+=Elinside.size();
	if(iftest!=0){
		parse(Ifinside);
		dids=true;
	}

	if(!dids){
		for(int jj=0;jj<ElIfinside.size();++jj){
			if(elifTest[jj]!=0){
				parse(ElIfinside[jj]);
				dids=true;
				break;
			}
		}
	}

	if(!dids && !Elinside.empty()){
		parse(Elinside);
		dids=true;
	}
	if(gotif) return retend;
	return start;
}


bool ScriptParse::parse(const Vector<std::string> &inset)
{
	Vector<std::string> pset=paramStrip(inset);
	std::string tmP;
	//cout<<"***"<<std::endl;
	//pset.print(cout, "\n");
	for(int i=0;i<pset.size();++i)
	{
		//cout<<pset[i]<<" | "<<int(pset[i].find("if("))<<std::endl;
		tmP=removeWhite(pset[i]);
		//cout<<"I: "<<i<<std::endl;
		//pset.print(cout, "\n");
		if(tmP.substr(0,5)=="loop("){
			i=parseLoops(pset, i);
		}else if(tmP.substr(0,3)=="if("){
			i=parseIfs(pset, i);
		}else{
			decide(pset[i]);
		}
	}
	return true;
}

bool ScriptParse::parse()
{	return parse(data);	}

//this is the master decider function
// and snaps it to the proper above function
bool ScriptParse::decide(std::string inSqe)
{
	if(inSqe.find("=")<inSqe.size()
	   && inSqe.find("==")!=inSqe.find("=")
	   && inSqe.find("!=")!=inSqe.find("=")-1 ){
		addVar(inSqe);
	}else if(inSqe.find("print(")<inSqe.size()){
		doPrint(inSqe);
	}

	return true;
}


//this should take in a string like "A=B" and add the
//proper line to the sequence vector
void ScriptParse::addVar(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("=")<tmS.size()
	   && inSqe.find("==")!=inSqe.find("=")
	   && inSqe.find("!=")!=inSqe.find("=")-1 )
	{
		ps=parse_param(tmS, '=');
		if(ps.size()<=1)
		{
			BLEXCEPTION(std::string(" Bad 'variable' usage for \"")+inSqe+"\""
				"\n should be \"A=(b*...)\" ")
		}
		std::string lhs=ps[0];
		std::string rhs=ps[1];

		if(ps.size()>=1){
			myParse.parse(rhs);
			myParse.addVar(lhs, myParse());
		}
	}
}

//this should take in a string like "A=B" and add the
//proper line to the sequence vector
void ScriptParse::addGlobalVar(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("=")<tmS.size()
	   && inSqe.find("==")!=inSqe.find("=")
	   && inSqe.find("!=")!=inSqe.find("=")-1 )
	{
		ps=parse_param(tmS, '=');
		if(ps.size()<=1)
		{
			BLEXCEPTION(std::string(" Bad 'variable' usage for \"")+inSqe+"\""
				"\n should be \"A=(b*...)\" ")
		}
		std::string lhs=ps[0];
		std::string rhs=ps[1];

		if(ps.size()>=1){
			myParse.parse(rhs);
			myParse.addGlobalVar(lhs, myParse());
		}
	}
}

//prints the argument
void ScriptParse::doPrint(std::string inSqe)
{
	std::string out;
	Vector<std::string> ps;
	std::string tmS=removeWhite(inSqe);
	if(tmS.find("print(")<tmS.size())
	{
		tmS=getInside(tmS);
		if(tmS.size()<1)
		{
			BLEXCEPTION(std::string(" Bad 'print' usage for \"")+inSqe+"\""+
				"\n  should be \"print(thing)\" ")
		}

		myParse.parse(tmS);
		std::cout<<tmS<<"="<<myParse()<<std::endl;
	}
}
END_BL_NAMESPACE



#endif
