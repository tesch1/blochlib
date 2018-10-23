/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-4-01
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

/********************************

	Parameter grabber class....can suck out pararmeters from a file...
	that exsists in specific chunks

	this class has various options for turning off and on certin input
	file format bits...explained in code

	an input file can be like so... where 'moo' is a section...

	moo{
		pars1 34
		par2 54
	}

	moo2{
		par1 567
		par 34
	}


	or like so (with no section and '=' parameter separators)

	par=90
	par=78

	or combos...

********************************/


#ifndef _params_h_
#define _params_h_ 1


#include <list>
#include <string>
#include <map>
#include <fstream>
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"


BEGIN_BL_NAMESPACE

class Parameters{

	private:
		std::map<std::string, Vector<std::string> > filechop; //this hold each file 'section'


	public:
		int MaxLineLength;
		std::string secopen; //default is '{'
		std::string secclose; //default is '}'

		char paramsep; //default is ' ' (a white space)

		bool interactive; //true-->if parmeter not found will ask for it in console

		Parameters():
			MaxLineLength(10000),
			secopen("{"), secclose("}"), paramsep('\0'),
			interactive(true)
		{};

		Parameters(const Parameters &cp):
			filechop(cp.filechop), MaxLineLength(cp.MaxLineLength),
			secopen(cp.secopen), secclose(cp.secclose), paramsep(cp.paramsep),
			interactive(cp.interactive)
		{}

		//assumes that 'in' is an entire 'file'
		Parameters(const Vector<std::string> &in,
			std::string open="{", std::string cl="}", char psep='\0', bool interac=true);

		//assumes the input vector is a section 'flag'
		Parameters(std::string flag, Vector<std::string> &in,
			std::string open="{", std::string cl="}", char psep='\0',bool interac=true);

		//assumes the input std::string is a file name to be read..
		Parameters(std::string file,
			std::string open="{", std::string cl="}", char psep='\0',bool interac=true);

		//will read the file according to the 'rules'
		Parameters(std::ifstream &readme,
			std::string open="{", std::string cl="}", char psep='\0',bool interac=true);

		//slurps an entire file into a VEctor<std::string>
		Vector<std::string> SuckFile(std::string &infile);
		Vector<std::string> SuckFile(std::ifstream &infile);
		void read(std::string &infile);
		void read(std::ifstream &infile);

		//std::string FindFlag(Vector<std::string> &in);
		Vector<std::string> GetInside(std::string flag,
								Vector<std::string> &from,
								std::string seco="{",
								std::string secc="}");


		//get a double parameter 'par' from the section 'sec'
		double getParamD(std::string par, std::string sec="", bool required=true, double def=0);

		//get a int parameter 'par' from the section 'sec'
		int getParamI(std::string par, std::string sec="", bool required=true, int def=0);

		//get a std::string parameter 'par' from the section 'sec'
		std::string getParamS(std::string par, std::string sec="", bool required=true, std::string def="");

		//get a char parameter 'par' from the section 'sec'
		char getParamC(std::string par, std::string sec="", bool required=true, char def='\0');

		//get a coord<> parameter 'par' from the section 'sec'
		coord<> getParamCoordD(std::string par, std::string sec="", char sep=',', bool required=true, coord<> def=coord<>(0.0,0.0,0.0));

		//get a coord<int> parameter 'par' from the section 'sec'
		coord<int,3> getParamCoordI(std::string par, std::string sec="", char sep=',', bool required=true,  coord<int> def=coord<int>(0,0,0));

		//get a Vector<double> parameter 'par' from the section 'sec'
		Vector<double> getParamVectorD(std::string par, std::string sec="", char sep=',', bool required=true,  Vector<double> def=Vector<double>());

		//get a coord<int> parameter 'par' from the section 'sec'
		Vector<int> getParamVectorI(std::string par, std::string sec="", char sep=',', bool required=true,  Vector<int> def=Vector<int>());

		//get a coord<int> parameter 'par' from the section 'sec'
		Vector<coord<> > getParamVectorCoordD(std::string par, std::string sec="", char sep=',', bool required=true,  Vector<coord<> > def=Vector<coord<> >());

		//returns the section 'sec' from the map if it exsits.
		Vector<std::string> &section(std::string sec);

		//returns the section 'sec' from the map if it exsits.
		bool addSection(std::string sec);

		//will read the file
		void operator=( std::ifstream &in);
		//will read the file
		void operator=(const Parameters &in);

		inline bool empty(){	return filechop.empty();		}

		//this gives one the ability to set (i.e. change or create) a parameter in a 'seciton'
		void setParam(std::string par,  int toset, std::string sec="");
		void setParam(std::string par,  double toset, std::string sec="");
		void setParam(std::string par,  std::string toset, std::string sec="");
		void setParam(std::string par,  char toset, std::string sec="");

		void setParam(std::string par,  Vector<double> toset, std::string sec="");
		void setParam(std::string par,  Vector<int> toset, std::string sec="");

		template<int N>
		void setParam(std::string par,  coord<double, N> toset, std::string sec="")
		{
			Vector<std::string> tmp;
			Vector<std::string> *vsec=&section(sec);
			for(int i=0;i<vsec->size();i++)
			{
				if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
				else{ tmp=parse_param(vsec->get(i), paramsep); }
				if(tmp.size()>0){
					if(tmp[0]==par)
					{
						if(paramsep=='\0')	tmp[1]=" ";
						else tmp[1]=paramsep;
						for(int j=0;j<N;++j){
							tmp[1]+=dbtost(toset[j]);
							if(j!=N-1) tmp[1]+=",";
						}
						(*vsec)[i]=collapsVS(tmp);
						return;
					}

				}
			}
			//the param as not in the list..so add it
			std::string in=par;
			if(paramsep=='\0')	in+=" ";
			else in+=paramsep;

			for(int j=0;j<N;++j){
				in+=dbtost(toset[j]);
				if(j!=N-1) in+=",";
			}
			vsec->push_back(in);
		}

		template<int N>
		void setParam(std::string par,  coord<int,N> toset, std::string sec="")
		{
			Vector<std::string> tmp;
			Vector<std::string> *vsec=&section(sec);
			for(int i=0;i<vsec->size();i++)
			{
				if(paramsep=='\0'){ tmp=parse_param(vsec->get(i));	}
				else{ tmp=parse_param(vsec->get(i), paramsep); }
				if(tmp.size()>0){
					if(tmp[0]==par)
					{
						if(paramsep=='\0')	tmp[1]=" ";
						else tmp[1]=paramsep;
						for(int j=0;j<N;++j){
							tmp[1]+=itost(toset[j]);
							if(j!=N-1) tmp[1]+=",";
						}
						(*vsec)[i]=collapsVS(tmp);
						return;
					}

				}
			}
			//the param as not in the list..so add it
			std::string in=par;
			if(paramsep=='\0')	in+=" ";
			else in+=paramsep;

			for(int j=0;j<N;++j){
				in+=itost(toset[j]);
				if(j!=N-1) in+=",";
			}
			vsec->push_back(in);

		}


		//writes the parameters to a file
		bool print(std::ostream &out);
		bool print(std::string out);
		bool print(const char *out);

};

std::ostream &operator<<(std::ostream &oo, Parameters &out);

END_BL_NAMESPACE


#endif

