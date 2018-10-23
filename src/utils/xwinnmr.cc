/*
 * Copyright (c)2000-2002 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-04-02
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
 	A Bruker XWINNMR file ('directory') reader....

 */


#ifndef XWINNMR_STREAM_CC_
#define XWINNMR_STREAM_CC_ 1



#include "utils/endians.h"
#include "utils/xwinnmr.h"
#include <iostream>

BEGIN_BL_NAMESPACE

using namespace BL_endians; //for the endian conversions


/***** Finds parameters in a generic Bruker 'acqus' file  ***/

//sets the current position in the
// file to that of the desired parameter
//returns true if it finds the parameter
//a parameter in the procpar looks like

/*

##$SW_h= 1124.10071942446
##$TD= 8192

*/
//this class gets a data chunk from
// the Varian 'proppar' file
//the procpar file contains the nessesary parameters
// for a the varian GUI to function properly
// it also contains the data (like sweep widths) and other
// data parameters NOT in the binary fid file

void XWINNMRacqu::FileErr() //a file read error
{
	std::cerr<<std::endl<<"Error: XWINNMRacqu..File read error"<<std::endl;
	std::cerr<<" file cannot be opened or read"<<std::endl;
}

void XWINNMRacqu::NotFoundErr() //parameter not found
{
	std::cerr<<std::endl<<"Error: XWINNMRacqu..Parameter Not Found"<<std::endl;
	std::cerr<<" parameter cannot be found in file"<<std::endl;
}

void XWINNMRacqu::wrongTypeErr() //parameter not found
{
	std::cerr<<std::endl<<"Error: XWINNMRacqu..wrong type"<<std::endl;
	std::cerr<<" the data element desired is NOT the same it is in the file"<<std::endl;
}


//simply sees if we can open the file
bool XWINNMRacqu::testopen(std::string fname)
{
	//first try to opent the file if it is not opened
	std::ifstream ffile(fname.c_str());
	if(ffile.fail()){
		ffile.close();
		return false;
	}else{
		ffile.close();
		return true;
	}


}

//opens the file for real, with an error message
// if it fails...it will 'kill' the old one
// first
bool XWINNMRacqu::open(std::string fname, bool showWarn)
{
	//first try to opent the file if it is not opened
	if(!file.is_open()){
		file.open(fname.c_str());
		if(file.fail()){
			if(showWarn) FileErr();
			return false;
		}else{
			return true;
		}
	}else{
		file.close();
		file.open(fname.c_str());
		if(file.fail()){
			if(showWarn) FileErr();
			return false;
		}else{
			return true;
		}
	}
	return true;
}

/*

##$SW_h= 1124.10071942446
##$TD= 8192

*/
//thus the parameter name is one line before the
// data element and the name is the first in thing
// in that line
bool XWINNMRacqu::findPar(std::string par, Vector<std::string> &parsed)
{
	//first try to opent the file if it is not opened
	if(!file.is_open()){
		file.open(fname.c_str());
		if(file.fail()){
			FileErr();
			return false;
		}
	}

	//first give the proper name to the finding entity
	std::string realPar="##$"+par;

	//reset the file pointer to the begining
	file.seekg(std::ios::beg);
	char liner[1024];
	while(!file.eof() && !file.fail())
	{
		file.getline(liner, 1024);
		parsed=parse_param(liner, '=');
		if(parsed.size()>0){
			if(parsed[0]==realPar){
				return true;
			}
		}
	}
	NotFoundErr();
	return false;
}

XWINNMRacqu::XWINNMRacqu(std::string fnamein)
{
	fname=fnamein;
}



//close the file upon destruction
XWINNMRacqu::~XWINNMRacqu()
{
	if(file.is_open())	file.close();
}



//because we do not what the data type
// for the data element is, we must write many of these
// this function and assume the user knows what the data type is
//this function gets a parameter called 'par'
// from the file, and asigned it to the data type
// 'data'.  it returns false if not found and true if found
// it also returns false if it knows the data types do not match

//the basic data element for a single number is
/*
##$DIGMOD= 1
*/

//we test for a 'string' type
// whcih bruker leaves as a '<' '>' item
bool XWINNMRacqu::get(std::string par, int &data)
{
	data=0; //set the par to some default

	Vector<std::string> parsed;

	//make sure the parameter is there
	if(findPar(par, parsed)){
		//get the same line that 'findPar' found
		//parsed should already have what we want

		if(parsed.size()>1){
			if(parsed[1].find('<') >= parsed[1].size()){
				data=std::atoi(parsed[1].c_str());
			}else{
				wrongTypeErr();
				return false;
			}
			return true;
		}
	}else{
		return false;
	}
    return false;
}

//the basic data element for a single number is
/*
##$O1= 2198.28676756606
*/

//behaves like the 'int' getter
bool XWINNMRacqu::get(std::string par, double &data)
{
	data=0; //set the par to some default

	Vector<std::string> parsed;

	//make sure the parameter is there
	if(findPar(par, parsed)){
		//get the same line that 'findPar' found
		//parsed should already have what we want

		if(parsed.size()>1){
			//we can test for the 'single-parameter' atribute
			if(parsed[1].find('<') >= parsed[1].size()){
				data=std::atof(parsed[1].c_str());
			}else{
				wrongTypeErr();
				return false;
			}
			return true;
		}
	}else{
		return false;
	}
    return false;
}

//the basic data element for a list of items is
/*
##$L= (0..31)
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
*/

//the first line tells us its name (an other varian internal unknown stuff)
//the second line tells us there is '3' parameters and then the  value of them
bool XWINNMRacqu::get(std::string par, Vector<double> &data)
{
	data.resize(0); //set the par to some default

	Vector<std::string> parsed;
	char liner[1024];

	//make sure the parameter is there
	if(findPar(par, parsed))
	{
		//get the next line after the parameter is found
		file.getline(liner, 1024);
		parsed=parse_param(liner);

		//test for string type
		data.resize(0);
		if(parsed.size()>=1){
			if(parsed[0].find('<')<=parsed[0].size()){
				wrongTypeErr();
				return false;
			}
			//check to make sure the second item is
			// a std::string (it will be enclosed in '<')
			for(int kk=0;kk<parsed.size();++kk){
				if(!parsed[kk].empty()){
					data.push_back(std::atof(parsed[kk].c_str()));
				}
			}
			return true;
		}
	}else{
		return false;
	}

    return false;
}
//the basic data element for a list of items is
/*
##$NUC1= <1H>
*/
//  The data element is ENCLOSED in < >
bool XWINNMRacqu::get(std::string par, std::string &data)
{
	data=""; //set the par to some default

	Vector<std::string> parsed;

	//make sure the parameter is there
	if(findPar(par, parsed)){
		//get the next line after the parameter is found

		//should be 2 items a '1 {#}'
		if(parsed.size()>1){
			//we can test for the 'single-parameter' atribute
			if(parsed[1].find('<')<parsed[1].size()){
				//remove the quotes
				data=parsed[1].substr(1);
				data=data.substr(0,data.size()-1);
			}else{
				wrongTypeErr();
				return false;
			}
			return true;
		}
	}else{
		return false;
	}

    return false;
}

/* THE STREAM *****************/
//The binary reader..to function properly
// it needs the ENTIRE directory
// thus is needs the acqus and/or acqus2


XWINNMRstream::XWINNMRstream(std::string DIR)
{	open(DIR); }

bool XWINNMRstream::open(std::string DIR, bool showWarn)
{
	//first we test if we can get the
	// 1D acqu file...if we cannot, then we are screwed
	isxwinnmr_=false;
	dir_=DIR;
	if(!acq1D.open(DIR+"/acqu",showWarn)){
		if(showWarn){
			std::cerr<<"Error: XWINNMRstream::open(directory)"<<std::endl;
			std::cerr<<" the 1D 'acqus' file cannot be opened"<<std::endl;
		}
		return false;
	}

	//now test to see if we have a a 2D file 'acqu2'
	if(!acq2D.testopen(DIR+"/acqu2")){
		is2D_=false;
	}else{
		acq2D.open(DIR+"/acqu2",showWarn);
		is2D_=true;
	}

	//now look for the 'fid' file if 1D
	// and 'ser' file if 2D

	if(is2D_){
		if(fp.is_open()) fp.close();
		std::string na=DIR+"/ser";
		fp.open(na.c_str(), std::ios::binary |  std::ios::in);
		if(fp.fail()){
			if(showWarn){
				std::cerr<<"Error: XWINNMRstream::open(directory)"<<std::endl;
				std::cerr<<" the 1D 'fid' file cannot be opened"<<std::endl;
			}
			return false;
		}

	}else{
		if(fp.is_open()) fp.close();
		std::string na=DIR+"/fid";
		fp.open(na.c_str(), std::ios::binary |  std::ios::in);
		if(fp.fail()){
			if(showWarn){
				std::cerr<<"Error: XWINNMRstream::open(directory)"<<std::endl;
				std::cerr<<" the 2D 'ser' file cannot be opened"<<std::endl;
			}
			return false;
		}
	}
	isxwinnmr_=true;
	return true;
}

void XWINNMRstream::close()
{	if(fp.is_open()) fp.close();	}


XWINNMRstream::~XWINNMRstream()
{	close();	}







END_BL_NAMESPACE


#endif




