/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 01-05-02
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
	A VNMR (from Varian) file reader and writer
	reads and writes VNMR 1D and 2D binary files

*/

#ifndef VNMR_STREAM_CC_
#define VNMR_STREAM_CC_ 1

#include "container/containers.h"
#include "utils/vnmrstream.h"
#include "utils/endians.h"
#include <iostream>
#include <algorithm>

BEGIN_BL_NAMESPACE

using namespace BL_endians;


/***** Finds parameters in a generic varian 'procpar' file  ***/

//this class gets a data chunk from
// the Varian 'proppar' file
//the procpar file contains the nessesary parameters
// for a the varian GUI to function properly
// it also contains the data (like sweep widths) and other
// data parameters NOT in the binary fid file

void VarianProcPar::FileErr() //a file read error
{
	std::cerr<<std::endl<<"Error: VarianProcPar..File read error"<<std::endl;
	std::cerr<<" file cannot be opened or read"<<std::endl;
}

void VarianProcPar::NotFoundErr() //parameter not found
{
	std::cerr<<std::endl<<"Error: VarianProcPar..Parameter Not Found"<<std::endl;
	std::cerr<<" parameter cannot be found in file"<<std::endl;
}

void VarianProcPar::wrongTypeErr() //parameter not found
{
	std::cerr<<std::endl<<"Error: VarianProcPar..wrong type"<<std::endl;
	std::cerr<<" the data element desired is NOT the same it is in the file"<<std::endl;
}
//sets the current position in the
// file to that of the desired parameter
//returns true if it finds the parameter
//a parameter in the procpar looks like

/*
sw1 1 1 2000000 0 0 2 1 0 1 64
1 1388.35861303
0
*/

//thus the parameter name is one line before the
// data element and the name is the first in thing
// in that line
bool VarianProcPar::findPar(std::string par)
{
	//first try to opent the file if it is not opened
	if(!file.is_open()){
		file.open(fname.c_str());
		if(file.fail()){
			FileErr();
			return false;
		}
	}

	//reset the file pointer to the begining
	file.seekg(std::ios::beg);
	Vector<std::string> parsed;
	char liner[1024];
	while(!file.eof() && !file.fail())
	{
		file.getline(liner, 1024);
		parsed=parse_param(liner);
		if(parsed.size()>0){
			if(parsed[0]==par){
				return true;
			}
		}
	}
	NotFoundErr();
	return false;
}

VarianProcPar::VarianProcPar(std::string fnamein)
{
	fname=fnamein;
	open(fnamein);

}

bool VarianProcPar::open(std::string fnamein, bool showWarn)
{

	if(file.is_open()) file.close();
	file.open(fnamein.c_str());
	if(file.fail()){
		if(showWarn) FileErr();
		return false;
	}
	return true;
}




//close the file upon destruction
VarianProcPar::~VarianProcPar()
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
nv 1 1 2000000 0 0 2 1 0 1 64
1 128
0
*/

//the first line tells us its name (an other varian internal unknown stuff)
//the second line tells us there is '1' parameter and then its value
//the thrird line tells us that this parameter is 'over'
bool VarianProcPar::get(std::string par, int &data)
{
	data=0; //set the par to some default

	Vector<std::string> parsed;
	char liner[1024];

	//make sure the parameter is there
	if(findPar(par)){
		//get the next line after the parameter is found
		file.getline(liner, 1024);
		parsed=parse_param(liner);

		//should be 2 items a '1 {#}'
		if(parsed.size()>1){
			//we can test for the 'single-parameter' atribute
			if(parsed[0]=="1"){
				if(parsed[1].find('\"') >= parsed[1].size()){
					data=std::atoi(parsed[1].c_str());
				}else{
					wrongTypeErr();
					return false;
				}
				return true;
			}
		}
	}else{
		return false;
    }
    return false;
}

//the basic data element for a single number is
/*
sw1 1 1 2000000 0 0 2 1 0 1 64
1 1388.35861303
0
*/

//the first line tells us its name (an other varian internal unknown stuff)
//the second line tells us there is '1' parameter and then its value
//the thrird line tells us that this parameter is 'over'
bool VarianProcPar::get(std::string par, double &data)
{
	data=0; //set the par to some default

	Vector<std::string> parsed;
	char liner[1024];

	//make sure the parameter is there
	if(findPar(par)){
		//get the next line after the parameter is found
		file.getline(liner, 1024);
		parsed=parse_param(liner);

		//should be 2 items a '1 {#}'
		if(parsed.size()>1){
			//we can test for the 'single-parameter' atribute
			if(parsed[0]=="1"){
				if(parsed[1].find('\"') >= parsed[1].size()){
					data=std::atof(parsed[1].c_str());
				}else{
					wrongTypeErr();
					return false;
				}
				return true;
			}
		}
	}else{
		return false;
    }
    return false;
}

//the basic data element for a list of items is
/*
sw1 1 1 2000000 0 0 2 1 0 1 64
3 1388.35861303 23434.9 45
0
*/

//the first line tells us its name (an other varian internal unknown stuff)
//the second line tells us there is '3' parameters and then the  value of them
bool VarianProcPar::get(std::string par, Vector<double> &data)
{
	data.resize(0); //set the par to some default

	Vector<std::string> parsed;
	char liner[1024];
	bool gotone=false;
	int maxct=0;
	int curct=1;

	//make sure the parameter is there
	if(findPar(par))
	{
		//while(!file.eof() && !file.fail())
		//{
		//get the next line after the parameter is found
			file.getline(liner, 1024);
			parsed=parse_param(liner);

			//if(!gotone)
			//{
				//should be 2 items a '1 {#}'
				if(parsed.size()>1){
					//we can test for the 'single-parameter' atribute
					maxct=atoi(parsed[0]);
					//check to make sure the second item is NOT
					// a std::string (it will be enclosed in '"')
					while(curct<=maxct){
						if(parsed[curct].find('\"')>=parsed[curct].size()){
							data.push_back(std::atof(parsed[curct].c_str()));
							curct++;
							gotone=true;
						}else{
							wrongTypeErr();
							return false;
						}
					}
				}
		/*	}else{
				if(maxct>curct)
				{
					//check to make sure the second item is NOT
					// a std::string (it will be enclosed in '"')
					if(parsed[0].find('\"')>=parsed[0].size() ){
						data.push_back(atof(parsed[0].c_str()));
						curct++;
					}else{
						wrongTypeErr();
						return false;
					}
				}else{
					break;
				}
			}
			*/
		//}

	}else{
		return false;
	}

    return false;
}
//the basic data element for a list of items is
/*
sw1 1 1 2000000 0 0 2 1 0 1 64
1 "1388.35861303"
0
*/
//the first line tells us its name (an other varian internal unknown stuff)
//the second line tells us there is '1' parameters and then the  value
//  The data element is ENCLOSED in QUOTES
bool VarianProcPar::get(std::string par, std::string &data)
{
	data=""; //set the par to some default

	Vector<std::string> parsed;
	char liner[1024];

	//make sure the parameter is there
	if(findPar(par)){
		//get the next line after the parameter is found
		file.getline(liner, 1024);
		parsed=parse_param(liner);

		//should be 2 items a '1 {#}'
		if(parsed.size()>1){
			//we can test for the 'single-parameter' atribute
			if(parsed[0]=="1"){
				if(parsed[1].find('\"')<parsed[1].size()){
					//remove the quotes
					data=parsed[1].substr(1);
					data=data.substr(0,data.size()-1);
				}else{
					wrongTypeErr();
					return false;
				}
				return true;
			}
		}
	}else{
		return false;
	}

    return false;
}

/***** The Binary Header Pareser for the Varian Binary fid File ***/

Head::Head()
{
  nblocks=0;
  ntraces=1;
  np=0;
  ebytes=sizeof(float);
  tbytes=0;
  nbheaders=0;
  bbytes=0;
  vers_id=0;
  status=S_DATA | S_32 | S_COMPLEX | S_NP | S_NI;
}

void Head::set(int ni,int np)
{
  nblocks=ni;
  ntraces=1;
  np=np*2;
  ebytes=sizeof(float);
  tbytes=np*ebytes;
  nbheaders=1;
  bbytes=ntraces*tbytes+nbheaders*sizeof(struct BlockHead);
  vers_id=65;
  status=S_DATA | S_32 | S_COMPLEX | S_NP | S_NI;
}

void Head::ByteSwap()
{
	int bige=BL_endians::AreWeBigEndian();
	if(!bige) BL_endians::ByteSwap(np);
	if(!bige) BL_endians::ByteSwap(ntraces);
	if(!bige) BL_endians::ByteSwap(nblocks);
	if(!bige) BL_endians::ByteSwap(ebytes);
	if(!bige) BL_endians::ByteSwap(tbytes);
	if(!bige) BL_endians::ByteSwap(bbytes);
	if(!bige) BL_endians::ByteSwap(vers_id);
	if(!bige) BL_endians::ByteSwap(status);
	if(!bige) BL_endians::ByteSwap(nbheaders);
}

std::ostream &operator<<(std::ostream &oo, Head &out)
{
	oo<<"VNMR Header Info: "<<std::endl;
	oo<<"Blocks: "<<out.nblocks<<std::endl
	  <<"Traces: "<<out.ntraces<<std::endl
	  <<"NumPoints*2: "<<out.np<<std::endl
	  <<"eBytes: "<<out.ebytes<<std::endl
	  <<"tBytes: "<<out.tbytes<<std::endl
	  <<"bBytes: "<<out.bbytes<<std::endl
	  <<"Version ID: "<<out.vers_id<<std::endl
	  <<"status: "<<out.status<<std::endl
	  <<"nbheaders: "<<out.nbheaders<<std::endl;
	return oo;
}

BlockHead::BlockHead()
{
  scale=0;

  status=S_DATA | S_32 | S_COMPLEX | MORE_BLOCKS |
            NP_CMPLX | NI_CMPLX ;

  index=0;
  mode=NP_PHMODE | NI_PHMODE;
  ctcount=1024;
  lpval=0;
  rpval=0;
  lvl=0;
  tlt=0;
}

void BlockHead::set()
{
  scale=0;

  status=S_DATA | S_32 | S_COMPLEX | S_FLOAT | MORE_BLOCKS |
            NP_CMPLX | NI_CMPLX;

  index=0;
  mode=NP_PHMODE | NI_PHMODE;
  ctcount=1024;
  lpval=0;
  rpval=0;
  lvl=0;
  tlt=0;
}

void BlockHead::ByteSwap()
{
	int bige=BL_endians::AreWeBigEndian();
	if(!bige) BL_endians::ByteSwap(scale);
	if(!bige) BL_endians::ByteSwap(status);
	if(!bige) BL_endians::ByteSwap(index);
	if(!bige) BL_endians::ByteSwap(mode);
	if(!bige) BL_endians::ByteSwap(ctcount);
	if(!bige) BL_endians::ByteSwap(lpval);
	if(!bige) BL_endians::ByteSwap(rpval);
	if(!bige) BL_endians::ByteSwap(lvl);
	if(!bige) BL_endians::ByteSwap(tlt);
}

std::ostream &operator<<(std::ostream &oo, BlockHead &out)
{
	oo<<"VNMR Block Header Info: "<<std::endl;
	oo<<"scale: "<<out.scale<<std::endl
	  <<"Status: "<<out.status<<std::endl
	  <<"Index: "<<out.index<<std::endl
	  <<"mode: "<<out.mode<<std::endl
	  <<"ctcount: "<<out.ctcount<<std::endl
	  <<"lpval: "<<out.lpval<<std::endl
	  <<"rpval: "<<out.rpval<<std::endl
	  <<"lvl: "<<out.lvl<<std::endl
	  <<"tlt: "<<out.tlt<<std::endl;
	return oo;
}



/****** VNMR stream for reading in the BINARY FIDs ****/

VNMRstream::VNMRstream(std::string in, std::ios::openmode inmode):
	iomode(inmode), fname(in), canwr(false)
{
	fp.open(in.c_str(), iomode);
	canwr=initFileHead();
}

VNMRstream::VNMRstream(const VNMRstream &cp):
	filehead(cp.filehead), blockhead(blockhead),
	iomode(cp.iomode),fname(cp.fname), canwr(cp.canwr)
{
	fp.open(fname.c_str(), iomode);
	canwr=initFileHead();
}

VNMRstream::~VNMRstream()
{ fp.close();	}


bool VNMRstream::initFileHead(bool showWarn)
{
	if(fp.fail())
	{
		if(showWarn)
			std::cerr<<std::endl<<"Error: VNMRstream.. cannot open file "<<fname<<std::endl;
		return false;
	}

	if(!fp.fail())
	{
		if(iomode & std::ios::in)
		{
			fp.read((char *)&filehead, sizeof(Head));
			filehead.ByteSwap();
			/*if(filehead.ntraces!=1)
			{
				std::cerr<<std::endl<<"Error: VNMRstream "<<std::endl;
				std::cerr<<" Can only process 1D files... "<<std::endl;
				std::cerr<<" Error: there are "<<filehead.ntraces<<" fids ... "<<std::endl;
				return false;
			}

			if (filehead.ebytes != 4)
			{
				std::cerr<<std::endl<<"Error: VNMRstream "<<std::endl;
				std::cerr<<"Error: this is only for 32bit data \n"
						 <<" (file header element bytes = "<<filehead.ebytes<<")"<<std::endl;
				return false;
			}*/
		}else if(iomode & std::ios::out){
			filehead.set(1,1024);
			blockhead.set();
			return true;
		}
	}else{
		return false;
	}
	return true;
}

//is the file a VNMR file?
bool VNMRstream::isVNMR()
{
	if(canwr){
		return (filehead.status & S_COMPLEX);
	}
	return false;
}

bool VNMRstream::open(std::string in, std::ios::openmode mode, bool showWarn)
{
	fname=in;
	iomode=mode;
	fp.open(in.c_str(), iomode);
	canwr=initFileHead(showWarn);
	return canwr;
}

bool VNMRstream::is2D()
{	return (filehead.nblocks>1);	}

void VNMRstream::close()
{	fp.close();	}

Vector<complex> VNMRstream::readone()
{	return readone(true);	}

Vector<complex> VNMRstream::readone(bool readhead)
{
	int bige=BL_endians::AreWeBigEndian();
	if(readhead){
		if(!fp.read((char *)&blockhead,sizeof(BlockHead)))
		{
			std::cerr<<"Error: unable to read datablock header"<<std::endl;
			canwr= false;
		}
	}

	char *data; data=new char[filehead.tbytes];

	fp.read(data,filehead.np*filehead.ebytes);
	int bytes=fp.gcount();

	if (bytes != filehead.tbytes) {
		std::cerr<<"Error: unable to read "<<filehead.tbytes
				 <<" bytes (only read "<<bytes<<" bytes)"<<std::endl;
	   canwr= false;
	}

	int i, j;

	bool isfloat=(blockhead.status & S_FLOAT) == S_FLOAT;
	float tmf;
	int NumOPoints=filehead.np;
	Vector<complex> listonums(NumOPoints/2, 0.0);

	if(isfloat){
		for (j=0,i=0;i<NumOPoints;i += 2,++j)
		{
			tmf=((float*)data)[i];
			if(!bige) ByteSwap(tmf);
			listonums[j].Re(tmf);
			tmf=((float*)data)[i+1];
			if(!bige) ByteSwap(tmf);
			listonums[j].Im( -tmf);
		}
	}else{
		long tmi;
		for (j=0,i=0;i<NumOPoints;i += 2,++j)
		{
			tmi=((int*)data)[i];
			if(!bige) ByteSwap(tmi);
			listonums[j].Re(tmi);
			tmi=((int*)data)[i+1];
			if(!bige) ByteSwap(tmi);
			listonums[j].Im(-tmi);
		}
	}
	delete [] data;
	return listonums;
}


bool VNMRstream::read(Vector<complex> &listonums, int index)
{
	if(is_open() && iomode & std::ios::in){
		fp.seekg( 0, std::ios::beg); //reset the file pointer
		canwr=initFileHead(); //read the file head

		if(canwr)
		{
			for(int i=0;i<index-1;++i){
				if(filehead.nbheaders>i){
					if(fp.seekg(sizeof(BlockHead), std::ios::cur)){
						std::cerr<<"Error: VNMRstream.read: unable to skip to index "<<index
							<<" only "<<i<<" in the file "<<std::endl;
						return false;
					}
				}
				if(fp.seekg(filehead.np*filehead.ebytes, std::ios::cur)){
					std::cerr<<"Error: VNMRstream.read: unable to skip to index "<<index
							<<" only "<<i<<" in the file "<<std::endl;
					return false;
				}
			}
			listonums=readone();

		}else{
			std::cerr<<std::endl<<"Error: VNMRstream.read(Vector<complex>) "<<std::endl;
			std::cerr<<" file cannot be read...."<<std::endl;
			return false;
		}
		return true;
	}else{
		std::cerr<<std::endl<<"Error: VNMRstream.read(Vector<complex>) "<<std::endl;
		std::cerr<<" file is not yet opened or io mode is not 'std::ios::in', cannot read...."<<std::endl;
		return false;
	}

}


bool VNMRstream::read(matrix &twoD)
{
	if(is_open() && iomode & std::ios::in){
		fp.seekg( 0, std::ios::beg); //reset the file pointer
		canwr=initFileHead(); //read the file head
		twoD.resize( filehead.np/2,filehead.nblocks, 0.0);

		Vector<complex> tmpv;

		for(int colct=0;colct<filehead.nblocks; ++colct)
		{
			tmpv=readone();
			if(canwr){
				twoD.putCol(colct, tmpv);
			}else{
				std::cerr<<std::endl<<"Error: VNMRstream.read(matrix) "<<std::endl;
				std::cerr<<" file cannot be read...."<<std::endl;
				return false;
			}
		}
		return true;
	}else{
		std::cerr<<std::endl<<"Error: VNMRstream.read(matrix) "<<std::endl;
		std::cerr<<" file is not yet opened or io mode is not 'std::ios::in', cannot read...."<<std::endl;
		return false;
	}
}


bool VNMRstream::writeone(const Vector<complex> &out)
{
	if(canwr){
		int NumOPoints=out.size()*2;
		float *tmf=new float[NumOPoints];
		fp.write((char *)&blockhead, sizeof(BlockHead));

		for (int j=0,i=0;i<NumOPoints;i += 2,++j)
		{
			tmf[i]=float(out[j].Re());
			tmf[i+1]=float(out[j].Im());
		}
		if(is_open() & iomode & std::ios::out){
			fp.write((char *)tmf, sizeof(float)*NumOPoints);
			return true;
		}else{
			return false;
		}
	}else{
		return false;
	}
}

bool VNMRstream::write(Vector<complex> &oneD)
{
	if(is_open() && iomode & std::ios::out){
		fp.seekg( 0, std::ios::beg); //reset the file pointer
		filehead.set(1, oneD.size());
		blockhead.set();
		fp.write((char *)&filehead, sizeof(Head));
		if(!writeone(oneD)){
			std::cerr<<std::endl<<"Error: VNMRstream.write(Vector<complex>) "<<std::endl;
			std::cerr<<" cannot write to file...."<<std::endl;
			return false;
		}
		return true;
	}else{
		std::cerr<<std::endl<<"Error: VNMRstream.write(Vector<complex>) "<<std::endl;
		std::cerr<<" file is not yet opened or io mode is not 'std::ios::out', cannot write...."<<std::endl;
		return false;
	}
}



bool VNMRstream::write(matrix &twoD)
{
	if(is_open() && iomode & std::ios::out){
		fp.seekg( 0, std::ios::beg); //reset the file pointer
		filehead.set(twoD.rows(), twoD.cols());
		blockhead.set();
		fp.write((char *)&filehead, sizeof(Head));
		for(int i=0; i<twoD.cols();++i){
			if(!writeone(twoD.col(i))){
				std::cerr<<std::endl<<"Error: VNMRstream.write(matrix) "<<std::endl;
				std::cerr<<" cannot write to file...."<<std::endl;
				return false;
			}
		}
		return true;
	}else{
		std::cerr<<std::endl<<"Error: VNMRstream.write(matrix) "<<std::endl;
		std::cerr<<" file is not yet opened or io mode is not 'std::ios::out', cannot write...."<<std::endl;
		return false;
	}
}

bool VNMRstream::write(complex *data, int d1)
{
	Vector<complex> dat(data, d1);
	return write(dat);
}

bool VNMRstream::write(complex **data, int d1, int d2)
{
	matrix dat(d1, d2, data);
	return write(dat);
}


VNMRstream &operator<<(VNMRstream &out, Vector<complex> &in)
{
	out.write(in);
	return out;
}

VNMRstream &operator<<(VNMRstream &out, matrix &in)
{
	out.write(in);
	return out;
}

VNMRstream &operator>>(VNMRstream &out, Vector<complex> &in)
{
	out.read(in);
	return out;
}

VNMRstream &operator>>(VNMRstream &out, matrix &in)
{
	out.read(in);
	return out;
}

/************ The TOTAL Varian Data Reader ******/


bool VNMRdata::open(std::string dir, bool showWarn)
{
	if(!VNMRstream::open(dir+"/fid", std::ios::in, showWarn)) return false;
	if(!VarianProcPar::open(dir+"/procpar",showWarn)) return false;
	return true;
}


void VNMRdata::getType()
{
	get("pss", pss);
	get("nv", nv);
	get("np", np);
	get("nf", nf);
	get("nv2", nv2);

	int curp=fp.tellg();
	fp.seekg(std::ios::beg);
	if(!canwr) canwr=initFileHead();
	fp.seekg(curp);

	dx = np / 2;
	dy = nv;
	dz = pss.size();
	dt = filehead.ntraces * filehead.nblocks / (dz * nv);

	if( nf == nv2)
	{
		theType=CSI;
		std::cout<<"CSI data set"<<std::endl;
		dx = np / 2;
		dy = nv;
		dz = nv2;
		dt = filehead.ntraces * filehead.nblocks / (dz * dy);
	}

	if ( pss.size()>1 && (nf == pss.size()*nv))
	{
		theType=Compressed;
		std::cout<<"Compressed - Compressed Multislice data set"<<std::endl;
		dx = np / 2;
		dy = nv;
		dz = pss.size();
		dt = 1;
	}

	if ( pss.size()==1 && (nv2 == filehead.nblocks ))
	{
		theType=ThreeD;
		std::cout<<"3D image set"<<std::endl;
		dx = np / 2;
		dy = nv;
		dz = nv2;
		dt = 1;
	}

	if ( pss[0]==0 && (nv == 0 ))
	{
		theType=TwoD;
		std::cout<<"2D spectra set"<<std::endl;
		dx = np / 2;
		dy = filehead.nblocks;
		dz = 1;
		dt = 1;
	}
}

void VNMRdata::readPetable(std::string fname)
{
	if(fname.size()>0){
		std::string petable= petablePath+"/"+fname;
		std::ifstream table(petable.c_str());
		if(table.fail()){
			std::cerr<<"Error: VNMRdata.read"<<std::endl;
			std::cerr<<" cannot find the 'petable' "<<petable<<std::endl;
			return;
		}
		//read in the file...We Must skip the first line
		char liner[256];

		int skline=1, ctline=0;
		while(!table.eof() && !table.fail()){
			table.getline(liner, 256);
			if(ctline>=skline){
				peTable.push_back(petableData(std::atoi(liner), ctline-1));
			}
			ctline++;
		}
		std::sort(peTable.data(), peTable.data()+peTable.size());
	}
}

void VNMRdata::read(Vector<matrix> &dat)
{
	//go to begining of the file
	fp.seekg(std::ios::beg);

	//set the chunk sizes
	getType();


	if(dt>1){
		std::cerr<<std::endl<<" Error: Cannot handle CSI type data yet"<<std::endl;
		return;
	}

	dat.resize(dz, matrix(dx,dy,0.0));
	int i, j,k=0;
	Vector<complex> tmpv;

	//int bige=BL_endians::AreWeBigEndian();
	std::cout<<"Data Size is: "<<dx<<" x "<<dy<<" x "<<dz<<std::endl;
	std::cout<<filehead<<blockhead<<std::endl;

	for(i=0;i<dz;++i){
		for(j=0;j<dy;++j){

			switch(theType){
				case ThreeD:
				case TwoD:
				case Compressed:
					tmpv=readone(k<filehead.nbheaders);
					if(canwr){
						dat[i].putRow(j,tmpv);
						++k;
					}else{
						std::cout<<"could not read the "<<i<<" "<<j<<" data rows"<<std::endl;
						return;
					}
			}
		}
	}

	//must correct for Varians weird ordering of the
	//fids...the 'new ordering' is stored in a file
	//called 'petable' whos name is in the procpar
	// file HOWEVER, this file is not with the normal
	// data set, so you need to find it...
	std::string petable;
	get("petable",petable);
	readPetable(petable);
	if(petable.size()>0 && peTable.size()==0){
		std::cerr<<"Error: VNMRdata.readPetable"<<std::endl;
		std::cerr<<" error reading the petable "<<petable<<std::endl;
		return;
	}

	Vector<matrix> olddat;
	olddat=dat;
	if(petable.size()>0)
	{
		for(i=0;i<peTable.size();++i){
			for(j=0;j<olddat.size();++j){
				for(k=0;k<dy;++k)
				dat[j](i,k)= olddat[j](peTable[i].index, k);
			}
		}
	}

}

END_BL_NAMESPACE



#endif
