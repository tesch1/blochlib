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

/***** Finds parameters in a generic Bruker 'acqus' file  ***/

//sets the current position in the
// file to that of the desired parameter
//returns true if it finds the parameter
//a parameter in the procpar looks like

/*

##$SW_h= 1124.10071942446
##$TD= 8192

*/

#ifndef XWINNMR_STREAM_H_
#define XWINNMR_STREAM_H_ 1



#include "container/containers.h"
#include "utils/endians.h"
#include <iostream>

BEGIN_BL_NAMESPACE

using namespace BL_endians; //for the endian conversions

class XWINNMRacqu{
	friend class XWINNMRstream;
	private:
		std::ifstream file; //this is the text file to read

		void FileErr(); //a file read error
		void NotFoundErr(); //parameter not found
		void wrongTypeErr(); //wrong parameter type

		//sets the current position in the
		// file to that of the desired parameter
		bool findPar(std::string par, Vector<std::string> &parsed);

		bool testopen(std::string fname);

	public:
		std::string fname; //the file name

		XWINNMRacqu(): fname(""){}
		XWINNMRacqu(std::string fname);

		~XWINNMRacqu();

		bool open(std::string fname, bool showWar=true);

		//because we do not what the data type
		// for the data element is, we must write many of these
		// this function and assume the user knows what the data type is
		//this function gets a parameter called 'par'
		// from the file, and asigned it to the data type
		// 'data'.  it returns false if not found and true if found
		// it also returns false if it knows the data types do not match

		bool get(std::string par, int &data);
		bool get(std::string par, double &data);
		bool get(std::string par, Vector<double> &data);
		bool get(std::string par, std::string &data);

};


//The binary reader..to function properly
// it needs the ENTIRE directory
// thus is needs the acqus and/or acqus2
class XWINNMRstream{
	friend class XWINNMRacqu;
	private:

		std::fstream fp; //FID data
		std::string dir_;

		bool isxwinnmr_; //is the dir an XWINNMR

		bool is2D_; //determined at the opening operations

	public:

		XWINNMRacqu acq1D; //the 1D params
		XWINNMRacqu acq2D; //the 2D params (may not exsist)

		XWINNMRstream():isxwinnmr_(false),is2D_(false){}
		XWINNMRstream(std::string DIR);

		~XWINNMRstream();

	//opend the 'directory' tree
	//returns false if fail OR not an XWINNMR file
		bool open(std::string DIR, bool showWarn=true);

	//is the file an XWINNMR directory
		inline bool isXWINNMR(){	return isxwinnmr_;	}

	//has the file been opend?
		inline bool is_open(){	return fp.is_open();	}

	//is the data file a 2D or 1D file
		inline bool is2D(){	return is2D_;	}

	//close the fstream
		void close();

	//read a 1D data block from the file
		template<class Ctype>
		bool read(Vector<Complex<Ctype> > &oned);

	//read a 2D data chunk (it will read 1D files with matrix dims 1xNpts)
		template<class Ctype>
		bool read(_matrix<Complex<Ctype>, FullMatrix > &twoD);


};


//read a 1D data block from the file
template<class Ctype>
bool XWINNMRstream::read(Vector<Complex<Ctype> > &oned)
{
	if(is2D_){
		std::cerr<<"Error: XWINNMRstream::read(Vector)"<<std::endl;
		std::cerr<<" data format is 2D!! not 1D"<<std::endl;
		return false;
	}
	//now we need to find the 'TD' points in
	// the acqu file
	int size1D;
	acq1D.get("TD", size1D);
	if(size1D<=0){
		std::cerr<<"Error: XWINNMRstream::read(Vector)"<<std::endl;
		std::cerr<<" Not points in 'fid' file..."<<std::endl;
		return false;
	}

	//now read all the pts in the FID
	oned.resize(size1D/2);
	int tmI=0;
	int ct=0;
    int bige=BL_endians::AreWeBigEndian();

	while(ct<size1D/2){
		if(fp.eof())
		{
			std::cerr<<"Error: XWINNMRstream::read(Vector)"<<std::endl;
			std::cerr<<" Not ENOUGH points in 'fid' file..."<<std::endl;
			return false;
		}
	//the first one is a real number
		fp.read((char *)&tmI, sizeof(int));
		if(!bige) BL_endians::ByteSwap(tmI); //bruker is Big endian...
		oned[ct].Re(tmI);

		if(fp.eof())
		{
			std::cerr<<"Error: XWINNMRstream::read(Vector)"<<std::endl;
			std::cerr<<" Not ENOUGH points in 'fid' file..."<<std::endl;
			return false;
		}

	//the first one is a imag number
		fp.read((char *)&tmI, sizeof(int));
		if(!bige) BL_endians::ByteSwap(tmI); //bruker is Big endian...
		oned[ct].Im(tmI);
		++ct;
	}
	return true;
}


//read a 2D data chunk (it will read 1D files with matrix dims 1xNpts)
template<class Ctype>
bool XWINNMRstream::read(_matrix<Complex<Ctype>, FullMatrix > &twoD)
{
	if(!is2D_){
		std::cerr<<"Error: XWINNMRstream::read(matrix)"<<std::endl;
		std::cerr<<" data format is 2D!! not 1D"<<std::endl;
		return false;
	}
	//now we need to find the 'TD' points in
	// the acqu file
	int size1D;
	acq1D.get("TD", size1D);
	if(size1D<=0){
		std::cerr<<"Error: XWINNMRstream::read(matrix)"<<std::endl;
		std::cerr<<" No points in 'ser' file..."<<std::endl;
		return false;
	}

	//now we need to find the 'TD' points in
	// the acqu file 2D
	int size2D;
	acq2D.get("TD", size2D);
	if(size2D<=0){
		std::cerr<<"Error: XWINNMRstream::read(matrix)"<<std::endl;
		std::cerr<<" No points in 'ser' file..."<<std::endl;
		return false;
	}

	//now read all the pts in the FID
	twoD.resize(size2D, size1D/2);
	int tmI=0;
	int bige=BL_endians::AreWeBigEndian();

	for(int i=0;i<size2D;++i){
		for(int j=0;j<size1D/2;++j){
			if(fp.eof())
			{
				std::cerr<<"Error: XWINNMRstream::read(matrix)"<<std::endl;
				std::cerr<<" Not ENOUGH points in 'ser' file..."<<std::endl;
				return false;
			}
		//the first one is a real number
			fp.read((char *)&tmI, sizeof(int));
			if(!bige) BL_endians::ByteSwap(tmI); //bruker is Big endian...
			twoD(i,j).Re(tmI);

			if(fp.eof())
			{
				std::cerr<<"Error: XWINNMRstream::read(matrix)"<<std::endl;
				std::cerr<<" Not ENOUGH points in 'ser' file..."<<std::endl;
				return false;
			}

		//the first one is a imag number
			fp.read((char *)&tmI, sizeof(int));
			if(!bige) BL_endians::ByteSwap(tmI); //bruker is Big endian...
			twoD(i,j).Im(tmI);
		}
	}
	return true;
}


END_BL_NAMESPACE


#endif




