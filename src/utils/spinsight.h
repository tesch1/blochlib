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

#ifndef SPINSIGHT_STREAM_H_
#define SPINSIGHT_STREAM_H_ 1



#include "container/containers.h"
#include "utils/endians.h"
#include <iostream>

BEGIN_BL_NAMESPACE

using namespace BL_endians; //for the endian conversions

class SpinSightacq{
	friend class SpinSightStream;
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

		SpinSightacq(): fname(""){}
		SpinSightacq(std::string fname);

		~SpinSightacq();

		bool open(std::string fname, bool showWarn=true);
		inline bool is_open(){	return file.is_open();	}

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
class SpinSightStream{
	private:

		std::fstream fp; //FID data
		std::string dir_;

		bool isSpinsight; //is the file a spinsight?
		bool is2D_; //determined at the opening operations

	public:

		SpinSightacq acq1D; //the 1D params
		SpinSightacq acq2D; //the 2D params (may not exsist)

		SpinSightStream():isSpinsight(false),is2D_(false){}
		SpinSightStream(std::string DIR);

		~SpinSightStream();

	//opend the 'directory' tree
		bool open(std::string DIR, bool showWarn=true);

	//is the dir a spinsight dir
		inline bool isSpinSight(){	return isSpinsight;	}

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
		bool read(_matrix<Complex<Ctype>, FullMatrix> &twoD);

};

//read a 1D data block from the file
template<class Ctype>
bool SpinSightStream::read(Vector<Complex<Ctype> > &oned)
{


	if(is2D_>1){
		std::cerr<<"Warning: SpinSightStream::read(Vector)"<<std::endl;
		std::cerr<<" data format is 2D!! not 1D"<<std::endl;
		std::cerr<<" not reading anything..."<<std::endl;
		return false;
	}
	//now we need to find the 'TD' points in
	// the acq file
	int size1D;
	bool got1D=acq1D.get("al", size1D);
	if(got1D && size1D<=0){
		std::cerr<<"Error: SpinSightStream::read(Vector)"<<std::endl;
		std::cerr<<" Not points in 'data' file..."<<std::endl;
		return false;
	}

	// get length of file:
	fp.seekg (0, std::ios::end);
	int NumOPoints = fp.tellg()/sizeof(int);
	fp.seekg (0, std::ios::beg);

	if(got1D){
		if(NumOPoints>size1D*2 || NumOPoints<size1D*2){
			std::cerr<<"Error: SpinSightStream  number of points in file: "<<NumOPoints/2<<std::endl;
			std::cerr<<" is not comparable with the number of points"<<std::endl;
			std::cerr<<" given by the 1D dimensions: "<<size1D<<std::endl;
			std::cerr<<" Not getting any data..."<<std::endl;
			return false;
		}
	}else{
		size1D=NumOPoints/2;
	}

	//now read all the pts in the FID
	oned.resize(size1D);
	int i=0, ct=0;
	long tmpint;
	//real bits
	int bige=BL_endians::AreWeBigEndian();
	while(i<size1D && !fp.eof()){
		fp.read((char *)&tmpint, sizeof(int));
		if(!bige) ByteSwap(tmpint); //spinsight is Big endian...
		oned[ct].Re(tmpint);
		++i; ++ct;
	}
	//imaginary bits
	ct=0;
	while(i<2*size1D && !fp.eof()){
		fp.read((char *)&tmpint, sizeof(int));
		if(!bige) ByteSwap(tmpint); //spinsight is Big endian...
		oned[ct].Im(tmpint);
		++i; ++ct;
	}


	return true;
}




//read a 2D data chunk (it will read 1D files with matrix dims 1xNpts)
// it will read a 1D as well, but with only one column in the
// matrix
template<class Ctype>
bool SpinSightStream::read(_matrix<Complex<Ctype>, FullMatrix> &twoD)
{
	//now we need to find the 'TD' points in
	// the acqu file
	int size1D;
	acq1D.get("al", size1D);
	if(size1D<=0){
		std::cerr<<"Error: SpinSightStream::read(matrix)"<<std::endl;
		std::cerr<<" No points in 'data' file..."<<std::endl;
		return false;
	}

	//now we need to find the 'TD' points in
	// the acqu file 2D
	int size2D=0;
	acq1D.get("al2", size2D);
	if(size2D<=0){
		std::cerr<<"Error: SpinSightStream::read(matrix)"<<std::endl;
		std::cerr<<" No points in 'data' file..."<<std::endl;
		return false;
	}

	// get length of file:
	fp.seekg (0, std::ios::end);
	int NumOPoints = fp.tellg()/sizeof(int);
	fp.seekg (0, std::ios::beg);

	if(NumOPoints>size1D*size2D*2 || NumOPoints<size2D*size1D*2){
		std::cerr<<"Error: SpinSightStream number of points in file: "<<NumOPoints<<std::endl;
		std::cerr<<" is not comparable with the number of points"<<std::endl;
		std::cerr<<" given by the 2D dimensions: "<<size1D<<" "<<size2D<<std::endl;
		std::cerr<<" Not getting any data..."<<std::endl;
		return false;
	}

	//now read all the pts in the FID
	twoD.resize(size2D,size1D);
	int tmpint=0,i=0;
	int bige=BL_endians::AreWeBigEndian();
	Vector<double> tmD(NumOPoints);
	while(i<NumOPoints && !fp.eof()){
		fp.read((char *)&tmpint, sizeof(int));
		if(!bige) ByteSwap(tmpint); //spinsight is Big endian...
		tmD[i]=tmpint;
		++i;
	}

	//now reorder the monster into out matrix
	//the input matrix is organized a Real block,
	// then the Imaginary block...
	//so the first FID in the data is from
	// Real(0..size2D)+Imag(size2D*size1D..size2D*size1D+size1D)
	//so the second FID in the data is from
	// Real(size2D..2*size2D)+Imag(size2D*size1D+size2D..size2D*size1D+2*size2D)
	for(i=0;i<size2D;++i){
		int reCT=i*size1D;
		int imCT=i*size1D + size2D*size1D;
		for(int j=0;j<size1D;++j){
			twoD(i, j) = complex(tmD(reCT), tmD(imCT));
			reCT++;imCT++;
		}
	}
	return true;
}


END_BL_NAMESPACE


#endif




