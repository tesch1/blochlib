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

#ifndef VNMR_STREAM_H_
#define VNMR_STREAM_H_


#include "container/containers.h"
#include "utils/endians.h"
#include <iostream>

BEGIN_BL_NAMESPACE


/* Stuff specific to Varian Data file */
#define S_FLOAT 0x8
#define S_COMPLEX	0x10
#define S_DATA 0x1
#define S_32 0x4
#define S_NP 0x800
#define S_NI 0x2000
#define NI_PHMODE	0x100
#define NP_PHMODE	0x1
#define NP_CMPLX 0x100
#define NI_CMPLX 0x400
#define MORE_BLOCKS	0x80

struct Head
{
 	long nblocks,ntraces,np,ebytes,tbytes,bbytes;
 	short vers_id,status;
 	long nbheaders;

	Head();
 	void set(int ni, int np);
	void ByteSwap();
};

std::ostream &operator<<(std::ostream &oo, Head &out);


struct BlockHead
{
 	short scale,status,index,mode;
 	long ctcount;
 	float lpval,rpval,lvl,tlt;

	BlockHead();
 	void set();
	void ByteSwap();
};
std::ostream &operator<<(std::ostream &oo, BlockHead &out);


//this class gets a data chunk from
// the Varian 'proppar' file
//the procpar file contains the nessesary parameters
// for a the varian GUI to function properly
// it also contains the data (like sweep widths) and other
// data parameters NOT in the binary fid file

class VarianProcPar{
	private:
		std::ifstream file; //this is the text file to read

		void FileErr(); //a file read error
		void NotFoundErr(); //parameter not found
		void wrongTypeErr(); //wrong parameter type

		//sets the current position in the
		// file to that of the desired parameter
		bool findPar(std::string par);

	public:
		std::string fname; //the file name

		VarianProcPar(): fname(""){}
		VarianProcPar(std::string fname);

		bool open(std::string ifname, bool showWarn=true);

		~VarianProcPar();

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

class VNMRstream{
	friend class VNMRdata;
	private:
		Head filehead;
		BlockHead blockhead;

		std::fstream fp;
		std::ios::openmode iomode;
		std::string fname;

		Vector<complex> readone();
		Vector<complex> readone(bool readhead);
		bool writeone(const Vector<complex> &out);

		bool initFileHead(bool showWarm=true);
		bool canwr;

	public:

		VNMRstream(){}
		VNMRstream(std::string in, std::ios::openmode inmode=std::ios::binary | std::ios::in);
		VNMRstream(const VNMRstream &cp);

		~VNMRstream();

	//open the file for reading or writing and read the header
	//of the file
		bool open(std::string fname, std::ios::openmode mode=std::ios::binary | std::ios::in, bool showWarn=true);

	//has the file been opend?
		inline bool is_open(){	return fp.is_open();	}

	//is the file a varian file?
		bool isVNMR();

	//is the data file a 2D or 1D file
		bool is2D();

	//close the fstream
		void close();

	//read a 1D data block from the file
		bool read(Vector<complex> &oned, int index=1);

	//read a 2D data chunk (it will read 1D files with matrix dims 1xNpts)
		bool read(matrix &twoD);

	//write a VNMR formated 1D file foma Vector
		bool write(Vector<complex> &oned);

	//write a VNMR formated 2D file froma matrix
		bool write(matrix &twoD);

	//write VNMR formated 1D file from a complex array
		bool write(complex *oneD, int d1);

	//write a VNMR formated 2D file from a complex array
		bool write(complex **twoD, int d1, int d2);

		friend VNMRstream &operator<<(VNMRstream &out, Vector<complex> &in);
		friend VNMRstream &operator<<(VNMRstream &out, matrix &in);

		friend VNMRstream &operator>>(VNMRstream &out, Vector<complex> &in);
		friend VNMRstream &operator>>(VNMRstream &out, matrix &in);

};


//This class encompasses BOTH the 'VNMRstream' and the 'VarianProcPar'
//to allow the reading of an entire dirctory strcuture

class VNMRdata:
 public VNMRstream,
 public VarianProcPar

{
	private:

	//data elements associated with determination
	//of the type of data we have 3D, 2D, compressed

		Vector<double> pss;
		int nv;
		int nv2;
		int np;
		int nf;

	//internal vars for reading in the
	// the data
		int dx, dy, dz, dt;

		enum dataType{CSI, Compressed, ThreeD,TwoD, OneD};

		dataType theType;

	//gets the
		void getType();
		void readPetable(std::string fname);

		class petableData{
			public:
				int which;
				int index;
				petableData(): which(0), index(0){}
				petableData(int w, int l): which(w), index(l){}

				friend inline bool operator<(const petableData &inw1,const petableData &inw)
				{	return inw1.which<inw.which;}
				friend inline bool operator>(const petableData &inw1,const petableData &inw)
				{	return inw1.which>inw.which;}
				friend inline bool operator<=(const petableData &inw1,const petableData &inw)
				{	return inw1.which<=inw.which;}
				friend inline bool operator>=(const petableData &inw1,const petableData &inw)
				{	return inw1.which>=inw.which;}
				friend inline bool operator==(const petableData &inw1,const petableData &inw)
				{	return inw1.which!=inw.which;}
				friend inline bool operator!=(const petableData &inw1,const petableData &inw)
				{	return inw1.which==inw.which;}
		};

		Vector<petableData> peTable;

	public:

		std::string petablePath;

		VNMRdata(){}
		VNMRdata(std::string dir):
			VNMRstream(dir+"/fid"),
			VarianProcPar(dir+"/procpar"),
			pss(0), nv(0), nv2(0), np(0), nf(0),
			theType(OneD),petablePath(".")
		{}

		bool open(std::string dir, bool showWarn=true);

		void read(Vector<matrix> &data);
	//read a 1D data block from the file
		inline bool read(Vector<complex> &oned, int index=1)
		{	return VNMRstream::read(oned, index);	}

	//read a 2D data chunk (it will read 1D files with matrix dims 1xNpts)
		inline bool read(matrix &twoD)
		{	return VNMRstream::read(twoD);	}

};



END_BL_NAMESPACE


#endif

