/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10-04-01
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

THe matlab file writer
	VERSION 5 ONLY...this is simple class that
	simply 'write' a matlab 5 binary file FOR a COMPLEX matrix
	sound easy??  HA, i do not think matlab could make it a
	bigger pain in the arse...



I'd like to thank to the makers of Gamma

  S.A. Smith, T.O. Levante, B.H. Meier and R.R. Ernst
  Computer Simulations in Magnetic Resonance:
  An Object-Oriented Programming Approach
  J. Magn. Reson., 106a, 75-105, (1994S, Series A)

http://gamma.magnet.fsu.edu/

for much of this setup...

*/

#ifndef _matlab5_h_
#define _matlab5_h_ 1

#include "utils/endians.h"
#include "container/matrix/matrix.h"
#include "container/rankType.h"
#include<string>
#include<fstream>

#ifndef ON_WINDOWS
 #include "blochconfig.h"
#endif

#if HAVE_UNISTD_H
	#include <unistd.h>
#endif

#include <time.h>


BEGIN_BL_NAMESPACE


#define VARTOSTRING(IN) std::string(#IN)

// Matlab 5 element types
extern const int miINT8;
extern const int miUINT8;
extern const int miINT16;
extern const int miUINT16;
extern const int miINT32;
extern const int miUINT32;
extern const int miSINGLE;
extern const int miDOUBLE;
extern const int miINT64;
extern const int miUINT64;
extern const int miMATRIX;

extern const int miCELL;
extern const int miSTRUCT;
extern const int miOBJ;
extern const int miCHARARR;
extern const int miSPARSEARR;
extern const int miDOUBLEARR;
extern const int miSINGLEARR;
extern const int miINT8ARR;
extern const int miUINT8ARR;
extern const int miINT16ARR;
extern const int miUINT16ARR;
extern const int miINT32ARR;
extern const int miUINT32ARR;


class Matlab5Header {
	public:
		static const int hdrsize=128;
		char txt[124];
		short version;
		char mm;
		char ii;
		int BigEndian_; //the 'input file' endian'ness'

		Matlab5Header();
		Matlab5Header(const Matlab5Header &in);

		inline int size()const	{	return hdrsize;	}

		bool SetHeader();

		inline void SetEndian()    { (BL_endians::AreWeBigEndian())?SetBigEnd():SetLittleEnd(); }
		inline void SetLittleEnd() { mm='I'; ii='M'; }
		inline void SetBigEnd()    { mm='M'; ii='I'; }
		inline int  BigEndian()    { if(mm=='M') return 1; return 0; }

		int write(std::ofstream &out);
		int read(std::ifstream &in);
};

class Matlab5Tag {
	public:
		char chars[8];
		int nbytes;
		int dtype;
		bool iscompressed;
		int bigendian;

		Matlab5Tag();
		Matlab5Tag(int bigend);
		Matlab5Tag(const Matlab5Tag &in);
		Matlab5Tag(int Type, int nbyts, int bigend=0, int cmpres=0);
		Matlab5Tag &operator=(const Matlab5Tag &in);

		~Matlab5Tag(){}

		inline int size() const {	return 8;	}
		inline int bytes() const	{	return nbytes;	}
		inline int Bytes()	const	{	return nbytes;	}

		int read(std::ifstream &in);
		int write(std::ofstream &out);

		std::string DataType() const;
		std::string DataSymbol() const;

		void IsCompressed();
};

class Matlab5SubData
{
	public:
		Matlab5Tag tag;
		char *data;
		int bigendian;

		Matlab5SubData();
		Matlab5SubData(int bigend);
		~Matlab5SubData();

		int read(std::ifstream &in);
		int write(std::ofstream &out);
};


class Matlab5ArrayFlag : public Matlab5SubData{
	public:
		int firstb;
		int Class;


		Matlab5ArrayFlag();
		Matlab5ArrayFlag(int bigend);
		Matlab5ArrayFlag(int cmplx, int Clas, int bigend, int global=0,int logical=0);

		template<class numt, class structure>
		Matlab5ArrayFlag(const _matrix<numt, structure>& mx, int bigend=0) :
			Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miUINT32,8,bigend);
			if(BiggerType(numt, complex) ) firstb = 8;
			else	    firstb = 0;
			Class     = miDOUBLEARR;
			data    = new char[8];
			data[2] = char(firstb);
			data[3] = char(Class);
		}

		template<class structure>
		Matlab5ArrayFlag(const _matrix<char, structure>& mx, int bigend=0) :
			Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miUINT32,8,bigend);
			firstb = 0;
			Class     = miCHARARR;
			data    = new char[8];
			data[2] = char(firstb);
			data[3] = char(Class);
		}

		inline int size()const{	return 2*tag.size();	}
		inline bool IsComplex(){ return firstb>=8?true:false;	}
		std::string ClassType() const;

		int write(std::ofstream &out);
		int read(std::ifstream &in);
};

class Matlab5ArrayName : public Matlab5SubData{
	public:

		std::string name;
		int numchars;

		Matlab5ArrayName();
		Matlab5ArrayName(int bigend);
		Matlab5ArrayName(const std::string &inn, int bigend=0);

		int size() const;

		int size(std::string inn) const;

		int write(std::ofstream &out);
		int read(std::ifstream &in);
};

class Matlab5ArrayDims : public Matlab5SubData{
	public:

		int *dims;
		int ndims;

		Matlab5ArrayDims();
		Matlab5ArrayDims(int bigend);

		template<class numt, class structure>
		Matlab5ArrayDims(const _matrix<numt, structure> &in, int bigend=0):
			Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miINT32, 8, bigend);
			dims= new int[2];
			dims[0]=in.rows();
			dims[1]=in.cols();
			ndims=2;
		}

		~Matlab5ArrayDims();

		int size() const;

		template<class numt, class structure>
		int size(const _matrix<numt, structure> &mx) const
		{			return 2*sizeof(long)+tag.size();	}

		template<class numt>
		int size(const Vector<numt> &mx) const
		{			return 2*sizeof(long)+tag.size();	}


		int write(std::ofstream &out);
		int read(std::ifstream &in);
};


class Matlab5RealData : public Matlab5SubData
{
	public:

		Matlab5RealData();
		Matlab5RealData(int bigend);

		template<class numt, class structure>
		Matlab5RealData(const _matrix<numt, structure> &in, int bigend=0) :
		 	Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miDOUBLE, in.rows()*in.cols(), bigend);
		}


		template<class structure>
		Matlab5RealData(const _matrix<char, structure> &in, int bigend=0) :
		 	Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miUINT8, in.rows()*in.cols(), bigend);
		}

		inline int size() const {	return tag.size();	}

		template<class numt, class structure>
		int size(const _matrix<numt, structure> &in) const
		{	return sizeof(double)*in.rows()*in.cols()+tag.size();	}

		template<class structure>
		int size(const _matrix<char, structure> &in) const
		{	return sizeof(char)*in.rows()*in.cols()+tag.size();	}




//		~Matlab5RealData();

		template<class numt, class structure>
		int write(const _matrix<numt, structure> &in, std::ofstream &fp)
		{
			long type = miDOUBLE;
			fp.write((char*)&type, sizeof(long));
			type = in.rows()*in.cols()*sizeof(double);
			fp.write((char*)&type, sizeof(long));
			int i,j;
			double d;
			for(j=0; j<in.cols(); j++)
			{
				for(i=0; i<in.rows(); i++)
				{
					d = Re(in(i,j));
					fp.write((char*)&d, sizeof(double));
				}
			}
			return 1;
		}


		template<class structure>
		int write(const _matrix<char, structure> &in, std::ofstream &fp)
		{
			long type = miUINT8;
			fp.write((char*)&type, sizeof(long));
			type = in.rows()*in.cols()*sizeof(char);
			fp.write((char*)&type, sizeof(long));
			int i,j;
			char d;
			for(j=0; j<in.cols(); j++)
			{
				for(i=0; i<in.rows(); i++)
				{
					d = in(i,j);
					fp.write((char*)&d, sizeof(char));
				}
			}
			return 1;
		}

		//a 'dummy' reader
		//template<class numt, class structure>
		int read(std::ifstream &fp)
		{
			if(!tag.read(fp)) return 0;
			if(!tag.iscompressed)	fp.seekg(tag.nbytes, std::ios::cur);
			return 1;
		}

		template<class numt, class structure>
		int read(const _matrix<numt, structure> &in, std::ifstream &fp)
		{
			if(!tag.read(fp)) return 0;
			int nr = in.rows();
			int nc = in.cols();;
			if(!BiggerType(numt, double))
			{
				std::cerr<<std::endl<<"Error: Matlab5ImagData.read()"<<std::endl;
				std::cerr<<" attempted to read 'double' data into a non double matrix..."<<std::endl;
				std::cerr<<" cannot read.."<<std::endl;
				return 0;
			}
			int swap = 0;
			if(tag.bigendian != BL_endians::AreWeBigEndian()) swap++;
			double d;
			char c;
			int i, j, ct=0;
			for(j=0; j<nc; j++){
				for(i=0; i<nr; i++){
					switch(tag.dtype)
					{
						case miUINT8:  fp.read(&c,sizeof(char)); d=double(c); ct++; break;
						case miDOUBLE:
						default:
							fp.read((char*)&d,sizeof(double));
							if(swap) ByteSwap(d);
							break;
					}
					in.put(i,j,in.get(i,j)+d);
				}
			}
		//	cout<<endl<<"CT:"<<ct<<"CT%8: "<<8-ct%8<<endl;
			if(ct){
				if(ct%8 !=0){
					for(int i=0;i<(8-ct%8);i++){	fp.read((char *)&c, sizeof(char));	}
				}
			}
		//	cout<<"MATRIX: "<<in<<endl;
			return 1;
		}


};


class Matlab5ImagData : public Matlab5SubData
{
	public:

		Matlab5ImagData();
		Matlab5ImagData(int bigend);

		template<class numt, class structure>
		Matlab5ImagData(const _matrix<numt, structure> &in, int bigend=0) :
		 	Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miDOUBLE, in.rows()*in.cols(), bigend);
		}

		template<class structure>
		Matlab5ImagData(const _matrix<char, structure> &in, int bigend=0) :
		 	Matlab5SubData(bigend)
		{
			tag=Matlab5Tag(miINT8, in.rows()*in.cols(), bigend);
		}


//		~Matlab5ImagData();

		inline int size() const {	return tag.size();	}

		template<class numt, class structure>
		int size(const _matrix<numt, structure> &in) const
		{	return sizeof(double)*in.rows()*in.cols()+tag.size();	}

		template<class structure>
		int size(const _matrix<char, structure> &in) const
		{	return sizeof(char)*in.rows()*in.cols()+tag.size();	}


		template<class numt, class structure>
		int write(const _matrix<numt, structure> &in, std::ofstream &fp)
		{
			long type = miDOUBLE;
			fp.write((char*)&type, sizeof(long));
			type = in.rows()*in.cols()*sizeof(double);
			fp.write((char*)&type, sizeof(long));
			int i,j;
			double d;
			for(j=0; j<in.cols(); j++)
			{
				for(i=0; i<in.rows(); i++)
				{
					d = Im(in(i,j));
					fp.write((char*)&d, sizeof(double));
				}
			}
			return 1;
		}


		//a 'dummy' reader
		template<class numt, class structure>
		int read(std::ifstream &fp)
		{
			if(!tag.read(fp)) return 0;
			if(!tag.iscompressed)	fp.seekg(tag.nbytes, std::ios::cur);
			return 1;
		}

		template<class numt, class structure>
		int read(const _matrix<numt, structure> &in, std::ifstream &fp)
		{
			if(!tag.read(fp)) return 0;
			int nr = in.rows();
			int nc = in.cols();;
			if(!BiggerType(numt, complex))
			{
				std::cerr<<std::endl<<"Error: Matlab5ImagData.read()"<<std::endl;
				std::cerr<<" attempted to read complex data into a non complex vector..."<<std::endl;
				std::cerr<<" cannot read.."<<std::endl;
				return 0;
			}
			int swap = 0;
			if(tag.bigendian != BL_endians::AreWeBigEndian()) swap++;
			double d;
			char c;
			int i, j,ct=0;
			for(j=0; j<nc; j++){
				for(i=0; i<nr; i++){
					switch(tag.dtype)
					{
						case miUINT8:  fp.read(&c,sizeof(char)); d=double(c); ct++; break;
						case miDOUBLE:
						default:
							fp.read((char*)&d,sizeof(double));
							if(swap) ByteSwap(d);
							break;
					}
					in.put(i,j,in.get(i,j)+complex(0,d));
				}
			}
			if(ct){
				if(ct%8 !=0){
					for(int i=0;i<(8-ct%8);i++){	fp.read((char *)&c, sizeof(char));	}
				}
			}
			return 1;
		}


};


/************** A TOTAL Matlab matrix 'chunk' ***********/

class Matlab5MatData {
	public:
		Matlab5Tag tag;
		Matlab5ArrayFlag  arrayflag;
		Matlab5ArrayDims  arraydim;
		Matlab5ArrayName  arrayname;
		Matlab5RealData  rdata;
		Matlab5ImagData  idata;
		int bigendian;

		Matlab5MatData();
		Matlab5MatData(int bigend);

		template<class numt, class structure>
		Matlab5MatData(const _matrix<numt, structure> &mx,std::string name, int bigend=0) :
			tag(miMATRIX, size(mx,name)), arrayflag(mx, bigend),
			arraydim(mx), arrayname(name), rdata(mx), idata(mx),
			bigendian(bigend)
		{}

		template<class numt, class structure>
		int size(const _matrix<numt, structure>& mx, const std::string& name) const
		{
			int nbts = arrayflag.size();
			nbts+= arraydim.size(mx);
			nbts += arrayname.size(name);
			nbts += rdata.size(mx);
			if(BiggerType(numt, complex))	nbts += idata.size(mx);
			return nbts;
		}


		int size() const;

		template<class numt, class structure>
		int write(std::ofstream& fp,const _matrix<numt, structure>& mx,    const std::string& N="fid")
		{
			int err = tag.write(fp);

			Matlab5ArrayFlag tma(mx, bigendian);
			err *= tma.write(fp);

			Matlab5ArrayDims tmd(mx, bigendian);
			err *= tmd.write(fp);

			Matlab5ArrayName tmn(N,bigendian);
			err *= tmn.write(fp);

			err *= rdata.write(mx,fp);
			if(tma.IsComplex()) err *= idata.write(mx,fp);
			return err;
		}

		template<class numt, class structure>
		int read(std::ifstream& fp, _matrix<numt, structure>& mx)
		{
			//int pos=fp.tellp();
			tag.bigendian=bigendian;
			int err  = tag.read(fp);
			arrayflag.bigendian=bigendian;
			err *= arrayflag.read( fp);
			arraydim.bigendian=bigendian;
			err *= arraydim.read( fp);
			if(arraydim.ndims>=2)	mx.resize(arraydim.dims[0], arraydim.dims[1]);
			mx.fill(ZeroType<numt>::zero());
			arrayname.bigendian=bigendian;
			err *= arrayname.read(fp);
			//cout<<"err: "<<err<<" DAT: "<<arrayname.name<<endl<<mx<<endl;
			rdata.bigendian=bigendian;
			err *= rdata.read(mx,fp);
			idata.bigendian=bigendian;
		//	cout<<"read: "<<arrayflag.IsComplex()<<endl;
			if(arrayflag.firstb){ err *= idata.read(mx,fp);	}
			if(!err)
			{
				std::cerr<<std::endl<<"Error: Matlab5MatData(read)"<<std::endl;
				std::cerr<<" error in reading Matrix in..."<<std::endl;
				return (0);
			}
			//fp.seekp(pos, ios::cur);
			return err;

		}


		void skip(std::ifstream& fp);
		void whos(std::ostream& ostr, std::ifstream& fp);

};


class Matlab5{

	private:
		std::string fname;
		std::ifstream ifp;
		std::ofstream ofp;
		std::ios::openmode iomode;
		int fsize;
		bool canwrt;
		Matlab5Header outhdr;


	public:

		Matlab5():
			iomode(std::ios::binary | std::ios::out), fsize(0),canwrt(false)
		{}

		Matlab5(std::string fname, std::ios::openmode mode=(std::ios::binary | std::ios::out));

		bool SetFile(std::string fname, std::ios::openmode mode=std::ios::binary | std::ios::out);
		bool ModeCheck( std::ios::openmode imode);
		int FileSize(std::ofstream &fp);
		int FileSize(std::ifstream &fp);

		void whos(std::ostream &out=std::cout);

		bool open(std::string fname, std::ios::openmode mode=std::ios::binary | std::ios::out);
		inline bool is_open(){	return ifp.is_open()||ofp.is_open();	}

		void close();

		//this returns a general variable (Always in a 'highest' bit count for that container)
		//you must 'cast' it down to the proper format
		complex GetVar(const std::string &vname);
		matrix GetMatrix(const std::string  &name);
		Vector<complex> GetVector(const std::string &name);
		std::string GetString(const std::string &name);

		template<class NumT>
		void get(const std::string &name, NumT &in)
		{
			complex out=GetVar(name);
			if(!BiggerType(NumT, complex) && Im(out) !=0.0){
				std::cerr<<std::endl<<"Warning: Matlab5.get(name, var)"<<std::endl;
				std::cerr<<" you desire a non complex number, yet the number retrieved"<<std::endl;
				std::cerr<<" is a complex number...takeing the real part."<<std::endl;
				std::cerr<<" for variable '"<<name<<"'"<<std::endl;
			}
			in=NumT(Re(out));
		}
//getters

		void get(const std::string &name, complex &in);
		void get(const std::string &name, std::string &in);

		template<class NumT>
		void get(const std::string &name, Vector<NumT> &in)
		{
			in=Re(GetVector(name));
		}
		void get(const std::string &name, Vector<complex> &in);

		template<class NumT, class Structure>
		void get(const std::string &name, _matrix<NumT, Structure> &in)
		{
			in=Re(GetMatrix(name));
		}

		template<class Structure>
		void get(const std::string &name, _matrix<complex, Structure> &in)
		{
			 in=GetMatrix(name);
		}

//writers
		void initwrite(std::ofstream &out);

		template<class NumT>
		void put(std::string name,NumT out)
		{
			if(canwrt){
				initwrite(ofp);
				_matrix<NumT, FullMatrix> mx(1,1,out);
				Matlab5MatData mdat(mx,name, BL_endians::AreWeBigEndian());
				if(!mdat.write(ofp, mx,name))
				{
					std::cerr<<"Warning:: Matlab5::put(complex)"<<std::endl;
					std::cerr<<" Error in writing to file....."<<std::endl;
				}
			}else{
				std::cerr<<"Warning:: Matlab5::put(complex)"<<std::endl;
				std::cerr<<" file Cannot be written '"<<name<<"' to (check 'ios' flags an permissions)"<<std::endl;
			}
		}

		template<class NumT>
		void put(std::string name, const Vector<NumT> &mx )
		{
			if(mx.size()==0)
			{
				std::cerr<<"Warning:: Matlab5::put(Vector)"<<std::endl;
				std::cerr<<" Empty Vector  '"<<name<<std::endl;
			}
			if(canwrt){
				initwrite(ofp);
				_matrix<NumT, FullMatrix> mxo(mx.size(), 1);
				mxo.putCol(0,mx);
				Matlab5MatData mdat(mxo,name, BL_endians::AreWeBigEndian());
				if(!mdat.write(ofp, mxo,name))
				{
					std::cerr<<"Warning:: Matlab5::put(Vector)"<<std::endl;
					std::cerr<<" Error in writing '"<<name<<"' to file....."<<std::endl;
				}
			}else{
				std::cerr<<"Warning:: Matlab5::put(Vector)"<<std::endl;
				std::cerr<<" file Cannot be written to (check 'ios' flags and/or permissions)"<<std::endl;
			}
		}

		template<class NumT, class structure>
		void put( std::string name,const _matrix<NumT, structure> &mx)
		{
			if(mx.rows()==0||mx.cols()==0)
			{
				std::cerr<<"Warning:: Matlab5::put(Matrix)"<<std::endl;
				std::cerr<<" Empty Matrix '"<<name<<"' "<<std::endl;
			//	return;
			}

			if(canwrt){
				initwrite(ofp);
				Matlab5MatData mdat(mx,name, BL_endians::AreWeBigEndian());
				if(!mdat.write(ofp, mx,name))
				{
					std::cerr<<"Warning:: Matlab5::put(Matrix)"<<std::endl;
					std::cerr<<" Error in writing '"<<name<<"' to file....."<<std::endl;
				}
			}else{
				std::cerr<<"Warning:: Matlab5::put(matrix)"<<std::endl;
				std::cerr<<" file Cannot be written to (check 'ios' flags an permissions)"<<std::endl;
			}
		}

		void put(std::string name, const std::string &mx );
		void put(std::string name,const char *mx );

		template<class Numt, int N>
		void put(std::string name, const Vector<coord<Numt, N> > &mx)
		{
			Vector<Numt> tmz(mx.size(), ZeroType<Numt>::zero());
			for(int i=0;i<N;++i)
			{
				std::string outname=name+itost(i);
				for(int j=0;j<mx.size();j++) tmz(j)=mx(j)(i);
				put(outname, tmz);
			}
		}

		template<class Numt, int N>
		void put(std::string name, const coord<Numt, N> &mx)
		{
			double tmz;
			for(int i=0;i<N;++i)
			{
				std::string outname=name+itost(i);
				tmz=mx(i);
				put(outname, tmz);
			}
		}

		template<class Numt, class Structure, int N>
		void put(std::string name, const _matrix<coord<Numt, N>, Structure> &mx)
		{
			_matrix<Numt, Structure> tmz(mx.rows(),mx.cols(), ZeroType<Numt>::zero());
			for(int i=0;i<N;++i)
			{
				std::string outname=name+itost(i);
				for(int j=0;j<mx.rows();j++) {
					for(int k=0;k<mx.cols();k++)
					{
						tmz(j,k)=mx(j,k)(i);
					}
				}
				put(outname, tmz);
			}
		}
};

typedef Matlab5 matstream;

END_BL_NAMESPACE


#endif
