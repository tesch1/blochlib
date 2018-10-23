/* wavestream.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 10.20.02
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
 	wavestream.h-->writes and reads a WAV audio file...
 */

#ifndef _BL_WAVE_FILE_H_
#define _BL_WAVE_FILE_H_ 1

#include "container/Vector/Vector.h"
#include "container/matrix/matrix.h"
#include <stack>
#include <vector>

BEGIN_BL_NAMESPACE



//let this guy know about the entire Riff File
class RiffFile;

/** The Riff file format data chunk **/

class RiffData
{
	public:

		char name[5];

	// the length, read from the second chunk header entry
		long int size;

	//the type of data chunk
		char subType[5];

	//the start of the chunk
	// offset in the file (use via ifstream::tellg())
		long int start;  // the file offset in bytes of the chunk contents

	//the start of the 'next' chunk
	//i.e. the end of this chunk
		long int end;

	// initialize at the file's current read position, and mark the file as bad
	// if there's an error.
		RiffData():
			size(0), start(0), end(0)
		{};

		RiffData(RiffFile& file);

	//determin if another data chunk is before or
	// after this one
		inline bool operator < (const RiffData& other) const
		{ return start < other.start; };

	//determin if another data chunk is the same one as
	// this one
		inline bool operator == (const RiffData& other) const
		{
			return (std::strcmp(name,other.name)==0) &&
				   (size == other.size) &&
			       (std::strcmp(subType,other.subType)==0) &&
				   (start == other.start);
		}
};


/***
THe 'Riff' file...consisting of a set
of RiffData elements
**/
class RiffFile {

	private:
		std::fstream fp;

		long int totalSize;

	//use the 'stack' container
	// in the STL to hold all the various
	// RiffData chunks in last-in-first-out order
		std::stack<RiffData, std::vector<RiffData> > data;

	public:
		RiffFile(){}
		RiffFile(std::string fname);
		RiffFile(const char * fname);
		//RiffFile(const RiffFile &rhs);

		RiffFile& operator=(RiffFile &rhs);

	//the master opener
	// gets the header, and set up the file for reading
	// returns 'false' if failure
		bool open(const char *fname, bool showWarn=true);

	//the closer
		inline bool close(){	if(fp.is_open()) fp.close(); return true;	}

		~RiffFile();

	//is the stack empty?
		bool empty() const;

	//empty the current stack
	// and reset the file to the begining
		bool clear();

	//add us a new element to the stack
		inline bool push()
		{	return push(0);	}

	//add us a new element to the stack
		bool push(const char* cType );

	//remove the first element on the stack
		bool pop();

	//the size of the stack
		long int currentSize() const;

		const char* currentName() const;
		const char* currentSubType() const;

		bool getNextItem(std::string& type, std::string& value);

		inline std::fstream &file()
		{ return fp; };

		std::ostream &print(std::ostream &oo, bool expand=true, int indent=0);

	protected:

		bool readItem(std::string& type, std::string& value);
};

std::ostream &operator<<(std::ostream &oo, RiffFile &out);


/**** THE WAVE FILE which is a specail RiffFile...****/
//NOTE:: if the file is opened for writing...then the header
//to the wave file is written ON CLOSE (or Destruction)
// this is because the number of data points in the file
// needs to be know BEFORE the header can be read...
// As a result the 'write' commands, simply fill
// a Buffer, which is then written upon close as well
// The buffer is a vector of 'shorts' that are converted
// to 'chars' if bitsPerChannel==8
class wavestream
{
	private:
		int bufferLen_;
		char *buffer_;
		unsigned long BufferStep; //initial buffer lengths
		unsigned long numSteps; //number of buffer length steps
		bool delBuffer; //to delete our buffer?

	protected:
		RiffFile riffFile;

		unsigned short formatType;
		unsigned short numChannels;
		unsigned long sampleRate;
		unsigned long bytesPerSecond;
		unsigned short bytesPerSample;
		unsigned short bitsPerChannel;
		unsigned long dataLength;

		const char* error;
	// true if we already wrote the header
		bool wroteHeader;

		std::ios::openmode iomode;

	public:
		wavestream();
		wavestream(const char *fname, std::ios::openmode iom);
		wavestream( wavestream &rhs);
		wavestream &operator=( wavestream &rhs);

		~wavestream();

		bool open(const char* name, std::ios::openmode);
		bool reset();
		bool close();

		unsigned short getFormatType() const;
		void setFormatType(unsigned short type);

		bool isCompressed() const;

		unsigned short getNumChannels() const;
		void setNumChannels(unsigned short num);

		unsigned long getSampleRate() const;
		void setSampleRate(unsigned long rate);

		unsigned long getBytesPerSecond() const;
		void setBytesPerSecond(unsigned long bytes);

		unsigned short getBytesPerSample() const;
		void setBytesPerSample(unsigned short bytes);

		short int getBitsPerChannel() const;
		void setBitsPerChannel(unsigned short bits);

		unsigned long getNumSamples() const;
		void setNumSamples(unsigned long num);

		float getNumSeconds() const;

		long int getDataLength() const;
		void setDataLength(long int numBytes);

		bool formatMatches(const wavestream& other);

		void copyFormat(const wavestream& other);

		void setupFormat(int sampleRate = 44100, short bitsPerChannel = 16, short channels = 1);

		inline std::fstream &file()
		{ return riffFile.file(); }

		inline RiffFile &getRiffFile()
		{ return riffFile; }

//the header writer
		bool writeHeader();

//read the raw...
		bool readRaw(char* buffer, unsigned long numBytes = 1);
		bool writeRaw(char* buffer, unsigned long numBytes = 1);

//Single bit outs
		bool read(unsigned char& sample);
		bool read(char& sample);
		bool read(short& sample);
		bool read(int& sample);
		bool read(float& sample);
		bool read(double& sample);
		bool read(Complex<float>& sample);
		bool read(Complex<double>& sample);


	//large chunks out...
		bool read(unsigned char *sample, unsigned long ct=1);
		bool read(char *sample, unsigned long ct=1);
		bool read(short *sample, unsigned long ct=1);
		bool read(int *sample, unsigned long ct=1);
		bool read(float* sample, unsigned long ct=1);
		bool read(double* sample, unsigned long ct=1);

	//vector readers
		bool read(Vector<unsigned char> &sample);
		bool read(Vector<char> &sample);
		bool read(Vector<short> &sample);
		bool read(Vector<int> &sample);
		bool read(Vector<float> &sample);
		bool read(Vector<double> &sample);
		bool read(Vector<Complex<float> > &sample);
		bool read(Vector<Complex<double> > &sample);

	//MATRIX readers...these will either give a 2xN or 1xN depending
	// on the number of channels!!
		bool read(_matrix<unsigned char, FullMatrix> &sample);
		bool read(_matrix<char, FullMatrix> &sample);
		bool read(_matrix<short, FullMatrix> &sample);
		bool read(_matrix<int, FullMatrix> &sample);
		bool read(_matrix<float, FullMatrix> &sample);
		bool read(_matrix<double, FullMatrix> &sample);
		bool read(_matrix<Complex<float>, FullMatrix> &sample);
		bool read(_matrix<Complex<double>, FullMatrix> &sample);



	//normal data type writers
		bool write(unsigned char sample);
		bool write(char sample);
		bool write(short sample);
		bool write(int sample);
		bool write(float sample);
		bool write(double sample);
		bool write(Complex<float> sample);
		bool write(Complex<double> sample);

	//pointer writers
		bool write(unsigned char *sample, unsigned long ct=1);
		bool write(char *sample, unsigned long ct=1);
		bool write(short *sample, unsigned long ct=1);
		bool write(int *sample, unsigned long ct=1);
		bool write(float *sample, unsigned long ct=1);
		bool write(double *sample, unsigned long ct=1);

	//Vector writers...
		bool write(Vector<unsigned char> sample);
		bool write(Vector<char> sample);
		bool write(Vector<short> sample);
		bool write(Vector<int> sample);
		bool write(Vector<float> sample);
		bool write(Vector<double> sample);
		bool write(Vector<Complex<float> > sample);
		bool write(Vector<Complex<double> > sample);

	//MATRIX readers...these will either give a 2xN or 1xN depending
	// on the number of channels!!
		bool write(_matrix<unsigned char, FullMatrix> sample);
		bool write(_matrix<char, FullMatrix> sample);
		bool write(_matrix<short, FullMatrix> sample);
		bool write(_matrix<int, FullMatrix> sample);
		bool write(_matrix<float, FullMatrix> sample);
		bool write(_matrix<double, FullMatrix> sample);
		bool write(_matrix<Complex<float>, FullMatrix> sample);
		bool write(_matrix<Complex<double>, FullMatrix> sample);

//find out 'sections'
		bool getFirstItem(std::string& type, std::string& value);
		bool getNextItem(std::string& type, std::string& value);

//copy one Wave into another
		bool copyFrom(wavestream& other);

		inline const char* getError() const
		{ return error; };

		bool operator==(const wavestream &rhs);

		inline void clearError(){ error = 0; };

		std::ostream &print(std::ostream &oo);



};


std::ostream & operator<<(std::ostream & oo, wavestream &out);


/*** The out streams to write data elements to a file ***/

wavestream &operator<<(wavestream &oo, unsigned char out);
wavestream &operator<<(wavestream &oo, Vector<unsigned char> out);
wavestream &operator<<(wavestream &oo, _matrix<unsigned char, FullMatrix> out);
wavestream &operator<<(wavestream &oo, char out);
wavestream &operator<<(wavestream &oo, Vector<char> out);
wavestream &operator<<(wavestream &oo,  _matrix<char, FullMatrix> out);
wavestream &operator<<(wavestream &oo, short out);
wavestream &operator<<(wavestream &oo, Vector<short> out);
wavestream &operator<<(wavestream &oo,  _matrix<short, FullMatrix> out);
wavestream &operator<<(wavestream &oo, int out);
wavestream &operator<<(wavestream &oo, Vector<int> out);
wavestream &operator<<(wavestream &oo,  _matrix<int, FullMatrix> out);
wavestream &operator<<(wavestream &oo, float out);
wavestream &operator<<(wavestream &oo, Vector<float> out);
wavestream &operator<<(wavestream &oo,  _matrix<float, FullMatrix> out);
wavestream &operator<<(wavestream &oo, double out);
wavestream &operator<<(wavestream &oo, Vector<double> out);
wavestream &operator<<(wavestream &oo,  _matrix<double, FullMatrix> out);
wavestream &operator<<(wavestream &oo, Complex<float> out);
wavestream &operator<<(wavestream &oo, Vector<Complex<float> > out);
wavestream &operator<<(wavestream &oo,  _matrix<Complex<float> , FullMatrix> out);
wavestream &operator<<(wavestream &oo, Complex<double> out);
wavestream &operator<<(wavestream &oo, Vector<Complex<double> > out);
wavestream &operator<<(wavestream &oo,  _matrix<Complex<double> , FullMatrix> out);

wavestream &operator>>(wavestream &oo, char &out);
wavestream &operator>>(wavestream &oo, Vector<char> &out);
wavestream &operator>>(wavestream &oo, _matrix<char, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, unsigned char &out);
wavestream &operator>>(wavestream &oo, Vector<unsigned char> &out);
wavestream &operator>>(wavestream &oo, _matrix<unsigned char, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, short &out);
wavestream &operator>>(wavestream &oo, Vector<short> &out);
wavestream &operator>>(wavestream &oo, _matrix<short, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, int &out);
wavestream &operator>>(wavestream &oo, Vector<int> &out);
wavestream &operator>>(wavestream &oo, _matrix<int, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, float &out);
wavestream &operator>>(wavestream &oo, Vector<float> &out);
wavestream &operator>>(wavestream &oo, _matrix<float, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, double &out);
wavestream &operator>>(wavestream &oo, Vector<double> &out);
wavestream &operator>>(wavestream &oo, _matrix<double, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, Complex<float> &out);
wavestream &operator>>(wavestream &oo, Vector<Complex<float> > &out);
wavestream &operator>>(wavestream &oo, _matrix<Complex<float>, FullMatrix> &out);
wavestream &operator>>(wavestream &oo, Complex<double> &out);
wavestream &operator>>(wavestream &oo, Vector<Complex<double> > &out);
wavestream &operator>>(wavestream &oo, _matrix<Complex<double>, FullMatrix> &out);


END_BL_NAMESPACE

#endif
