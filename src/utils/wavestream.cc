
/* wavestream.cc ********/


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
 	wavestream.cc-->writes and reads a WAV audio file...
 */
#ifndef _BL_WAVE_FILE_CC_
#define _BL_WAVE_FILE_CC_ 1

#include "utils/endians.h"
#include "utils/wavestream.h"


BEGIN_BL_NAMESPACE


/***
some constant data lists
 named to avoid any silly confusion in other programs
**/

// These are the tags inside a Riff File
// and their translation in english
struct BL__RiffTag__ {
	std::string typeName;  // four-letter name
	std::string realName;  // English name
};

const int BL__numExtraTypes__ = 24;
const BL__RiffTag__ BL__RiffDataTypes[BL__numExtraTypes__] = {
	{ "DISP", "Display name" },
	{ "IARL", "Archival location" },
	{ "IART", "Artist" },
	{ "ICMS", "Commissioned" },
	{ "ICMT", "Comments" },
	{ "ICOP", "Copyright" },
	{ "ICRD", "Creation date" },
	{ "ICRP", "Cropped" },
	{ "IDIM", "Dimensions" },
	{ "IDPI", "Dots Per Inch" },
	{ "IENG", "Engineer" },
	{ "IGNR", "Genre" },
	{ "IKEY", "Keywords" },
	{ "ILGT", "Lightness" },
	{ "IMED", "Medium" },
	{ "INAM", "Name" },
	{ "IPLT", "Palette setting" },
	{ "IPRD", "Product" },
	{ "ISBJ", "Subject" },
	{ "ISFT", "Software" },
	{ "ISHP", "Sharpness" },
	{ "ISRC", "Source" },
	{ "ISRF", "Source Form" },
	{ "ITCH", "Technician" },
};

/***************RiffData******/

RiffData::RiffData(RiffFile& thefile)
{

//get the endian ness of the machine
	int bige=BL_endians::AreWeBigEndian();

// read the data name (4 chars)
	thefile.file().read(name, 4*sizeof(char));
	name[4] = '\0'; //set the termination

// read the data size
	thefile.file().read((char *)&size, sizeof(long));

// reverse the endianism of the chunk size.
// if the platform is NOT little endian
	if(bige) BL_endians::ByteSwap(size);

// if this is a RIFF or LIST chunk, read its subtype.
	if (std::strcmp(name, "RIFF") == 0
		|| std::strcmp(name, "LIST") == 0)
	{

	// read the data name (4 chars)
		thefile.file().read(subType, 4*sizeof(char));
		subType[4] = '\0';

	// subtract the subtype from the size of the data.
		size -= 4;
	} else{
		*subType = '\0';
	}

// the chunk starts end the name and size.
	start = thefile.file().tellg();

// the next chunk starts end this one,
	end = start + size;
//but starts on a word boundary.
	if (end % 2)	end++;
}

/***************RiffFile******/

RiffFile::RiffFile(std::string fname)
{	open(fname.c_str());	}

RiffFile::RiffFile(const char *fname)
{	open(fname);	}

bool RiffFile::open(const char *fname, bool showWarn)
{
	if(fp.is_open()) fp.close();

	fp.open(fname, std::ios::in | std::ios::binary);
	if (fp.fail() || !clear()) {
		fp.close();

		if(showWarn)
		{
			std::cerr<<std::endl<<"Error: RiffFile(std::string) "<<std::endl;
			std::cerr<<" Cannot open file....'"<<fname<<"'"<<std::endl;
		}
		return false;
	}
	return true;
}




RiffFile &RiffFile::operator=(RiffFile &rhs)
{
	data=rhs.data;
	totalSize=rhs.totalSize;
	return *this;
}


RiffFile::~RiffFile()
{	if (fp.is_open())	fp.close();	}


bool RiffFile::empty()const {	return data.empty();	}


bool RiffFile::clear()
{

// clear the data stack
	while (!data.empty())	data.pop();

// rewind to the start of the file
	if (fp.is_open())	fp.seekg(std::ios::beg);

// look for a valid RIFF header
	RiffData top(*this);

//we should not be at the end of the file at this point
// and the first name should be 'RIFF' in the file
	if (fp.eof() || std::strcmp(top.name, "RIFF"))	return false;

// put it on the stack,
// and reset the 'total' size
// as the get pointer.
	totalSize = top.size;
	data.push(top);

//everything worked fine
	return true;
}

bool RiffFile::push(const char* cType)
{

// cannot do anything if we've not start properly
	if (data.empty())	return false;

// first, go to the start of the current chunk,
// if we're looking for a named
// chunk.
	if (cType)
	{
		fp.seekg(data.top().start);
		if(fp.fail()) return false;
	}

// read data until one matches or we exhaust this chunk
	while (!fp.eof() && fp.tellg() < data.top().end)
	{
		RiffData chunk(*this);

		if (!fp.eof())
		{
		// see if the subchunk type matches

			if (!cType || strcmp(chunk.name, cType) == 0){
			// push the chunk, and succeed
				data.push(chunk);
				return true;
			}else{

			//  go to the next one.
				fp.seekg(chunk.end);
				if(fp.eof())	return false; //gone too far
			}
		}
	}

// couldn't find it;
// put ourselves back from whence we start
// and return error.
	fp.seekg(data.top().start);
	return false;
}

bool RiffFile::pop()
{

// if only the top level chunk (or not even that),
// then we cannot pop upwards (as there is nothing to go to).
	if (data.size() < 2)	return false;

// move the file ptr to the place end this current
// one we are about to pop
	fp.seekg(data.top().end);

// Pop up the stack.
	data.pop();
	return true;
}

long RiffFile::currentSize() const
{
	return (!empty())?(data.top().size):0;
}

const char* RiffFile::currentName() const
{
	return (!empty())?(data.top().name):0;
}

const char* RiffFile::currentSubType() const
{
	return (!empty() && data.top().subType[0])?(data.top().subType):0;
}

bool RiffFile::getNextItem(std::string& type, std::string& value)
{

// if the current chunk is LIST/INFO,
// then try to read another sub data element.
	if (std::strcmp(currentName(), "LIST") == 0	&&
	    std::strcmp(currentSubType(), "INFO") == 0)
	{
		if (push()) {

		//got one we know
			if (readItem(type, value))		return true;
		//unknown type..try to get the next one
			else	return getNextItem(type, value);
		} else {
			// got to the end of the LIST/INFO chunk.
			// move back out and continue trying to find the one we want.
			pop();
			return getNextItem(type, value);
		}
// not in a LIST/INFO chunk,
// so look for the next DISP or LIST/INFO.
	} else {
		push();

	// DISP chunk: read and pop back out.
		if (std::strcmp(currentName(), "DISP") == 0) {
			return readItem(type, value);

	// LIST/INFO chunk: read first element
		} else if (std::strcmp(currentName(), "LIST") == 0 &&
			       std::strcmp(currentSubType(), "INFO") == 0)
		{
			return getNextItem(type, value);

	// Some other chunk  Pop and move on.
		}else{
			return pop()?getNextItem(type, value):false;
		}
	}
}

// reads extra data from the current chunk, and pops out of it.
bool RiffFile::readItem(std::string& type, std::string& value)
{
// see if it is in our list of types above
	bool found = false;

	for (int i = 0; i < BL__numExtraTypes__; i++) {

	//got one
		if (std::strcmp(currentName(), BL__RiffDataTypes[i].typeName.c_str()) == 0)
		{
			type = BL__RiffDataTypes[i].realName;
			found = true;
		}
	}

// DISP data skip four bytes before the display name starts.
	if (std::strcmp(currentName(), "DISP") == 0) {
		char tmC[5];
		fp.read(tmC, 4*sizeof(char));
	}

// read the value, if we got one we know
	if (found)
	{
		int c;
		fp.read((char *)&c, 1);
		value = "";
		while (c != '\0' && c != EOF)
		{
			fp.read((char *)&c, 1);
			value += char(c);
		}
	}

// whether we recognize it or not, pop out.
	pop();
	return found;
}

std::ostream & RiffFile::print(std::ostream &oo, bool expand, int indent)
{
	for (int i = 0; i < indent; i++)oo << ' ';

	oo << "Chunk type: " << currentName();

	if (currentSubType())
		oo << " (" << currentSubType() << ")";

	oo << ", " << currentSize() << " bytes" << std::endl;

// show all SUB DATA BITS
	if (std::strcmp(currentName(), "RIFF") == 0 ||
	   (expand &&
	   std::strcmp(currentName(), "LIST") == 0))
	{
		while (push())
		{
			print(oo, expand, indent+2);
			pop();
		}
	}
	return oo;
}

std::ostream &operator<<(std::ostream &oo, RiffFile &out)
{	return out.print(oo, true);	}


/**** WAVE FILE *****/

/********** Setters of the 'options' in the
*********** wave file...out litte 'Riff' tag
***********/
unsigned short wavestream::getFormatType() const
{ return formatType; };

void wavestream::setFormatType(unsigned short type)
{ formatType = type; };

bool wavestream::isCompressed() const
{ return formatType != 1; };

unsigned short wavestream::getNumChannels() const
{ return numChannels; };

void wavestream::setNumChannels(unsigned short num)
{ numChannels = num;  };

unsigned long wavestream::getSampleRate() const
{ return sampleRate; };

void wavestream::setSampleRate(unsigned long rate)
{ sampleRate = rate;  };

unsigned long wavestream::getBytesPerSecond() const
{ return bytesPerSecond; };

void wavestream::setBytesPerSecond(unsigned long bytes)
{ bytesPerSecond = bytes;  };

unsigned short wavestream::getBytesPerSample() const
{ return bytesPerSample; };

void wavestream::setBytesPerSample(unsigned short bytes)
{ bytesPerSample = bytes;  };

short int wavestream::getBitsPerChannel() const
{ return bitsPerChannel; };

void wavestream::setBitsPerChannel(unsigned short bits)
{ bitsPerChannel = bits;  };

unsigned long wavestream::getNumSamples() const
{
	return (getBytesPerSample())?
		getDataLength() / getBytesPerSample(): 0;
};

void wavestream::setNumSamples(unsigned long num)
{ setDataLength(num * getBytesPerSample()); };

float wavestream::getNumSeconds() const
{
	return getBytesPerSecond()?
	float(getDataLength()) / getBytesPerSecond(): 0;
};

long int wavestream::getDataLength() const
{ return dataLength; };

void wavestream::setDataLength(long int numBytes)
{ dataLength = numBytes; };

/******* Master Setter of the wave parameters *****/
void wavestream::setupFormat(int sampleRate, short bitsPerChannel, short channels)
{
	setFormatType(1);
	setNumChannels(channels);
	setSampleRate(sampleRate);
	setBytesPerSample((unsigned short)((bitsPerChannel >> 3) * channels));
	setBytesPerSecond(sampleRate * getBytesPerSample());
	setBitsPerChannel(bitsPerChannel);
	setNumSamples(0);
}

/**** Other nessesary functions ****/

// constants for the WAVE format
// length of fmt contents
const int BL__fmtChunkLength = 16;

// from "WAVE" to sample data
const int BL__waveHeaderLength = 4 + 8 + BL__fmtChunkLength + 8;


wavestream::wavestream():
	formatType(0),
	numChannels(0),
	sampleRate(0),
	bytesPerSecond(0),
	bytesPerSample(0),
	bitsPerChannel(0),
	dataLength(0),
	error(0)
{
	numSteps=1;
	bufferLen_=0;
	delBuffer=false;
}

wavestream::wavestream(const char *fname, std::ios::openmode iom):
	formatType(0),
	numChannels(0),
	sampleRate(0),
	bytesPerSecond(0),
	bytesPerSample(0),
	bitsPerChannel(0),
	dataLength(0),
	error(0)
{
	open(fname, iom);
}

//if the iomode is to write, we need to write the
// header then close the file...otherwise (if read)
// we only need to close the file.
bool wavestream::close()
{
	if(file().is_open())
	{
		if(iomode & std::ios::out){
			//go to the begining of the file
			file().seekg(std::ios::beg);
			//write the header
			//calculate the total number of smaples (and bytes) in our buffer
			setDataLength(bufferLen_);
			if(!writeHeader()){
				std::cerr<<"wavestream::close(..)"<<std::endl;
				std::cerr<<"could not write Header.. "<<std::endl;
				if(delBuffer) delete [] buffer_;
				riffFile.close();
				return false;
			}

		file().write(buffer_, bufferLen_*sizeof(char));
		//delete our buffer and reset the numSteps counter
			if(delBuffer) delete [] buffer_;
			numSteps=1;
			bufferLen_=0;
			delBuffer=false;

		}

		return riffFile.close();
	}
	return true;
}

wavestream::wavestream(wavestream& other)
{
	formatType = other.formatType;
	numChannels = other.numChannels;
	sampleRate = other.sampleRate;
	bytesPerSecond = other.bytesPerSecond;
	bytesPerSample = other.bytesPerSample;
	bitsPerChannel = other.bitsPerChannel;
	riffFile=other.riffFile;
	buffer_=other.buffer_;
	bufferLen_=other.bufferLen_;
}

wavestream &wavestream::operator=( wavestream& other)
{
	formatType = other.formatType;
	numChannels = other.numChannels;
	sampleRate = other.sampleRate;
	bytesPerSecond = other.bytesPerSecond;
	bytesPerSample = other.bytesPerSample;
	bitsPerChannel = other.bitsPerChannel;
	riffFile=other.riffFile;
	bufferLen_=other.bufferLen_;
	if(delBuffer) delete [] buffer_;
	numSteps=other.numSteps;

	if(bufferLen_ >0){
		buffer_=new  char[bufferLen_];
		std::memcpy(buffer_, other.buffer_, bufferLen_*sizeof(char));
	}

	return *this;
}

wavestream::~wavestream()
{
	close();
}

bool wavestream::operator==(const wavestream& other)
{
	return formatType == other.formatType
		&& numChannels == other.numChannels
		&& sampleRate == other.sampleRate
		&& bytesPerSecond == other.bytesPerSecond
		&& bytesPerSample == other.bytesPerSample
		&& bitsPerChannel == other.bitsPerChannel;
}

void wavestream::copyFormat(const wavestream& other)
{
	formatType = other.formatType;
	numChannels = other.numChannels;
	sampleRate = other.sampleRate;
	bytesPerSecond = other.bytesPerSecond;
	bytesPerSample = other.bytesPerSample;
	bitsPerChannel = other.bitsPerChannel;
}

bool wavestream::reset()
{
	if(file().is_open()){
		if(iomode & std::ios::out){
			file().seekg(std::ios::beg);
			return true;
		}else if(std::ios::in){
			if(!riffFile.clear() ||
			   !riffFile.push("data"))
			{
				std::cerr<<"wavestream::reset(..)"<<std::endl;
				std::cerr<<"could not find the data on 'reset'"<<std::endl;
				return false;
			}else{
				return true;
			}
		}
	}
	return true;
}

bool wavestream::open(const char* fname, std::ios::openmode iomo)
{
	iomode=iomo;
	if(file().is_open()) file().close();

	if(iomode & std::ios::in){
		try {
		// open the RIFF file
			delBuffer=false;
			riffFile.open(fname);
			if(file().fail()) throw error = "Couldn't open file";

		// read the header information
			if (std::strcmp(riffFile.currentName(), "RIFF") ||
				std::strcmp(riffFile.currentSubType(), "WAVE") ||
				!riffFile.push("fmt "))
			{
				throw error = "Couldn't find RIFF, WAVE, or fmt";
			}

			unsigned long FmtSize = riffFile.currentSize();
			char* Chunk = new char[FmtSize];
			try {
				file().read(Chunk, FmtSize*sizeof(char));
				if(file().fail())
					throw error = "Error reading format chunk";

				riffFile.pop();

				// set the format attribute members
				formatType = *((short*) Chunk);
				numChannels = *((short*) (Chunk + 2));
				sampleRate = *((long*) (Chunk + 4));
				bytesPerSecond = *((long*) (Chunk + 8));
				bytesPerSample = *((short*) (Chunk + 12));
				bitsPerChannel = *((short*) (Chunk + 14));

				// position at the data chunk
				if (!riffFile.push("data"))
					throw error = "Couldn't find data chunk";

				// get the size of the data chunk
				dataLength = riffFile.currentSize();

				delete[] Chunk;
			} catch (...) {
				std::cerr<<"wavestream::open(..)"<<std::endl;
				std::cerr<<error<<std::endl;

				delete[] Chunk;
				throw error;
			}
		} catch (...) {
			std::cerr<<"wavestream::open(..)"<<std::endl;
			std::cerr<<error<<std::endl;
			close();
			return false;
		}
		return true;
	}else if(iomode & std::ios::out){
		file().open(fname, iomode);
		if(file().fail()){
			std::cerr<<"wavestream::open(..)"<<std::endl;
			std::cerr<<"could not open file.. "<<std::endl;
			return false;
		}
		numSteps=1;
		BufferStep=10*1024;
		buffer_=new char[BufferStep];
		bufferLen_=0;
		delBuffer=true;

	//	if(!writeHeader()){
	//		std::cerr<<"wavestream::open(..)"<<std::endl;
	//		std::cerr<<"could not write Header.. "<<std::endl;
	//		return false;
	//	}
		return true;
	}
	return false;
}

/******* Aux Functions for finding chunks in the Riff File ***/
bool wavestream::getFirstItem(std::string& type, std::string& value)
{
	return (riffFile.clear() && riffFile.getNextItem(type, value));
}

bool wavestream::getNextItem(std::string& type, std::string& value)
{
	return (riffFile.getNextItem(type, value));
}


// READERS AND WRITERS

//writes the header
bool wavestream::writeHeader()
{
// move to the start of the file
	file().seekg(std::ios::beg);

//no file/failure
	if(file().fail())	return false;

	// write the file header
	unsigned long wholeLength = BL__waveHeaderLength + dataLength;
	unsigned long chunkLength = BL__fmtChunkLength;

	file()<<"RIFF";
	file().write((char *)&wholeLength, sizeof(wholeLength));
std::cout<<"Total Wav File size: "<<wholeLength<<std::endl;
	file()<<"WAVE";
	file()<<"fmt ";
	file().write((char *)&chunkLength, sizeof(chunkLength));
	file().write((char *)&formatType, sizeof(formatType));
	file().write((char *)&numChannels, sizeof(numChannels));
	file().write((char *)&sampleRate, sizeof(sampleRate));
	file().write((char *)&bytesPerSecond, sizeof(bytesPerSecond));
	file().write((char *)&bytesPerSample, sizeof(bytesPerSample));
	file().write((char *)&bitsPerChannel, sizeof(bitsPerChannel));
	file()<<"data";
	file().write((char *)&dataLength, sizeof(dataLength));

	wroteHeader = true;

	return true;
}

//RAW  data reader
bool wavestream::readRaw(char* buffer, unsigned long int numBytes)
{
	file().read(buffer, numBytes);

	if (file().fail()) {
		std::cerr<<"wavestream::readRaw(..)"<<std::endl;
		std::cerr<<"could not read the data bits... "<<std::endl;
		return false;
	}

	return true;
}

//read single float
bool wavestream::read(float& sample)
{
	double fSample;
	bool retval = read(fSample);
	sample = fSample;
	return retval;
}


//read single double
bool wavestream::read(double& sample)
{
	bool retval = false;

//need to do some fancy bit shifting
// if we want a double, as the data
// is stored as 'chars' so we need for 8 bits per chanel
// and a 'short' for 16 bits per chanel
	if (getBitsPerChannel() == 8) {
		unsigned char cSample;
		retval = read(cSample);
		sample = double(cSample) / ((1 << (8 - 1)) - 1) - 1;
	} else if (getBitsPerChannel() == 16) {
		short sSample;
		retval = read(sSample);
		sample = double(sSample) / ((1 << (16 - 1)) - 1);
	} else{
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for doubles should be 8 or 16"<<std::endl;
	}

	return retval;
}

//float array
bool wavestream::read(float *sample, unsigned long ct)
{
	bool tm=true;
	for(unsigned long i=0;i<ct;++i){
	  tm &=read(sample[i]);
  	}
  	return tm;
}

//double array
bool wavestream::read(double *sample, unsigned long ct)
{
	bool tm=true;
	for(unsigned long i=0;i<ct;++i){
	  tm &= read(sample[i]);
  	}
  	return tm;

}

//read single char
bool wavestream::read(char& sample)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return readRaw(&sample,1);
};


//read single char
bool wavestream::read(unsigned char& sample)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return readRaw((char*) &sample,1);
};

//read short single
bool wavestream::read(short& sample)
{
	if (getBitsPerChannel() != 16) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for short should be 16"<<std::endl;
		return false;
	}

	return readRaw((char*) &sample, 2);
};

//read 'int' single
bool wavestream::read(int &samples)
{
	short tmp;
	if(read(tmp)){
	 	samples=tmp;
	 	return true;
	}
	return false;
}

//read 'complex float' single
bool wavestream::read(Complex<float> &samples)
{
	//short tmp;
	if(read(samples.Re()) && getNumChannels()==2){
	 	if(read(samples.Im())){
	 		return true;
		}
	}
	return false;
}

//read 'complex double' single
bool wavestream::read(Complex<double> &samples)
{
	//short tmp;
	if(read(samples.Re())&& getNumChannels()==2){
	 	if(read(samples.Im())){
	 		return true;
		}
	}
	return false;
}


//read char array
bool wavestream::read(unsigned char* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return readRaw((char*) samples, getNumChannels() * count);
}

//read char array
bool wavestream::read(char* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return readRaw((char*) samples, getNumChannels() * count);
}

//read 'short' samples
bool wavestream::read(short* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 16) {
		std::cerr<<"wavestream::read(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for short should be 16"<<std::endl;
		return false;
	}

	return readRaw((char*) samples, 2 * getNumChannels() * count);
}

//read 'int' samples
bool wavestream::read(int* samples, unsigned long int count)
{
	bool tm=true;
	short tmp;
	for(unsigned int i=0;i<count;++i){
		tm &= read(tmp);
		if(tm) samples[i]=tm;
	}
	return tm;
}

/** VECTOR READERS **/
bool wavestream::read(Vector<unsigned char> &sample)
{
	sample.resize(getNumSamples());
	return read(sample.data(), (unsigned long)sample.size());
}

bool wavestream::read(Vector<char> &sample)
{
	sample.resize(getNumSamples());
	return read((unsigned char *)sample.data(), (unsigned long)sample.size());
}

bool wavestream::read(Vector<short> &sample)
{
	sample.resize(getNumSamples());
	return read(sample.data(), (unsigned long)sample.size());
}

bool wavestream::read(Vector<int> &sample)
{
	sample.resize(getNumSamples());
	return read(sample.data(), (unsigned long)sample.size());
}

bool wavestream::read(Vector<float> &sample)
{
	sample.resize(getNumSamples());
	bool tm=true;
	for(int i=0;i<sample.size();++i){
		tm&= read(sample(i));
	}
	return tm;
}

bool wavestream::read(Vector<double> &sample)
{
	sample.resize(getNumSamples());
	bool tm=true;
	for(int i=0;i<sample.size();++i){
		tm&= read(sample(i));
	}
	return tm;
}

/** MAtrix Readers **/
//MATRIX readers...these will either give a 2xN or 1xN depending
// on the number of channels!!
bool wavestream::read(_matrix<unsigned char, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

bool wavestream::read(_matrix<char, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

bool wavestream::read(_matrix<short, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

bool wavestream::read(_matrix<int, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

bool wavestream::read(_matrix<float, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

bool wavestream::read(_matrix<double, FullMatrix> &sample)
{
	sample.resize(1, getNumSamples());
	bool tm=read(sample.data(), (unsigned long)sample.cols());
	if(getNumChannels()==2 && tm){
		sample.reshape(2,getNumSamples()/2);
	}
	return tm;
}

//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::read(_matrix<Complex<float>, FullMatrix> &sample)
{
	_matrix<float, FullMatrix> tmpD;
	bool tm=read(tmpD);
	if(getNumChannels()==2 && tm){
		sample.resize(1, getNumSamples()/2);
		sample.putCol(1, tmpD.row(0)+Complex<float>(0,1)*tmpD.row(1));
	}else if(tm){
		sample=tmpD;
	}
	return tm;
}

//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::read(_matrix<Complex<double>, FullMatrix> &sample)
{
	_matrix<double, FullMatrix> tmpD;
	bool tm=read(tmpD);
	if(getNumChannels()==2 && tm){
		sample.resize(1, getNumSamples()/2);
		sample.putCol(1, tmpD.row(0)+Complex<double>(0,1)*tmpD.row(1));
	}else if(tm){
		sample=tmpD;
	}
	return tm;
}

//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::read(Vector<Complex<float> > &sample)
{
	_matrix<float, FullMatrix> tmpD;
	bool tm=read(tmpD);
	if(getNumChannels()==2 && tm){
		sample.resize(getNumSamples()/2);
		sample=tmpD.row(0)+Complex<float>(0,1)*tmpD.row(1);
	}else if(tm){
		sample=tmpD.row(0);
	}
	return tm;
}

//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::read(Vector<Complex<double> > &sample)
{
	_matrix<double, FullMatrix> tmpD;
	bool tm=read(tmpD);
	if(getNumChannels()==2 && tm){
		sample.resize(getNumSamples()/2);
		sample=tmpD.row(0)+Complex<double>(0,1)*tmpD.row(1);
	}else if(tm){
		sample=tmpD.row(0);
	}
	return tm;
}

/***********************
****************************
***************WRITES **/


//RAW
bool wavestream::writeRaw(char* INbuffer, unsigned long int numBytes)
{
	while(bufferLen_+numBytes>=numSteps*BufferStep)
	{
		char *tmPBuf; tmPBuf=new char[bufferLen_]; //tmp buffer

	//copy exsisting data into temp
		std::memcpy(tmPBuf, buffer_, bufferLen_*sizeof(char));

	//delete our old buffer
		delete [] buffer_;

	//create a new one of length ((numSteps+1)*BufferStep)
		buffer_ = new char[(numSteps+1)*BufferStep];

	//copy the data back to our master buffer
		std::memcpy(buffer_,tmPBuf,bufferLen_*sizeof(char));

	//advnace the buf steps
		numSteps++;

	//delete our tmp buffer
		delete [] tmPBuf;
	}

//now place the new data into out buffer
	std::memcpy(&buffer_[bufferLen_], INbuffer, numBytes*sizeof(char));

//advance the bufferlen
	bufferLen_+=numBytes;
	return true;
}

//	file().write(buffer, numBytes);
//	if(file().fail())
//	{
//		std::cerr<<"wavestream::write(..)"<<std::endl;
//		std::cerr<<"Cannot write to file... "<<std::endl;
//		return false;
//	}
//	setDataLength(getDataLength() + numBytes);

//	return true;
//}

//unsigned char array
bool wavestream::write(unsigned char* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return writeRaw((char*) samples, getNumChannels() * count);
}

//unsigned char array
bool wavestream::write( char* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return writeRaw((char*) samples, getNumChannels() * count);
}

//write 'int' array
bool wavestream::write(int* samples, unsigned long int count)
{
	bool tm=true;
	short tmp;
	for(unsigned long int i=0;i<count;++i){
		tmp=samples[i];
		tm &= write(tmp);
	}
	return tm;
}

//write 'int' single
bool wavestream::write(int samples)
{
	short tmp=samples;
	return write(tmp);
}

//unsigned char single
bool wavestream::write(unsigned char sample)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return writeRaw((char*) &sample);
};

bool wavestream::write( char sample)
{
	if (getBitsPerChannel() != 8) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for char should be 8"<<std::endl;
		return false;
	}

	return writeRaw( &sample);
};

//write short single
bool wavestream::write(short sample)
{
	if (getBitsPerChannel() != 16) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for short should be 16"<<std::endl;
		return false;
	}

	return writeRaw((char*) &sample, 2);
};
//short array
bool wavestream::write(short* samples, unsigned long int count)
{
	if (getBitsPerChannel() != 16) {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for short should be 16"<<std::endl;
		return false;
	}

	return writeRaw((char*) samples, 2 * getNumChannels() * count);
}

//float single
bool wavestream::write(float sample)
{
	return write(double(sample));
}

//coplex float single
bool wavestream::write(Complex<float> sample)
{
	write(Re(sample));
	return write(Im(sample));
}

//coplex double single
bool wavestream::write(Complex<double> sample)
{
	write(Re(sample));
	return write(Im(sample));
}

//double single
bool wavestream::write(double sample)
{
	if (getBitsPerChannel() == 8)
		return write((unsigned char)((sample + 1) * ((1 << (8 - 1)) - 1)));
	else if (getBitsPerChannel() == 16)
		return write(short(sample * ((1 << (16 - 1)) - 1)));
	else {
		std::cerr<<"wavestream::write(..)"<<std::endl;
		std::cerr<<"Sample size mismatch (bitsPerChannel): ";
		std::cerr<<" for doubles should be 8 or 16"<<std::endl;
		return false;
	}
}

//float single
bool wavestream::write(float *sample, unsigned long ct)
{
	bool tm=true;
	for(unsigned long i=0;i<ct;++i){
	  tm &= write(double(sample[i]));
  	}
  	return tm;
}

//double single
bool wavestream::write(double *sample, unsigned long ct)
{
	bool tm=true;
	for(unsigned long i=0;i<ct;++i){
	  tm &= write(sample[i]);
  	}
  	return tm;

}

/*** VECTOR writers ****/
bool wavestream::write(Vector<unsigned char> sample)
{	return write(sample.data(), (unsigned int)sample.size()); }

bool wavestream::write(Vector<char> sample)
{	return write((unsigned char *)sample.data(), (unsigned long)sample.size()); }

bool wavestream::write(Vector<short> sample)
{	return write(sample.data(), (unsigned long)sample.size()); }

bool wavestream::write(Vector<int> sample)
{	return write(sample.data(), (unsigned long)sample.size()); }

bool wavestream::write(Vector<float> sample)
{	return write(sample.data(), (unsigned long)sample.size()); }

bool wavestream::write(Vector<double> sample)
{	return write(sample.data(), (unsigned long)sample.size()); }

//double complex
bool wavestream::write(Vector<Complex<double> > sample)
{
	bool tm=true;
	setNumChannels(2);
	for(int i=0;i<sample.size();++i){
	  tm &= write(Re(sample[i]));
	  tm &= write(Im(sample[i]));
  	}
  	return tm;

}

//double single
bool wavestream::write(Vector<Complex<float> > sample)
{
	bool tm=true;
	setNumChannels(2);
	for(int i=0;i<sample.size();++i){
	  tm &= write(Re(sample[i]));
	  tm &= write(Im(sample[i]));
  	}
  	return tm;

}

/** MAtrix Readers **/
//MATRIX readers...these will either give a 2xN or 1xN depending
// on the number of channels!!
bool wavestream::write(_matrix<unsigned char, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}

bool wavestream::write(_matrix<char, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}

bool wavestream::write(_matrix<short, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}

bool wavestream::write(_matrix<int, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}

bool wavestream::write(_matrix<float, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}


bool wavestream::write(_matrix<double, FullMatrix> sample)
{
	if(sample.rows()==2){
		setNumChannels(2);
		return write(sample.data(), (unsigned long)sample.cols());
	}else if(sample.rows()==1){
		setNumChannels(1);
		return write(sample.data(), (unsigned long)sample.cols());
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write"<<std::endl;
		std::cerr<<" output matrix must be 1xN or 2xN..."<<std::endl;
		return false;
	}
	return false;
}

//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::write(_matrix<Complex<float>, FullMatrix> sample)
{
	bool tm=true;
	if(sample.rows()==1){
		setNumChannels(2);
		for(int i=0;i<sample.cols();++i){
			tm&= write(Re(sample(0,i)));
			tm&= write(Im(sample(0,i)));
		}
		return tm;
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write(complex matrix)"<<std::endl;
		std::cerr<<" output matrix must be 1xN..."<<std::endl;
		return false;
	}
	return tm;
}


//the complex types will not get 'reshaped' instead
// the Right channel is the real, the Left Channel is in the
// imag, (if channels==2, otherwise only the 'real' is filled)
bool wavestream::write(_matrix<Complex<double>, FullMatrix> sample)
{
	bool tm=true;
	if(sample.rows()==1){
		setNumChannels(2);
		for(int i=0;i<sample.cols();++i){
			tm&= write(Re(sample(0,i)));
			tm&= write(Im(sample(0,i)));
		}
		return tm;
	}else{
		std::cerr<<std::endl<<"Error: wavestream::write(complex matrix)"<<std::endl;
		std::cerr<<" output matrix must be 1xN..."<<std::endl;
		return false;
	}
	return tm;
}

std::ostream & wavestream::print(std::ostream &oo)
{
	oo
		<< "Format:           " << getFormatType()
		<< (isCompressed()? " (compressed)" : " (PCM)") << std::endl
		<< "Channels:         " << getNumChannels() << std::endl
		<< "Sample rate:      " << getSampleRate() << std::endl
		<< "Bytes per second: " << getBytesPerSecond() << std::endl
		<< "Bytes per sample: " << getBytesPerSample() << std::endl
		<< "Bits per channel: " << getBitsPerChannel() << std::endl
		<< "Bytes:            " << getDataLength() << std::endl
		<< "Samples:          " << getNumSamples() << std::endl
		<< "Seconds:          " << getNumSeconds() << std::endl
		<< "File pointer:     " << (file()? "good" : "null") << std::endl;

	return oo;
}

std::ostream & operator<<(std::ostream & oo, wavestream &out)
{	return out.print(oo);	}

//Writes

wavestream &operator<<(wavestream &oo, char out)
{	oo.write(out);	return oo; }
wavestream &operator<<(wavestream &oo, Vector< char> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix< char, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo,unsigned char out)
{	oo.write(out);	return oo; }
wavestream &operator<<(wavestream &oo, Vector<unsigned char> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<unsigned char, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, short out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<short> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<short, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, int out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<int> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<int, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, float out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<float> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<float, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, double out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<double> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<double, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, Complex<float> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<Complex<float> > out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<Complex<float>, FullMatrix> out)
{	oo.write(out);	return oo;}

wavestream &operator<<(wavestream &oo, Complex<double> out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, Vector<Complex<double> > out)
{	oo.write(out);	return oo;}
wavestream &operator<<(wavestream &oo, _matrix<Complex<double>, FullMatrix> out)
{	oo.write(out);	return oo;}


//READS
wavestream &operator>>(wavestream &oo,unsigned char &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<unsigned char> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<unsigned char, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, char &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<char> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<char, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, short &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<short> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<short, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, int &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<int> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<int, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, float &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<float> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<float, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, double &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<double> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<double, FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, Complex<float> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<Complex<float> > &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<Complex<float> , FullMatrix> &out)
{	oo.read(out);	return oo;}

wavestream &operator>>(wavestream &oo, Complex<double> &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, Vector<Complex<double> > &out)
{	oo.read(out);	return oo;}
wavestream &operator>>(wavestream &oo, _matrix<Complex<double> , FullMatrix> &out)
{	oo.read(out);	return oo;}


END_BL_NAMESPACE

#endif
