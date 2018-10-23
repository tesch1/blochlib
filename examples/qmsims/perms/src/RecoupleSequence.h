


/*
this program calculates the tensorial components of various permutations
on the standard recoupling sequence Post-C7
*/

#ifndef _RecouplePermutions__h_
#define _RecouplePermutions__h_ 1

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

class PulseTrainIterator;

//Generic Pulse train container
class PulseTrain{
	public:
		Vector<double> amplitudes;
		Vector<double> phases;
		Vector<double> times;

		typedef PulseTrainIterator iterator;

		PulseTrain(){}
		PulseTrain(const Vector<double> &amp, const Vector<double> &phases,const Vector<double> &times);

		void operator=(const PulseTrain &rhs);

		inline double beginTime(int i){	return times(i);	}
		inline double endTime(int i){	return times(i+1);	}
		inline double phase(int i){		return phases(i);	}
		inline double amplitude(int i){	return amplitudes(i);	}
		inline double finalTime(){	return times(times.size()-1);	}

		inline int size() const {	return amplitudes.size();	}

		matrix Pulse(SpinSys &in, std::string on, int which);

//		matrix propagator(SolidSys &in, Rotations &rots, std::string on);

};

class PulseTrainIterator
{
	private:
		PulseTrain *data;
		int curpos;
		//

	public:
		PulseTrainIterator():
			data(NULL)
		{}

		PulseTrainIterator(PulseTrain &in):
			data(&in), curpos(0)
		{}

		~PulseTrainIterator(){	data=NULL;	}

		void operator=(PulseTrainIterator &rhs);

		inline double beginTime(){	return data->times(curpos);	}
		inline double endTime(){	return data->times(curpos+1);	}
		inline double phase(){		return data->phases(curpos);	}
		inline double amplitude(){	return data->amplitudes(curpos);	}
		inline double angle(){	return PI2*amplitude()*(endTime()-beginTime());	}

		inline void operator++(){	++curpos;	}
		inline void operator++(int){	++curpos;	}

		inline void reset(){	curpos=0;	}
		operator bool(){	return curpos<data->size();	}

		matrix Pulse(SpinSys &in, std::string on);

};


//this class holds the basic info for a recoupling sequence
// and also generates them.....
// the pulse amplitudes are always the same for these sequences
// only the angle changes (i.e. the application time)
class Recouple :
	public PulseTrain
{
	friend class RecoupleTrain;

	private:
		int rtype_;
		int baseSym_; // phases propagate as 2Pi*X/N..this is N
		int fractSym_;	//phases propagate as 2Pi*X/N..this is X

		double amp_; //pulse amplitude IN HZ
		double sphase_; //starting phase
		double stime_; //starting time
	public:

		enum{
			R=	0x000001, // angle=(360, 360) phases=(0, -0; phi1, -phi1;...)
			Rp=	0x000002, // angle=(90, 270, -90, -270) phase=(0, 180, -0, -180; phi1, phi1+180, -phi1, -phi1-180;...)
			C=	0x000004, // angle=(360, 360) phases=(0, 0+Pi; phi1, phi1+Pi;...)
			Cp=	0x000008,  // angle=(90, 360, 270) phases(0, 0+Pi,0; phi1, phi1+Pi, phi1 ...)
			Cpp=0x000010  // angle=(90, 360, 270) phases(0, 0+Pi,0; phi1, phi1+Pi, phi1 ...)
		}RecoupleTypes;

		typedef PulseTrainIterator iterator;

		static int StringToRtype(std::string rtype);
		static std::string RtypeToString(int rtype);

		Recouple(){}
		Recouple(int rtype, int baseSym, int fractSym,  double amplitude, double phaseoff=0.0, double start=0.0);
		Recouple(std::string rtype, int baseSym, int fractSym, double amplitude, double phaseoff=0.0, double start=0.0);
		Recouple(Parameters &pset, std::string sec="");

		void operator=(const Recouple &rhs);

		void setRecoupleType(int rtype);
		void setBaseSym(int baseSym);
		void setFractionalSym(int fracsym);

		void setSequence(int rtype, int basesym, int fracsym);

		void setAmplitude(double amp);
		void setStartPhase(double sphase);
		void setStartTime(double stime);

		void setPulseParams(double amp, double sphase, double stime);

		void read(Parameters &pset, std::string sec="");
		void init(double amp=50000.0, double startp=0.0, double startt=0.0);

		inline int size() const {	return amplitudes.size();	}

		inline double endTime() const {	return times[size()];	}
		inline double beginTime() const {	return times[0];	}

		friend ostream &operator<<(ostream &oo,const Recouple &out);

};


/*******************************************************************/
/**** Post -C7 train  componets   **/
/**** 1) the data struct **/
/**** 2) the iterator **/
/*******************************************************************/

class RecoupleTrainIter;
/*******************************************************************/
/**** Post -C7 train Itterator  **/
/*******************************************************************/

class RecoupleTrain{
	public:

		typedef RecoupleTrainIter iterator;

		Vector<Recouple> train;

		RecoupleTrain(){}

		RecoupleTrain(int len);
		RecoupleTrain(const Vector<Recouple> &inp );
		RecoupleTrain(Parameters &pset, std::string sec="");

		void setBeginTime(double inbtime);

		inline int size() const {	return train.size();	}

		inline Recouple &element(int i){	return train[i];	}
		inline Recouple element(int i) const {return train[i];	}

		inline double beginTime(){	return train[0].beginTime();	}
		inline double endTime(){	return train[size()-1].endTime();	}

		inline Recouple &operator[](int i){	return train[i];	}
		inline Recouple operator[](int i)const{	return train[i];	}
		inline Recouple &operator()(int i){	return train[i];	}
		inline Recouple operator()(int i)const{	return train[i];	}

//		matrix propagator(SolidSys &in, Rotations &rots, std::string on);

};

ostream &operator<<(ostream &oo,const RecoupleTrain &out);

/*******************************************************************/
/**** Post -C7 train Itterator  **/
/*******************************************************************/

class RecoupleTrainIter{
	private:
		RecoupleTrain *data;
		Recouple::iterator curRecouple;

		int curpos;


	public:
		RecoupleTrainIter(RecoupleTrain &in);
		void operator=(RecoupleTrainIter &in);

		void operator++();
		void operator++(int);

		void reset();

		inline operator bool(){	return curpos<data->size();	}

		inline double beginTime(){	return curRecouple.beginTime();	}
		inline double endTime(){	return curRecouple.endTime();	}
		inline double phase(){		return curRecouple.phase();	}
		inline double amplitude(){	return curRecouple.amplitude();	}
		inline double angle(){	return curRecouple.angle();	}

		matrix Pulse(SpinSys &in, std::string on);
};


#endif
