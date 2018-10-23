


#ifndef _Pulse_h_
#define _Pulse_h_

#include <string.h>
#include <math.h>
#include "container/Vector/Vector.h"
#include "container/grids/coords.h"
#include <iostream>

#ifdef __WIN32__
 #include "blochconfigCW.h"
#else
 #include "blochconfig.h"
#endif

#ifdef HAVE_CLIMITS
 #include <climits>
#else
 #include <limits.h>
#endif

BEGIN_BL_NAMESPACE



/* Pulses:: two classes here
	SPulse --> contains the data for a pulse on ONe spin type
		name-->"1H" or "13C", etc
		amplitude--> amplitude of the pulse in Hz!
		phase --> pulse phase...(RADIANS)!!
		offset--> pulse frequency offset in Hz from the resonant condition
		tbeing --> time the pulse begings
		tend--> time the pulse ends

	Pulse --> maintinas a list of 'SPulse' for pulseing different pulses

  this class will end up being a template parameter in the the 'Bloch' class
  so there will be two distingct Bloch classes one that takes a Pulse as a template
  parameter and one that does not (i.e. no pulse) this is for speed purposes
  as getting all the data from the pulses can be time consuming (especcialy
  if performing the lookup on over 9000 spins)

 */

 class SPulse {
	 private:
	 	std::string name_;
	 	double amp_;
	 	double phase_;
	 	double sph_; double cph_;	//the sin and cos of the phases.
	 	double off_;
	 	double tb_;
	 	double te_;

	 	coord<> tmm_;		//tmp coord to store Xpart, Ypart, Zpart in...avoids construction many times.

	 	bool applyall_;		//apply to ALL the times requested...(ignore the time flags)
	 						//default will be TRUE!!

	 public:
	 	enum{Start=INT_MIN, End=INT_MAX};

		SPulse():
	 		name_("1H"),amp_(0),phase_(0), sph_(0),cph_(0),
	 		off_(0), tb_(Start), te_(End),applyall_(true){}

	 	SPulse(const SPulse &cp):
	 		name_(cp.name_),amp_(cp.amp_), phase_(cp.phase_),sph_(cp.sph_),cph_(cp.cph_),
	 		off_(cp.off_), tb_(cp.tb_), te_(cp.te_),applyall_(cp.applyall_){}

	 	SPulse(std::string nm, double amp=0.0,double phase=0.0, double offset=0.0, double tbegin=Start, double tend=End):
	 		name_(nm), amp_(amp), phase_(phase),sph_(sin(phase)),cph_(cos(phase)),
	 		off_(offset), tb_(tbegin), te_(tend),
	 		applyall_(tbegin==Start && tend==End)
	 	{}

	 	~SPulse(){}

	 	SPulse &operator=(const SPulse &rhs);

	 	inline std::string symbol()			{	return name_;	}
	 	inline void SetSymbol(std::string nm)	{	name_=nm;		}

	 	inline double amplitude()			{	return amp_;	}
	 	inline double amplitude(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return amp_;
			return 0.;
		}
		inline void setAmplitude(double in)	{	amp_=in;		}

	 	inline double phase()			{	return phase_;	}
	 	inline double phase(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return phase_;
			return 0.;
		}
		inline void setPhase(double in)	{	phase_=in;	sph_=sin(in); cph_=cos(in);	}


	 	inline double offset()				{	return off_;	}
	 	inline double offset(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return off_;
			return 0.;
		}
	 	inline void setOffset(double in)		{	off_=in;		}

	 	inline double Xpart()			{	return amp_*cph_;	}
	 	inline double Xpart(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return amp_*cph_;
			return 0.;
		}
	 	inline double Ypart()			{	return amp_*sph_;	}
	 	inline double Ypart(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return amp_*sph_;
			return 0.;
		}

	 	inline double Zpart()			{	return off_;	}
	 	inline double Zpart(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	return off_;
			return 0.;
		}

		inline coord<> CoordPulse()
		{			return coord<>(amp_*cph_, amp_*sph_, off_);		}

		inline coord<> coordPulse()
		{	return CoordPulse();	}

	 	inline coord<> CoordPulse(double intime)
	 	{
			if(applyall_ || (intime>=tb_ && intime<=te_))	tmm_(amp_*cph_, amp_*sph_, off_);
			return tmm_;
		}
		inline coord<> coordPulse(double intime)
	 	{	return CoordPulse(intime);	}

	 	inline double beginTime()			{	return tb_;		}
	 	inline void setBeginTime(double in)	{	tb_=in;	applyall_=false; }

		inline double endTime()				{	return te_;		}
	 	inline void setEndTime(double in)		{	te_=in;	applyall_=false;	}

		inline void setTime(double tb, double te)	{	tb_=tb; te_=te; applyall_=false;	}

		double timeForAngle(double angle);

	 	void print(std::ostream &oo);
	 	void print();
};

std::ostream &operator<<(std::ostream &oo, SPulse &out);

//A NULL Pulse class...i.e. does nothing just a class to make the template
//parameters happy in 'bloch'

class NullPulse{};
typedef NullPulse NoPulse;

// Pulse class

class Pulse{
	private:
		Vector<SPulse> data_;

	public:
		Pulse(){}

		Pulse(Pulse &cp):
			data_(cp.data_)
		{}

		Pulse(SPulse &in);
		Pulse(Vector<SPulse> &in);
		Pulse(std::string nm, double amp=0.0, double phase=0.0, double offset=0.0, double tbegin=SPulse::Start, double tend=SPulse::End);

		~Pulse(){}

		inline int Exists(const std::string &nm) const
		{
			int i=0;
			for(i=0;i<data_.size();i++)
			{
				if(data_(i).symbol()==nm) return i;
			}
			return -1;
		}

		inline bool ToPulse(const std::string &nm) const
		{
			int i=0;
			for(i=0;i<data_.size();i++)
			{
				if(data_(i).symbol()==nm) return true;
			}
			return false;
		}

		inline int Exists(int i) const
		{
			if(i<data_.size() && i>=0) return i;
			return -1;
		}


	//assignments
		void operator=( SPulse &rhs);
		void operator=( Pulse &rhs);

//Symbol/ name functions..................
		std::string symbol(int nm);
		void setSymbol(int tes, std::string nm);
		inline void setSymbol( std::string nm,int tes)
		{	setSymbol(tes, nm);	}

//Amplitude Functions........................
		double amplitude(const std::string &nm);
		double amplitude(int nm);
		double amplitude(double time, int nm);
		double amplitude(double time, std::string nm);

		void setAmplitude(double amp, std::string nm);
		void setAmplitude(double amp, int nm);

//Phase functions..
		double phase(const std::string &nm);
		double phase(int nm);
		double phase(double time, int nm);
		double phase(double time, const std::string &nm);

		void setPhase(double amp,const std::string &nm);
		void setPhase(double amp, int nm);

// Offset functions................
		double offset(const std::string &nm);
		double offset(int nm);
		double offset(double time, int nm);
		double offset(double time, const std::string &nm);

		void setOffset(double off, int nm);
		void setOffset(double off, const std::string &nm);

//X parts, Yparts, and Zparts for a given time and spin..
		inline double Xpart(double tin, const std::string &nm)
		{
			int i=Exists(nm);
			if(i!=-1)
			{
				return data_(i).Xpart(tin);
			}
			return 0;
		}

		inline double Ypart(double tin, const std::string &nm)
		{
			int i=Exists(nm);
			if(i!=-1)
			{
				return data_(i).Ypart(tin);
			}
			return 0;
		}

		inline double Zpart(double tin, const std::string &nm)
		{
			int i=Exists(nm);
			if(i!=-1)
			{
				return data_(i).Zpart(tin);
			}
			return 0;
		}

		inline coord<> CoordPulse(double tin, const std::string &nm){
			int i=Exists(nm);
			if(i!=-1) return data_(i).CoordPulse(tin);
			return ZeroType<coord<> >::zero();
		}


//Time functions......................
		inline void setTime(double tb, double te, const std::string &nm);
		inline void setTime(double tb, double te, int nm);
		inline void setBeginTime(double tb, std::string nm);
		inline double beginTime(std::string nm);
		inline double beginTime(int nm);
		inline void setBeginTime(double tb, int nm);
		inline void setEndTime(double te, std::string nm);
		inline void setEndTime(double te, int nm);
		inline double endTime(std::string nm);
		inline double endTime(int nm);

//operator and other functions

		inline int size()	const {	return data_.size();	}
		double timeForAngle(double angle,const std::string &nm);

		SPulse &operator()(int i) ;
		SPulse &operator[](int i) ;

		Pulse &operator+(const SPulse &in);
		Pulse &operator+=(const SPulse &in);

		Pulse &operator+=(const Pulse &in);
		Pulse &operator+(const Pulse &in);

		void add(SPulse in);
		void add(std::string nm, double amp,double phase=0.0, double offset=0.0, double tbegin=SPulse::Start, double tend=SPulse::End);



//print functions
		void print(std::ostream &oo);
		void print();
};


std::ostream &operator<<(std::ostream &oo, Pulse &out);

END_BL_NAMESPACE



#endif



