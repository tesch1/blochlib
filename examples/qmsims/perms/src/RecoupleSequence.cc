

/*
this program calculates the tensorial components of various permutations
on the standard recoupling sequence Post-C7
*/

#ifndef _RecouplePermutions__cc_
#define _RecouplePermutions__cc_ 1


#include "RecoupleSequence.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;


/** Pulse Train **/
PulseTrain::PulseTrain(const Vector<double> &amp,
						const Vector<double> &phs,
						const Vector<double> &tim)
{
	if(amp.size()!= phs.size() != tim.size()+1){
		std::cerr<<std::endl<<"Error: PulseTrain(amps, phases, times)"<<std::endl;
		std::cerr<<" amp.size()==phases.size()==times.size()+1 "<<std::endl;
		exit(-1);
	}
	amplitudes=amp;
	phases=phs;
	times=tim;
}

void PulseTrain::operator=(const PulseTrain &rhs)
{
	if(&rhs==this) return;
	amplitudes=rhs.amplitudes;
	phases=rhs.phases;
	times=rhs.times;
}


matrix PulseTrain::Pulse(SpinSys &sys, std::string on, int which)
{
	matrix ret(sys.F0());
	RunTimeAssert(which<size());

	if(amplitude(which)==0.0)	return ret;

	for(int i=0;i<sys.size();++i){
		if(sys.symbol(i)==on){
			ret+=amplitudes(which)*cos(phases(which))*sys.Ix(i)+
			     amplitudes(which)*sin(phases(which))*sys.Iy(i);
		}
	}
	return ret;
}


/*****************************************************/
/** PulseTrain iterator **/
/*****************************************************/
matrix PulseTrainIterator::Pulse(SpinSys &sys, std::string on)
{
	matrix ret(sys.F0());

	if(amplitude()==0.0)	return ret;
	for(int i=0;i<sys.size();++i){
		if(sys.symbol(i)==on){
			ret+=amplitude()*cos(phase())*sys.Ix(i)+amplitude()*sin(phase())*sys.Iy(i);
		}
	}
	return ret;
}


void PulseTrainIterator::operator=(PulseTrainIterator &rhs)
{
	if(&rhs==this)	return;
	data=rhs.data;
	curpos=rhs.curpos;
}


/*****************************************************/
/*** Recouple ***/
/*****************************************************/

int Recouple::StringToRtype(std::string rtype)
{
	if(rtype=="R")	return R;
	if(rtype=="Rp") return Rp;
	if(rtype=="Cp")	return Cp;
	if(rtype=="C")	return C;
	if(rtype=="Cpp")	return Cpp;
	return 0;
}

std::string Recouple::RtypeToString(int rtype)
{
	if(rtype==R)	return "R";
	if(rtype==Rp) 	return "Rp";
	if(rtype==Cp)	return "Cp";
	if(rtype==C)	return "C";
	if(rtype==Cpp)	return "Cpp";
	return "";
}

void Recouple::init(double amp, double startp, double startt)
{
	amp_=amp;
	sphase_=startp;
	stime_=startt;
	double root;

	if(rtype_==R){
		amplitudes.resize(2*baseSym_,amp_);
		phases.resize(2*baseSym_);
		times.resize(2*baseSym_+1);
		int ct=0;
		times[0]=stime_;
		for(int i=0;i<2*baseSym_;i+=2){
			root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root;
			phases[i+1]=-root;

			times[i+1]=times[i]+1.0/amplitudes[i];
			times[i+2]=times[i+1]+1.0/amplitudes[i+1];
			++ct;
		}

	}else if(rtype_==Rp){
		amplitudes.resize(4*baseSym_,amp_);
		phases.resize(4*baseSym_);
		times.resize(4*baseSym_+1);
		int ct=0;
		times[0]=stime_;
		for(int i=0;i<4*baseSym_;i+=4){
			root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root;
			phases[i+1]=root+Pi;
			phases[i+2]=2.*Pi-root;
			phases[i+3]=2.*Pi-root-Pi;


			times[i+1]=times[i]+0.25*1.0/amplitudes[i];
			times[i+2]=times[i+1]+0.75*1.0/amplitudes[i+1];
			times[i+3]=times[i+2]+0.25*1.0/amplitudes[i+2];
			times[i+4]=times[i+3]+0.75*1.0/amplitudes[i+3];
			++ct;
		}


	}else if(rtype_==C){
		amplitudes.resize(2*baseSym_,amp_);
		phases.resize(2*baseSym_);
		times.resize(2*baseSym_+1);
		int ct=0;
		times[0]=stime_;
		for(int i=0;i<2*baseSym_;i+=2){
			root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root;
			phases[i+1]=root+PI;

			times[i+1]=times[i]+1.0/amplitudes(i);
			times[i+2]=times[i+1]+1.0/amplitudes(i+1);
			++ct;
		}
	}else if(rtype_==Cp){
		amplitudes.resize(3*baseSym_,amp_);
		phases.resize(3*baseSym_);
		times.resize(3*baseSym_+1);
		int ct=0;
		times[0]=stime_;
		for(int i=0;i<3*baseSym_;i+=3){
			root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root ;
			phases[i+1]=root+PI;
			phases[i+2]=root;

			times[i+1]=times[i]+0.25*1.0/amplitudes(i);
			times[i+2]=times[i+1]+1.0/amplitudes(i+1);
			times[i+3]=times[i+2]+0.75*1.0/amplitudes(i+2);
			++ct;
		}

	}else if(rtype_==Cpp){
		int i=0;
		/*
		For 'odd' BaseSym_ we need to 'chop' the
		middle pulse putting the 'bar' set on the frount
		 and the non-bar on the end
		 for baseSym==7, if we define a unit as
		 (90, p_n)(270, p_n+pi)=N (where n=0...7)
		 and
		 (90, p_n+pi)(270, p_n)=Nbar (where n=0...7)
		 then the resulting 'permutation' is
		 3bar,2,2Bar,1,1Bar,0,0Bar,6,6bar,5,5Bar,4,4Bar,3

		 for baseSym_==8 we would have
		 3,3bar,2,2bar,1,1bar,0,0bar,7,7bar,6,6bar,5,5bar,4,4bar

		So out 'generic' Cpp sequence is then
		if baseSym is Even
		Cpp-->n=[BaseSym/2-1...0]+
		        [BaseSym-1...BaseSym/2]
		if BaseSym is Odd
		Cpp-->n=[((BaseSym+1)/2)-1]bar +
		        [((BaseSym+1)/2)-2...0]+
		        [BaseSym-1...(BaseSym+1)/2)]+
		        [((BaseSym+1)/2)-1]
		*/
		if(baseSym_%2==0){ //even
			//our Sub unit here is (90, 360, 270)=(N, Nbar)
			//[BaseSym/2-1...0]
			amplitudes.resize(3*baseSym_,amp_);
			phases.resize(3*baseSym_);
			times.resize(3*baseSym_+1);
			times[0]=stime_;
			for(int ct=(baseSym_/2)-1; ct>=0; --ct)
			{
				root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
				phases[i]=root ;
				phases[i+1]=root+PI;
				phases[i+2]=root;

				times[i+1]=times[i]+0.25*1.0/amplitudes(i);
				times[i+2]=times[i+1]+1.0/amplitudes(i+1);
				times[i+3]=times[i+2]+0.75*1.0/amplitudes(i+2);
				i+=3;
			}

			//[BaseSym-1...BaseSym/2]
			for(int ct=baseSym_-1; ct>=baseSym_/2; --ct)
			{
				root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
				phases[i]=root ;
				phases[i+1]=root+PI;
				phases[i+2]=root;

				times[i+1]=times[i]+0.25*1.0/amplitudes(i);
				times[i+2]=times[i+1]+1.0/amplitudes(i+1);
				times[i+3]=times[i+2]+0.75*1.0/amplitudes(i+2);
				i+=3;
			}
		}else{ //odd
		//we add a 'extra' element due to the 'split' pulse
			amplitudes.resize(3*baseSym_+1,amp_);
			phases.resize(3*baseSym_+1);
			times.resize(3*baseSym_+2);
			times[0]=stime_;

			//our Sub unit here is (90, 270)=(N)
			//[((BaseSym+1)/2)-1]bar
			root=((double(baseSym_+1)/2.0)-1.0)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root ;
			phases[i+1]=root+PI;

			times[i+1]=times[i]+0.25*1.0/amplitudes(i);
			times[i+2]=times[i+1]+0.75*1.0/amplitudes(i+1);
		    i+=2;

		    //[((BaseSym+1)/2)-2...0]
			for(int ct=((baseSym_+1)/2)-2; ct>=0; --ct)
			{
				root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
				phases[i]=root ;
				phases[i+1]=root+PI;
				phases[i+2]=root;

				times[i+1]=times[i]+0.25*1.0/amplitudes(i);
				times[i+2]=times[i+1]+1.0/amplitudes(i+1);
				times[i+3]=times[i+2]+0.75*1.0/amplitudes(i+2);
				i+=3;
			}
		    //[BaseSym-1...(BaseSym+1)/2)]+
			for(int ct=baseSym_-1; ct>=(baseSym_+1)/2; --ct)
			{
				root=double(ct)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
				phases[i]=root ;
				phases[i+1]=root+PI;
				phases[i+2]=root;

				times[i+1]=times[i]+0.25*1.0/amplitudes(i);
				times[i+2]=times[i+1]+1.0/amplitudes(i+1);
				times[i+3]=times[i+2]+0.75*1.0/amplitudes(i+2);
				i+=3;
			}

			//[((BaseSym+1)/2)-1]
			root=((double(baseSym_+1)/2.0)-1.0)*(2.*Pi*double(fractSym_)/double(baseSym_))+sphase_;
			phases[i]=root+Pi;
			phases[i+1]=root;

			times[i+1]=times[i]+0.25*1.0/amplitudes(i);
			times[i+2]=times[i+1]+0.75*1.0/amplitudes(i+1);
		    i+=2;
		}
	}else{
		std::cerr<<std::endl<<"Error: Recouple::init()"<<std::endl;
		std::cerr<<" Valid entries are 'R', 'Rp', 'C', 'Cp', and 'Cpp' for 'rtype'"<<std::endl;
		exit(-1);
	}
}

void Recouple::read(Parameters &pset, std::string sec)
{
	pset.addSection(sec);
	amp_=pset.getParamD("amplitude", sec);
	sphase_=pset.getParamD("startphase", sec, false);
	stime_=pset.getParamD("starttime", sec, false);
	rtype_=StringToRtype(pset.getParamS("rtype", sec));
	baseSym_=pset.getParamI("baseSym", sec);
	fractSym_=pset.getParamI("fractSym", sec);

	init(amp_, sphase_, stime_);
}

Recouple::Recouple(int rtype, int baseSym, int fractSym, double amp, double phaseoff, double startt)
	: rtype_(rtype), baseSym_(baseSym), fractSym_(fractSym)
{	init(amp, phaseoff, startt);	}

Recouple::Recouple(std::string rtype, int baseSym, int fractSym, double amp, double phaseoff, double startt)
	: rtype_(StringToRtype(rtype)), baseSym_(baseSym), fractSym_(fractSym)
{	init(amp, phaseoff, startt);	}

Recouple::Recouple(Parameters &pset, std::string sec)
{	read(pset,sec);	}

void Recouple::operator=(const Recouple &rhs)
{
	rtype_=rhs.rtype_;
	baseSym_=rhs.baseSym_;
	fractSym_=rhs.fractSym_;
	amp_=rhs.amp_;
	sphase_=rhs.sphase_;
	stime_=rhs.stime_;
	PulseTrain::operator=(rhs);
}

void Recouple::setRecoupleType(int rtype)
{
	rtype_=rtype;
	init(amp_, sphase_, stime_);
}

void Recouple::setBaseSym(int baseSym)
{
	baseSym_=baseSym;
	init(amp_, sphase_, stime_);
}

void Recouple::setFractionalSym(int fractSym)
{
	fractSym_=fractSym;
	init(amp_, sphase_, stime_);
}


void Recouple::setSequence(int rtype, int baseSym, int fractSym)
{
	rtype_=rtype;
	baseSym_=baseSym;
	fractSym_=fractSym;
	init(amp_, sphase_, stime_);
}


void Recouple::setAmplitude(double amp)
{	init(amp, sphase_, stime_);	}

void Recouple::setStartPhase(double sphase)
{	init(amp_, sphase, stime_);	}

void Recouple::setStartTime(double stime)
{	init(amp_, sphase_, stime);	}

void Recouple::setPulseParams(double amp, double sphase, double stime)
{	init(amp, sphase, stime);	}



ostream &operator<<(ostream &oo, const Recouple &out)
{
	oo<<"Recoupling Sequence::"
	  <<out.RtypeToString(out.rtype_)
	  <<",n="<<out.fractSym_
	  <<",N="<<out.baseSym_
	  <<" --start phase:"<<out.phases[0]
	  <<" amplitude: "<<out.amplitudes[0]
	  <<" start time: "<<out.times[0];

	return oo;
}



/*****************************************************/
/*** Recouple train ***/
/*****************************************************/
RecoupleTrain::RecoupleTrain(int len)
{
	train.resize(len);
}

RecoupleTrain::RecoupleTrain(const Vector<Recouple> &inp ):
	train(inp)
{}

RecoupleTrain::RecoupleTrain(Parameters &pset, std::string sec)
{
	int le=pset.getParamI("trainlength", sec);
	std::string sbase=pset.getParamS("secbase", sec, false, "sec");
	train.resize(le);
	Parameters pset2(pset.section(sec));
	for(int i=0;i<le;++i){
		std::string secn=sbase+itost(i);
		train[i]=Recouple(pset2, secn);
	}
}

void RecoupleTrain::setBeginTime(double inbtime)
{
	double lastbtime=inbtime;
	for(int i=0;i<train.size();++i){
		train[i].setStartTime(lastbtime);
		//for(int j=0;j<train[i].size()+1;++j){
			//if(j==0)	train[i].times[j]=lastbtime;
		//	train[i].times[j]+=lastbtime-train[i].beginTime();
		//}
		lastbtime=train[i].endTime();
	}
}

ostream &operator<<(ostream &oo,const RecoupleTrain &out)
{
	for(int i=0;i<out.size();++i){
		oo<<out[i]<<endl;
	}
	return oo;
}


/*****************************************************/
/*** RecoupleTrain Iterator **/
/*****************************************************/
void RecoupleTrainIter::reset(){
	curpos=0;
	if(data->size()>0){
		Recouple::iterator moo(data->element(0));
		curRecouple=moo;
	}
}


RecoupleTrainIter::RecoupleTrainIter(RecoupleTrain &in):
	data(&in)
{	reset();	}


void RecoupleTrainIter::operator=(RecoupleTrainIter &rhs)
{
	if(&rhs==this) return;
	data=rhs.data;
	curpos=rhs.curpos;
	curRecouple=rhs.curRecouple;
}

void RecoupleTrainIter::operator++()
{
	++curRecouple;
	if(!curRecouple){
		++curpos;
		Recouple::iterator moo(data->element(curpos));
		curRecouple=moo;
	}
}

void RecoupleTrainIter::operator++(int)
{	operator++();	}


matrix RecoupleTrainIter::Pulse(SpinSys &in, std::string on)
{	return curRecouple.Pulse(in, on);	}



#endif
