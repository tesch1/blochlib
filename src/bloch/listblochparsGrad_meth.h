


#ifndef _listblochparsGrad_meth_cc_
#define _listblochparsGrad_meth_cc_

#include "utils/utils.h"
#include "bloch/blochParams.h"
#include "bloch/listblochpars.h"
#include "utils/constants.h"


BEGIN_BL_NAMESPACE


template<class Eng>
double &ListBlochParamsGrad<Eng>::offset()
{
	if(grads_ && apply_){
		off_=ListBlochParams::offset()+(grads_->offset()*gammaGauss());
	}else{
		off_=ListBlochParams::offset();
	}
	return off_;
}

template<class Eng>
void const ListBlochParamsGrad<Eng>::IterErr() 
{
	BLEXCEPTION(" itterated too far...")
}


template<class Eng>
ListBlochParamsGrad<Eng> &ListBlochParamsGrad<Eng>::operator=(const ListBlochParamsGrad<Eng> &rhs)
{
	if(this==&rhs) return *this;
	grads_=rhs.grads_;
	ListBlochParams::operator=(rhs);
	return *this;
}


template<class Eng>
void ListBlochParamsGrad<Eng>::operator++()
{
	ListBlochParams::operator++();
	if(grads_) ++(*grads_);
}



template<class Eng>
void ListBlochParamsGrad<Eng>::reset()
{
	if(grads_) grads_->reset();
	ListBlochParams::reset();
}



template<class Eng>
void ListBlochParamsGrad<Eng>::SetGrads(GradientGrid<Eng> &in) 
{
	if(in.size() != size()){
		BLEXCEPTION(" Gradient grid size must be the same as the number of spins...")
	}
	grads_=&in;
}


/* The entire list is saved like so...

----------
LISTBLOCHPARS   -->this is in ASCII
totalbyte Bo Temperature numberofspins
pars1 pars2 ... parsN
ENDLISTBLOCHPARS --->this is in ASCII
--------------
The first and last words allow us to find the parameters
if nested in a larger file of other data (binary or otherwise sooo DONOT use those words for
anything else...

*/

/*
bool ListBlochParams::write(std::fstream &oo)	//binary write
{
	if(!oo){
		cerr<<endl<<"Error: ListBlochPars::write(fstream)"<<endl;
		cerr<<" Bad file...."<<endl;
		return false;
	}
	oo<<endl<<"LISTBLOCHPARS\n";
	reset();
	long totalsize=0;
	while(*this){
		totalsize+=Vpars_(curpos_).binarySize();
		++(*this);
	}
	reset();
	oo.write(&totalsize, sizeof(long));
	oo.write(&Vpars_(0).Bo_, sizeof(double));
	oo.write(&Vpars_(0).Temperature_, sizeof(double));
	int nsp=Vpars_.size();
	oo.write(&nsp, sizeof(int));
	reset();
	while(*this){
		if(!Vpars_(curpos_).write(oo)) return false;
		++(*this);
	}
	reset();
	oo<<endl<<"ENDLISTBLOCHPARS\n";
	return true;

}

bool ListBlochParams::read(std::fstream &in)	//binart read
{
	if(!in){
		cerr<<endl<<"Error: ListBlochPars::read(fstream)"<<endl;
		cerr<<" Bad file...."<<endl;
		return false;
	}
	string temp="";
	char tm[1000];
	while(temp!="LISTBLOCHPARS"){
		in.getline(tm, 1000, '\n');
		temp=string(tm);
	}
	long totalsize=0;
	in.read(&totalsize, sizeof(long));
	in.read(&Vpars_(0).Bo_, sizeof(double));
	in.read(&Vpars_(0).Temperature_, sizeof(double));
	int nsp=1;
	in.read(&nsp, sizeof(int));
	Vpars_.resize(nsp);
	reset();
	long stepbin=0;
	while(*this ){
		Vpars_(curpos_).read(in);
		//stepbin+=Vpars_(curpos_).binarySize();
		++(*this);
	}
	reset();

	while(temp!="ENDLISTBLOCHPARS"){
		in.getline(tm, 1000, '\n');
		temp=string(tm);
	}
	return true;
}



std::ostream& operator<<(std::ostream &oo,ListBlochParams &out)
{
	oo<<"Bfield="<<out.Bo()<<endl;
	oo<<"Temperature="<<out.Temperature()<<endl;
	out.Vpars_.RecordSep="\n";
	oo<<out.Vpars_<<endl;
	return oo;
}
*/

END_BL_NAMESPACE


#endif



