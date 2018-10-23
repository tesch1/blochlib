

#include "blochlib.h"
#include "propreduce.h"

using namespace BlochLib;
using namespace std;

PropReduce::PropReduce(int bas, int fac, ostream *oo)
{	setParams(bas, fac, oo);	UseMe=0; Mults=0;}

PropReduce::~PropReduce(){	logf=NULL;	}

void PropReduce::setParams(int bas, int fac, ostream *oo)
{
	//find the greatest common divisor...
	int u = abs(bas);
    int v = abs(fac);
    int q,t;
   	while(v){
       q = int(floor( double(u)/double(v) ));
       t = u - v*q;
       u = v;
       v = t;
   }

    base=bas/u;
	factor=fac/u;
	logf=oo;


/*
the oringal sequence lists
if (base=9, and factor=7)
Self Reduction: 900000=0 1 2 3 4 5 6
Self Reduction: 1800000=7 8 0 1 2 3 4
Self Reduction: 2700000=5 6 7 8 0 1 2
Self Reduction: 3600000=3 4 5 6 7 8 0
Self Reduction: 4500000=1 2 3 4 5 6 7
Self Reduction: 5400000=8 0 1 2 3 4 5
Self Reduction: 6300000=6 7 8 0 1 2 3
Self Reduction: 7200000=4 5 6 7 8 0 1
Self Reduction: 8100000=2 3 4 5 6 7 8

if(base=7, factor=9)
Self Reduction: 900000=0 1 2 3 4 5 6 0 1
Self Reduction: 1800000=2 3 4 5 6 0 1 2 3
Self Reduction: 2700000=4 5 6 0 1 2 3 4 5
Self Reduction: 3600000=6 0 1 2 3 4 5 6 0
Self Reduction: 4500000=1 2 3 4 5 6 0 1 2
Self Reduction: 5400000=3 4 5 6 0 1 2 3 4
Self Reduction: 6300000=5 6 0 1 2 3 4 5 6

*/
	int ct=0, ct2=0;
//give the first datCPName a value
	dat.resize(base, Vector<int>(factor));
	datTag=100000*std::max(base, factor);
    for(int i=0;i<base*factor;++i){
		dat[ct][ct2]=i%base;
		++ct2;
		if(ct2>=factor){
			if(logf) *logf<<"Original SEQ: "<<"="<<dat[ct]<<std::endl;
			++ct;
			ct2=0;
		}
	}

//the size of the foward possible reductions
// have a max length of min(base,FACTOR)
// -if factor>base then dat[i].size() will be smaller then some of
// these reduction and will be meaningless nothing will be reduced
// -if base>factor then we will run out of integer marks
//
	int maxRed=std::min(base, factor);
/*
FOwards look like
if (base=9, and factor=7)
Foward reduction: 700=0 1
Foward reduction: 1400=0 1 2
Foward reduction: 2100=0 1 2 3
Foward reduction: 2800=0 1 2 3 4
Foward reduction: 3500=0 1 2 3 4 5
Foward reduction: 4200=0 1 2 3 4 5 6

if (base=7, and factor=9)
SAME AS ABOVE!

*/

	FowardName.resize(maxRed-1);
	FowardRed.resize(maxRed-1);

//create those lists...which should
	baseTag=100*maxRed;
	for(int i=1;i<maxRed;++i){
		FowardName[i-1]=i*baseTag;
		FowardRed[i-1].resize(i+1);
		for(int j=0;j<=i;++j)	FowardRed[i-1][j]=j;

		if(logf) *logf<<"Foward reduction: "<<FowardName[i-1]<<"="<<FowardRed[i-1]<<std::endl;
	}


/*

we need to generate ALL the backwards
becuase we can play tricks with some time shifting
later on

Backwards look like
if(base=9 and factor=7)
Back reduction: -700=7 8
Back reduction: -1400=6 7 8
Back reduction: -2100=5 6 7 8
Back reduction: -2800=4 5 6 7 8
Back reduction: -3500=3 4 5 6 7 8
Back reduction: -4200=2 3 4 5 6 7 8
Back reduction: -4200=1 2 3 4 5 6 7 8
Back reduction: -4200=0 1 2 3 4 5 6 7 8

if (base=7, and factor=9)
Back reduction: -700=5 6
Back reduction: -1400=4 5 6
Back reduction: -2100=3 4 5 6
Back reduction: -2800=2 3 4 5 6
Back reduction: -3500=1 2 3 4 5 6
Back reduction: -4200=0 1 2 3 4 5 6

*/
	backTag=100*base;
	BackName.resize(base-1);
	BackRed.resize(base-1);
	for(int i=1;i<maxRed;++i){
		BackName[i-1]=-i*baseTag;
		BackRed[i-1].resize(i+1);
		for(int j=base-i-1, k=0;j<base;++j,++k){	BackRed[i-1][k]=j;}
		if(logf) *logf<<"Back reduction: "<<BackName[i-1]<<"="<<BackRed[i-1]<<std::endl;
	}

/*
these 'special' ones are simply a combo of
the self sequences but with a few time reverse
propogators at the begning and end

like so
  U[0]' data[0] U[data[0][base-1]]'
  U[1]' data[1] U[data[1][base-1]]'
  ...

if(base=7 and factor=9)
NOT NESSESARY

if(base=9 and factor=7)
1 time shift Reduction: 1 2 3 4 5 6 7 -- U(0)' self(0) U(7)
1 time shift Reduction: 8 0 1 2 3 4 5 -- U(7)' self(1) U(5)
1 time shift Reduction: 6 7 8 0 1 2 3 -- U(5)' self(2) U(3)
1 time shift Reduction: 4 5 6 7 8 0 1 -- U(3)' self(3) U(1)
1 time shift Reduction: 2 3 4 5 6 7 8 -- U(1)' self(4) U(8)
1 time shift Reduction: 0 1 2 3 4 5 6 -- U(8)' self(5) U(6)
1 time shift Reduction: 7 8 0 1 2 3 4 -- U(6)' self(6) U(4)
1 time shift Reduction: 5 6 7 8 0 1 2 -- U(4)' self(7) U(2)
1 time shift Reduction: 3 4 5 6 7 8 0 -- U(2)' self(8) U(0)

we only need to do this IF base>factor otherwise
the foward and backwards take care of most everything

this carries an extra list becuase we can only use these
if it comes AFTER we have already generated the self(n)
otherwise it is useless

*/
	if(base>factor)
	{
		Shift1Order.resize(base);
		Shift1Red.resize(base, Vector<int>(factor));
		shiftTag=20000*std::max(base, factor);
		Shift1Name.resize(base);
		for(int i=0;i<base;++i)
		{
			//the name tag
			Shift1Name[i]=(i+1)*speTag;

			//the order
			Shift1Order[i]=i;

			//the end
			Shift1Red[i][dat[i].size()-1]=(dat[i][dat[i].size()-1]==(base-1))?0:(dat[i][dat[i].size()-1]+1);

			//file in the rest
			for(int k=0;k<dat[i].size()-1;++k) Shift1Red[i][k]=dat[i][k+1];

			if(logf) *logf<<"1 time shift reduction: "<<Shift1Name[i]<<"="<<Shift1Red[i]<<std::endl;
		}
	}




}


bool PropReduce::iteration(Vector<Vector<int> > &dat,
						   Vector<Vector<int> > &propRed,
				           Vector<Vector<int> > &subN,
				           Vector<int> &name,
				           int &numRep)
{
	//loops to find the matches
	bool gotanyTot=false;

//loop thought the basic data seq's
	for(int i=0;i<dat.size();++i){
		bool gotany=false;
		Vector<int> curU;
	//loop through each chunk within a sequence
		for(int M=0;M<dat[i].size();++M){
			bool got=false;
			int p=0;
		//reduce going backwards as we want to catch the largest chunks first
			for(p=subN.size()-1;p>=0;--p){
				if((subN[p].size()+M)<=dat[i].size() ){
					//cout<<p<<" "<<Range(M,M+subN[p].size()-1)<<endl;
					if(subN[p]==dat[i](Range(M,M+subN[p].size()-1))){
						got=true;
						numRep++;
						break;
					}
				}
			}

			//we got a match for a valid reduction
			if(got){
			//fill the new set up to the 'M' stop mark above
				for(int k=0;k<M;++k)	curU.push_back(dat[i][k]);
			//replace the reduced chunk with the name
				curU.push_back(name[p]);
			//fill the new set up from the stop mark to the end
				for(int k=subN[p].size()+M;k<dat[i].size();++k){
					curU.push_back(dat[i][k]);
				}
				propRed[i]=curU;
				gotany=true;
				break;
			}
		}
	//if we did not get any copy the original seq back
		if(!gotany){
			for(int k=0;k<dat[i].size();++k)	curU.push_back(dat[i][k]);
			propRed[i]=(curU);
		}else{
			gotanyTot=true;
		}
	}
	return gotanyTot;
}


/*** foward reductions... ***/
void PropReduce::fowardReduce()
{
	int dum=0;
	Vector<Vector<int> > propRed(base,Vector<int>(0));
	while(iteration(dat, propRed, FowardRed, FowardName,dum)){
		dat=propRed;
	}

	int multi=0;
	for(int i=0;i<dat.size();++i){
		if(logf) *logf<<"Sequence "<<i<<": "<<dat[i]<<endl;
		multi+=dat[i].size();
	}
	if(logf) *logf<<" After Foward Reduction...Number of multiplications: "<<multi<<endl<<endl;
}

/***Back Reductions ***/
// the back reductions we do note get for free
// (like the forward ones whcih we have to calc
// from the exp(H) operation), so the number
// of back reductions used depends on the total multiplication
// saveings...so we need to go through the entire loops of
// back reductions...
void PropReduce::backReduce()
{
	Vector<Vector<int> > propRed(base,Vector<int>(0));
	Vector<Vector<int> > holdDat(dat.size());
	for(int i=0;i<dat.size();++i)	holdDat[i]=dat[i];

	Mults=1000000;
	int multi=0, dum=0;
	UseMe=0;
	Vector<Vector<int> > curBack;
	Vector<int> curName;
	for(int k=0;k<BackRed.size();++k){
		if(logf) *logf<<" Number of 'Back Reductions': "<<k<<endl;
		curBack=BackRed(Range(0,k));
		curName=BackName(Range(0,k));

		for(int i=0;i<dat.size();++i)	dat[i]=holdDat[i];


		while(iteration(dat, propRed, curBack, curName,dum)){
			dat=propRed;

			multi=curBack.size();
			for(int j=0;j<dat.size();++j){
				multi+=dat[j].size();
			}

			if(Mults>multi){
				UseMe=k;
				Mults=multi;
			}
		}
		for(int j=0;j<dat.size();++j){
			if(logf) *logf<<"Sequence "<<j<<": "<<dat[j]<<std::endl;
		}
		if(logf) *logf<<" After Back Reduction...Number of multiplications: "<<multi<<std::endl<<std::endl;
	}

//need to 'regen' the best one for displaying
	if(logf) *logf<<" Number of 'Back Reductions': "<<UseMe<<std::endl;
	curBack=BackRed(Range(0,UseMe));
	curName=BackName(Range(0,UseMe));

	for(int i=0;i<dat.size();++i)	dat[i]=holdDat[i];

	Vector<int> BackNeedToGen;
	while(iteration(dat, propRed, curBack, curName,dum)){
		dat=propRed;
	}
}

/*** Special Reductions ***/
void PropReduce::shift1Reduce()
{
	Vector<Vector<int> > propRed(base,Vector<int>(0));
	int rep=0;
	//loops to find the matches
	bool gotanyTot=true;

//loop thought the basic data seq's
	while(gotanyTot){
		gotanyTot=false; //need to reset to false...
		for(int i=0;i<dat.size();++i){
			bool gotany=false;
			Vector<int> curU;
		//loop through each chunk within a sequence
			for(int M=0;M<dat[i].size();++M){
				bool got=false;
				int p=0;
				
			//all special reduces are the same size
			//as the dat  but we need to record the ORDER in which things need
			// to be computed, as, this will be best suited to 
			// calculated group propogators out of order...
				for(p=0;p<Shift1Red.size();++p){
					if(Shift1Red[p]==dat[i]){ // && Shift1Order[p]<i
						got=true; break;
					}
				}
			
				//we got a match for a valid reduction
				if(got){
				//replace the reduced chunk with the name
					curU.push_back(Shift1Name[p]);
					propRed[i]=curU;
					Shift1Order[i] = p;
					gotany=true;
					break;
				}
			}
		//if we did not get any, copy the original seq back
			if(!gotany){
				for(int k=0;k<dat[i].size();++k)	curU.push_back(dat[i][k]);
				propRed[i]=(curU);
			}else{
				gotanyTot=true;
			}
		}
		dat=propRed;
	}

	Vector<Vector<int> > curBack=BackRed(Range(0, UseMe));
	int multi=curBack.size();
	for(int i=0;i<dat.size();++i){
		if(logf) *logf<<"Sequence "<<i<<": "<<dat[i]<<std::endl;
		multi+=dat[i].size();
	}
	//need to add 2 for every speacil one used as it takes 2 additional
	// muls to make the time shifted ones
	multi+=2*rep;
	int ttt=Mults-multi; //savings for 'specials'
	Mults-=ttt;
	if(logf) *logf<<" After First Time Shift Reductions...Number of multiplications: "<<Mults<<std::endl<<std::endl;
}


void PropReduce::reduce()
{

	//can only do this if factor>base
//	if(factor>base) selfReduce();
	fowardReduce();
	backReduce();
	if(base>factor)	shift1Reduce();
	if(logf) *logf<<endl<<" The Best Reduction is for using "<<UseMe+1<<" Back Reductions"<<std::endl;
	if(logf) *logf<<" For a grand total of "<<Mults<<" multipications"<<std::endl;
	if(logf) *logf<<" The total Sequence...."<<std::endl;
	for(int j=0;j<dat.size();++j){
		if(logf) *logf<<"Sequence "<<j<<": "<<dat[j]<<std::endl;
	}
}

//does only the sequence reduction
void PropReduce::sequenceReduce()
{
	if(logf) *logf << "**sequenceReduce**"<<endl;
	fowardReduce();
	backReduce();
	if(logf){
		*logf<<endl<<" The Best Reduction is for using "<<UseMe+1<<" Back Reductions"<<std::endl;
		*logf<<" For a grand total of "<<Mults<<" multipications"<<std::endl;
		*logf<<" The total Sequence...."<<std::endl;
		for(int j=0;j<dat.size();++j){
			*logf<<"Sequence "<<j<<": "<<dat[j]<<std::endl;
		}
	}
	
}

//does only time shiftng reductions
void PropReduce::timeShiftReduce()
{
	if(logf) *logf << "**timeShiftReduce**"<<endl;
	shift1Reduce();
	if(logf){
		*logf<<" For a grand total of "<<Mults<<" multipications"<<std::endl;
		*logf<<" The total Sequence...."<<std::endl;
		for(int j=0;j<dat.size();++j){
			*logf<<"Sequence "<<j<<": "<<dat[j]<<" ORDER:: "<<Shift1Order[j]<<std::endl;
		}
	}
	
}


//these functions will create the propogators
// from 3 input matrix lists..the first are the
// individual propogators ("0","1","2"...)
// the second the 'Foward' props ("0*1","0*1*2" ...)
// the third the, 'Back' props ("7*8", "6*7*8"...)
// the forth is the place to fill...
void PropReduce::
		generateProps(Vector<matrix> &indiv,
					  Vector<matrix> &Foward,
					  Vector<matrix> &Back,
					  Vector<matrix> &FillMe)
{
	if(indiv.size() != base){
		std::cerr<<"PropReduce::generateProps()"<<endl;
		std::cerr<<" Individual Matricies must have length 'base'"<<endl;
		exit(1);
	}

	if(Foward.size() != base-1){
		std::cerr<<"PropReduce::generateProps()"<<endl;
		std::cerr<<" Foward Matricies must have length 'base-1'"<<endl;
		exit(1);
	}

	if(FillMe.size() != base){
		std::cerr<<"PropReduce::generateProps()"<<endl;
		std::cerr<<" FillMe Matricies must have length 'base'"<<endl;
		exit(1);
	}

	if(Back.size() != UseMe+1){
		std::cerr<<"PropReduce::generateProps()"<<endl;
		std::cerr<<" Back Matricies must have the proper "<<endl;
		std::cerr<<" length from 'maxBackReduce()'"<<endl;
		exit(1);
	}

	for(int i=0;i<dat.size();++i){
		//cout<<"Propogator: "<<i<<" ||: ";
		for(int j=0;j<dat[i].size();j++){
			//cout<<" INdex: "<<dat[i][j];
			if(j==0){
				if(dat[i][j]>=baseTag && dat[i][j]!=speTag){
				//	cout<<" F reduced: "<<dat[i][j]/baseTag -1;
					FillMe[i]=Foward[dat[i][j]/baseTag -1];
				}else if(dat[i][j]<0){
				//	cout<<" B reduced: "<<-dat[i][j]/baseTag -1;
					FillMe[i]=Back[-dat[i][j]/baseTag -1];
				}else if(dat[i][j]==speTag){
				//	cout<<" S reduced: "<<base-2;
					FillMe[i]=adjoint(indiv[base-1])*Foward[base-2]*adjoint(indiv[0]);
				}else{
				//	cout<<" N reduced: "<<dat[i][j];
					FillMe[i]=indiv[dat[i][j]];
				}
			}else{
				if(dat[i][j]>=baseTag && dat[i][j]!=speTag){
				//	cout<<" F reduced: "<<dat[i][j]/baseTag -1;
					FillMe[i]=Foward[dat[i][j]/baseTag -1]*FillMe[i];
				}else if(dat[i][j]<0){
				//	cout<<" B reduced: "<<-dat[i][j]/baseTag -1;
					FillMe[i]=Back[-dat[i][j]/baseTag -1]*FillMe[i];
				}else if(dat[i][j]==speTag){
				//	cout<<" S reduced: "<<base-2;
					FillMe[i]=adjoint(indiv[base-1])*Foward[base-2]*adjoint(indiv[0])*FillMe[i];
				}else{
				//	cout<<" N reduced: "<<dat[i][j];
					FillMe[i]=indiv[dat[i][j]]*FillMe[i];
				}
			}
		}
	//	cout<<endl;
	}
}




