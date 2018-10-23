

#ifndef _Prop_Reduce_h_
#define _Prop_Reduce_h_ 1


#include "blochlib.h"

/**** This class should be used as follows...

	PropReduce myred(base fact, log);
	myred.reduce();

	for(int i=0;i<base;++i){
		< generate the ind and fow props... >
	}
	for(int i=0;i<myred.maxBackReduce();++i){
		< generate the back props... >
	}

	Vector<matrix> myred.generateProps(ind, fow, bac);
***/

class PropReduce{
	private:
		BlochLib::Vector<BlochLib::Vector<int> > BackRed;
		BlochLib::Vector<int> BackName;

		BlochLib::Vector<BlochLib::Vector<int> > FowardRed;
		BlochLib::Vector<int> FowardName;

		BlochLib::Vector<BlochLib::Vector<int> > Shift1Red;
		BlochLib::Vector<int> Shift1Order;
		BlochLib::Vector<int> Shift1Name;

		BlochLib::Vector<BlochLib::Vector<int> > dat;
		BlochLib::Vector<BlochLib::Vector<int> > datCP; //copies dat for the datReduce
		BlochLib::Vector<int> datCPName; //the names for datCP

		int Mults, UseMe, baseTag,backTag, datTag, speTag,shiftTag;

		bool iteration(BlochLib::Vector<BlochLib::Vector<int> > &dat,
					   BlochLib::Vector<BlochLib::Vector<int> > &propRed,
					   BlochLib::Vector<BlochLib::Vector<int> > &subN,
					   BlochLib::Vector<int> &name,
					   int &numRep);

	public:

		int base, factor;
		 std::ostream *logf;

	//constructors
		PropReduce(){}
		PropReduce(int bas, int factor, std::ostream *oo=0);

		~PropReduce();

		void setParams(int bas, int fact,  std::ostream *oo=0);

	//runs through the 'dat' real sequences and pulls out the repetitions
	//that are the foward sequences---0,1--0,1,2--0,1,2,3...
		void fowardReduce();

	//runs through the 'dat' real sequences and pulls out the repetitions
	//that are the backwards sequences---9,10--8,9,10--7,8,9,10...
		void backReduce();

	//runs through the 'dat' real sequences and pulls out the repetitions
	//that are the shift1 sequences---
	// that is squences that are a simple time shift from the 
	// one or so previous
		void shift1Reduce();

	//des all the reductions in sequences
		void reduce();
	
	//does only the sequence reduction
		void sequenceReduce();
		
	//does only time shiftng reductions
		void timeShiftReduce();

		inline int maxBackReduce() const { 	return UseMe+1;	}
		inline int maxFowardReduce() const	{	return FowardRed.size();	}

	//these functions will create the propogators
	// from 3 input matrix lists..the first are the
	// individual propogators ("0","1","2"...)
	// the second the 'Foward' props ("0*1","0*1*2" ...)
	// the third the, 'Back' props ("7*8", "6*7*8"...)
	// the forth is the place to fill...
		void generateProps(BlochLib::Vector<BlochLib::matrix> &indiv,
		 					BlochLib::Vector<BlochLib::matrix> &Foward,
		 					BlochLib::Vector<BlochLib::matrix> &Back,
		 					BlochLib::Vector<BlochLib::matrix> &FillMe);

};


#endif


