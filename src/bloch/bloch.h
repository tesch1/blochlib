

#ifndef _bloch_h_
#define _bloch_h_ 1


#include "mpi/mpi_config.h"

#include "utils/constants.h"
#include "container/grids/coords.h"
#include "container/matrix/matrix.h"

#include "bloch/listblochpars.h"
#include "bloch/pulse.h"
#include "bloch/bloch_basic.h"
#include "bloch/bloch_interac.h"
#include <string>

BEGIN_BL_NAMESPACE



/*
	Explanations::
		--There are 4 template parameters for a 'Bloch' function driver
			1) ParamSet-->the parameter set (contains the back end driver for gradients, grid, and spin type)
			2) Pulses-->A pulse or no pulse...
			3) Interactions_t --> OTHER interactions NOT included in the ParamSet (the detail of this class
				is explained below)
			4) T--> the precision of the calculation (default=double)
		--There are also 2 'pre-expressed' versions of 'Bloch' where only the ParamSet is needed
			1) A 'Basic' No pulse version (no external interactions)
			2) A 'Basic' Pulse version (no external interactions)
		--There are 2 'pre-expressed' versions of 'Bloch' that require a ParamSet AND external Intereactions
			1) One with no pulses
			2) one with pulses

	the classes above are for simple bloch sims (ones with the _basic_ interactions (offsets, T1, T2))
   the next series are for additional interactions (suseptibility, radiation damping, etc)

   a new class is set up to act as the container for the interactions such that each series of 'interactions'
   can be passed as one template parameter to the bloch functions
*/


//pre declaration of the Blcoh with Pulses
template<class ParamSet, class Pulses=NullPulse, class Interactions_t=NoInteractions, class T=double >
class Bloch: public BlochBasic<ParamSet>{};

//foward delcaration for the 'operator=()' below
template<class ParamSet>
class Bloch<ParamSet, Pulse,NoInteractions, double>;

//This is the special class that contains NO pulses...
template<class ParamSet >
class Bloch<ParamSet, NullPulse, NoInteractions, double>:
	public BlochBasic<ParamSet>
{

	public:
		Bloch():
			BlochBasic<ParamSet>() {}

		Bloch(int N):
			BlochBasic<ParamSet>(N){}

		Bloch(std::string sp):
			BlochBasic<ParamSet>(sp) {}

		Bloch(int N, std::string sp):
			BlochBasic<ParamSet>(N,sp) {}

		Bloch(int N, std::string sp, ParamSet &in):
			BlochBasic<ParamSet>(N,sp,in){}

		Bloch(ParamSet &in):
			BlochBasic<ParamSet>(in){}

		Bloch(Bloch &in):
			BlochBasic<ParamSet>(in){}

		template<class Pulse_T, class Inters_T>
		Bloch(Bloch<ParamSet, Pulse_T, Inters_T> &in):
			BlochBasic<ParamSet>(in)
		{}

		~Bloch(){}

		typedef typename BlochBasic<ParamSet>::iterator iterator;

		void operator=(const Bloch &rhs)
		{
			if(this==&rhs) return;
			BlochBasic<ParamSet>::operator=(rhs);
		}

		template<class otherInters>
		Bloch &operator=(const Bloch<ParamSet, Pulse, otherInters, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			return *this;
		}

};



//This is the special class that contains WITH pulses...
// Uses the class "Pulse" found in "Pulse.h"

template<class ParamSet>
class Bloch<ParamSet, Pulse,NoInteractions, double>: public BlochBasic<ParamSet>{
	private:
		Pulse *plist;

	public:

		typedef typename BlochBasic<ParamSet>::iterator iterator;


		Bloch():
			BlochBasic<ParamSet>(), plist(0) {}

		Bloch(int N):
			BlochBasic<ParamSet>(N), plist(0) {}

		Bloch(std::string sp):
			BlochBasic<ParamSet>(sp), plist(0) {}

		Bloch(int N, std::string sp):
			BlochBasic<ParamSet>(N,sp), plist(0) {}

		Bloch(int N, std::string sp, ParamSet &in):
			BlochBasic<ParamSet>(N,sp, in), plist(0) {}

		Bloch(int N, std::string sp, ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(N,sp,in)
		{
			plist=&pl;
		}

		Bloch(int N, ParamSet &in):
			BlochBasic<ParamSet>(N,in), plist(0) {}

		Bloch(int N, ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(N,in)
		{
			plist=&pl;
		}

		Bloch(ParamSet &in):
			BlochBasic<ParamSet>(in) {}

		Bloch(ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(in)
		{
			plist=&pl;
		}

		Bloch(Bloch &in):
			BlochBasic<ParamSet>(in), plist(in.plist) {}


		Bloch(Bloch<ParamSet, NullPulse,NoInteractions, double> &in):
			BlochBasic<ParamSet>(in), plist(0) {}

		template<class  OtherInters_T>
		Bloch(Bloch<ParamSet, NullPulse,OtherInters_T, double> &in):
			BlochBasic<ParamSet>(in), plist(0) {}

		template<class  OtherInters_T>
		Bloch(Bloch<ParamSet, Pulse,OtherInters_T, double> &in):
			BlochBasic<ParamSet>(in), plist(in.plist) {}

		~Bloch()
		{
			plist=0;
		}

		inline Pulse *Pulses()		const {	return plist;	}
		inline void setPulses(Pulse &in){	plist=&in;	}
		inline void setPulses(Pulse *in){	plist=in;	}

		Bloch &operator=(const Bloch &rhs){
			if(this==&rhs) return *this;
			BlochBasic<ParamSet>::operator=(rhs);
			plist=rhs.plist;
			return *this;
		}

		template<class otherInters>
		Bloch &operator=(const Bloch<ParamSet, NullPulse,otherInters, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			plist=0;
			return *this;
		}
/*
		inline void realPulse(double t, coord<> &iny, coord<> &dydt, const std::string &sym)
		{
			if(plist)
			{
				coord<> B1=(plist->CoordPulse(t,sym));
				dydt.x()-=(B1.z()*iny.y()-B1.y()*iny.z());
				dydt.x()-=(B1.x()*iny.z()-B1.z()*iny.x());
				dydt.x()-=(B1.y()*iny.x()-B1.x()*iny.y());
			}
		}
*/
		inline coord<> realPulse(double t, coord<> &iny, coord<> &dydt, const std::string &sym)
		{
			if(plist)
			{
				return  cross(iny, plist->CoordPulse(t,sym));
						/*coord<>(
						-(B1.z()*iny.y()-B1.y()*iny.z()),
						-(B1.x()*iny.z()-B1.z()*iny.x()),
						-(B1.y()*iny.x()-B1.x()*iny.y())
					);*/
			}
			return ZeroType<coord<> >::zero();
		}

		template<class Piter>
		inline rmatrix pulseJacobian(double t,Piter *onit)
		{
			static rmatrix J(3,3,0);
			static coord<> B1;
			if(plist)
			{
				B1=(plist->CoordPulse(t,onit->symbol()));
				J(0,0)=0.0; 	J(0,1)=-B1.z(); 	J(0,2)=B1.y();
				J(1,0)=B1.z();	J(1,1)=0.0;			J(1,2)=-B1.x();
				J(2,0)=-B1.y();	J(2,1)=B1.x();		J(2,2)=0.0;
			}
			return J;
		}


		void function(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt)
		{
			//int begin=0, end=this->parameters()->size(), div=1;
			//Range splitR=MPIsplitLoop(begin,end, div);

			//iterator myIt(*(this->parameters()), splitR);
			iterator myIt(*(this->parameters()));
			while(myIt)
			{
				dMdt[i]=(realPulse(t,M[i], dMdt[i], myIt.symbol()));
				++myIt;
			}

			//if(div>1) MPIreconstruct(dMdt, begin, end);
			//MPIscatter(dMdt);

			if(CalcVar_==BlochOps::Variational) //calc the variational parts if we want
			{
				variationalFunction( t,M,dMdt);
			}
		}

//function to calculate the variational parts of the diff eq
// this goes from N...N*3 in the 'M' vector
		void variationalFunction(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt)
		{
			int i=size();
			static rmatrix J(3,3,0); //our jaciboian matrix

		//	int begin=0, end=this->parameters()->size(), div=1;
		//	Range splitR=MPIsplitLoop(begin,end, div);

		//	iterator myIt(*(this->parameters()), splitR);
			iterator myIt(*(this->parameters()));
			while(myIt)
			{
				J=pulseJacobian(t,&myIt);//+myIt.jacobian();
				dMdt[i]=J*M[i];
				dMdt[i+1]=J*M[i+1];
				dMdt[i+2]=J*M[i+2];
				i+=3;
				++myIt;
			}

			//if(div>1) MPIreconstruct(dMdt, begin, end);
		}

//all along the diagonal...
		rmatrix magneticField( Vector<coord<> > &M)
		{
			if(plist)
			{
				rmatrix tmp(pars->size()*3, pars->size()*3, 0);
				iterator myIt(*pars);
				int ct=0, i=0;
				while(myIt)
				{
					tmp(ct, ct)=B1.x();
					tmp(ct+1, ct+1)=B1.y();
					tmp(ct+2, ct+2)=B1.z();
					++i; ct+=3;
					++myIt;
				}
				return tmp;//+BlochBasic<ParamSet>::MagenticField(M);
			}
			return rmatrix(pars->size()*3, pars->size()*3, 0); //BlochBasic<ParamSet>::MagenticField(M);
		}

		rmatrix evolutionMatrix(Vector<coord<> > &M)
		{
			if(plist)
			{
				rmatrix tmp(pars->size()*3, pars->size()*3, 0);
				iterator myIt(*pars);
				int ct=0, i=0;
				while(myIt)
				{
					coord<> B1=(plist->CoordPulse(t,sym));

					tmp(ct, ct)=(M[i].z()*B1.y()-M[i].y()*B1.z());
					tmp(ct+1, ct+1)=(M[i].x()*B1.z()-M[i].z()*B1.x());
					tmp(ct+2, ct+2)=(M[i].y()*B1.x()-M[i].x()*B1.y());
					++i; ct+=3;
					++myIt;
				}
				return tmp;//+BlochBasic<ParamSet>::evolutionMatrix(M);
			}
			return rmatrix(pars->size()*3, pars->size()*3, 0);
				//BlochBasic<ParamSet>::evolutionMatrix(M);
		}
};

//This is the special class that contains NO pulses...but OTHER interactions
template<class ParamSet, class Interaction_t>
class Bloch<ParamSet, NullPulse, Interaction_t, double>: public BlochBasic<ParamSet>{
	private:
		Interaction_t *otherints_;

	public:

		typedef typename BlochBasic<ParamSet>::iterator iterator;


		Bloch():
			BlochBasic<ParamSet>() {}

		Bloch(int N):
			BlochBasic<ParamSet>(N){}

		Bloch(std::string sp):
			BlochBasic<ParamSet>(sp) {}

		Bloch(int N, std::string sp):
			BlochBasic<ParamSet>(N,sp) {}

		Bloch(int N, std::string sp, ParamSet &in):
			BlochBasic<ParamSet>(N,sp,in){}

		Bloch(int N, std::string sp, ParamSet &in, Interaction_t &intin):
			BlochBasic<ParamSet>(N,sp,in)
		{
			otherints=&intin;
		}

		Bloch(ParamSet &in):
			BlochBasic<ParamSet>(in){}

		Bloch(ParamSet &in, Interaction_t &intin):
			BlochBasic<ParamSet>(in)
		{
			otherints_=&intin;
		}

		Bloch(Bloch &in):
			BlochBasic<ParamSet>(in)
		{
			otherints_=in.otherints_;
		}

		~Bloch(){ }



		inline void setInteractions(Interaction_t &in){	otherints_=&in;	}
		inline void setInteractions(Interaction_t *in){	otherints_=in;	}
		inline Interaction_t *interactions()	const{	return otherints_;	}


		Bloch &operator=(const Bloch &rhs)
		{
			if(this==&rhs) return *this;
			BlochBasic<ParamSet>::operator=(rhs);
			otherints_=rhs.otherints;
			return *this;
		}

		template<class otherInters>
		Bloch &operator=(const Bloch<ParamSet, Pulse, otherInters, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			return *this;
		}

		Bloch &operator=(const Bloch<ParamSet, Pulse, Interaction_t, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			otherints_=rhs.Interactions();
			return *this;
		}

		void function(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt)
		{
		//	int begin=0, end=this->size(), div=1;
		//	Range splitR=MPIsplitLoop(begin,end, div);
//std::cout<<std::endl<<std::endl<<"INter: Range: "<<splitR<<" rank: "<<MPIrank<<" size: "<<this->size()<<std::endl<<std::endl;

		//	iterator myIt(*(this->parameters()), splitR);
			iterator myIt(*(this->parameters()));
			coord<> tM(0);
			if(otherints_->TotalMag)
			{
				tM=sum(M(Range(0, myIt.size()-1)));
		//		tM=sum(M(splitR));
		//		MPIreduce(tM, Reduce::Add);
			}
			if(otherints_->PreCalc) otherints_->preCalculate(t,M, dMdt, this->parameters(), tM);
			//std::cout<<"-----"<<std::endl;
			while(myIt)
			{
				dMdt[myIt.curpos()]=(otherints_->function(t, M[myIt.curpos()], dMdt[myIt.curpos()],&myIt, tM));
				//std::cout<<dMdt[myIt.curpos()]<<" | "<<M[myIt.curpos()]<<std::endl;
				++myIt;
			}
		//	if(div>1) MPIreconstruct(dMdt, begin, end);
		//	MPIscatter(dMdt);

			if(CalcVar_==BlochOps::Variational) //calc the variational parts if we want
			{
				variationalFunction( t,M,dMdt, tM);
			}
		}

//function to calculate the variational parts of the diff eq
// this goes from N...N*3 in the 'M' vector
		void variationalFunction(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt, coord<> &tM)
		{
		//	int begin=0, end=this->size(), div=1;
		//	Range splitR=MPIsplitLoop(begin,end, div);

			int i=this->parameters()->size();
			int j=0;
			static rmatrix J(3,3,0); //our jaciboian matrix
			//iterator myIt(*(this->parameters()), splitR);
			iterator myIt(*(this->parameters()));
			//dMdt.fill(ZeroType<coord<> >::zero());
			while(myIt)
			{
				J=otherints_->jacobian(t, M[j], &myIt, tM);//+myIt.jacobian();
				dMdt[i]=J*M[i];
				dMdt[i+1]=J*M[i+1];
				dMdt[i+2]=J*M[i+2];
				i+=3;
				++j;
				++myIt;
			}

			//if(div>1) MPIreconstruct(dMdt, begin, end);
			//MPIscatter(dMdt);
		}

		rmatrix magneticField( Vector<coord<> > &M)
		{
			coord<> tM(0);
			if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
			return otherints_->magneticField(0.0,M, this->parameters(), 	tM);//+ BlochBasic<ParamSet>::magneticField(M);
		}

		rmatrix evolutionMatrix( Vector<coord<> > &M)
		{
			coord<> tM(0);
			if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
			return otherints_->evolutionMatrix(0.0,M, this->parameters(), 	tM);//+ BlochBasic<ParamSet>::evolutionMatrix(M);
		}

};


//This is the special class that contains WITH pulses... AND other interactions
// Uses the class "Pulse" found in "Pulse.h"

template<class ParamSet, class Interaction_t>
class Bloch<ParamSet, Pulse,Interaction_t, double>: public BlochBasic<ParamSet>{
	private:
		Pulse *plist;
		Interaction_t *otherints_;
	public:

		typedef typename BlochBasic<ParamSet>::iterator iterator;

		Bloch():
			BlochBasic<ParamSet>(), plist(0) {}

		Bloch(int N):
			BlochBasic<ParamSet>(N), plist(0) {}

		Bloch(std::string sp):
			BlochBasic<ParamSet>(sp), plist(0) {}

		Bloch(int N, std::string sp):
			BlochBasic<ParamSet>(N,sp), plist(0) {}

		Bloch(int N, std::string sp, ParamSet &in):
			BlochBasic<ParamSet>(N,sp, in), plist(0) {}

		Bloch(int N, std::string sp, ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(N,sp,in)
		{
			plist=&pl;
		}

		Bloch(int N, ParamSet &in):
			BlochBasic<ParamSet>(N,in), plist(0) {}

		Bloch(int N, ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(N,in)
		{
			plist=&pl;
		}

		Bloch(int N, ParamSet &in, Pulse &pl,Interaction_t &ints):
			BlochBasic<ParamSet>(N,in)
		{
			plist=&pl;
			otherints_=&ints;
		}

		Bloch(ParamSet &in):
			BlochBasic<ParamSet>(in) {}

		Bloch(ParamSet &in, Pulse &pl):
			BlochBasic<ParamSet>(in)
		{
			plist=&pl;
		}

		Bloch(ParamSet &in, Pulse &pl,Interaction_t &ints ):
			BlochBasic<ParamSet>(in)
		{
			plist=&pl;
			otherints_=&ints;
		}

		Bloch(Bloch &in):
			BlochBasic<ParamSet>(in), plist(in.plist),otherints_(in.otherints_)
		{}


		Bloch(Bloch<ParamSet, NullPulse, Interaction_t, double> &in):
			BlochBasic<ParamSet>(in), plist(0),otherints_(in.otherints_)
		{}

		~Bloch()
		{
			plist=0;
		}

		inline void setInteractions(Interaction_t &in){	otherints_=&in;	}
		inline void setInteractions(Interaction_t *in){	otherints_=in;	}

		inline void setPulses(Pulse &in){	plist=&in;	}
		inline void setPulses(Pulse *in){	plist=in;	}


		inline Pulse *pulses()	const	 {	return plist;	}
		inline Interaction_t *interactions()	const{	return otherints_;	}
		inline Pulse *Pulses()	const	 {	return plist;	}
		inline Interaction_t *Interactions()	const{	return otherints_;	}

		void operator=(const Bloch &rhs)
		{
			if(this==&rhs) return;
			BlochBasic<ParamSet>::operator=(rhs);
			otherints_=rhs.otherints_;
			plist=rhs.plist;
		}

		template<class otherInters>
		void operator=(const Bloch<ParamSet, NullPulse,otherInters, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			plist=0;
		}

		void operator=(const Bloch<ParamSet, NullPulse, Interaction_t, double> &rhs)
		{
			BlochBasic<ParamSet>::operator=(rhs);
			plist=0;
			otherints_=rhs.Interactions();
		}
/*
		inline void realPulse(double t, coord<> &iny, coord<> &dydt, const std::string &sym)
		{
			if(plist)
			{
				coord<> B1=(plist->CoordPulse(t,sym));
				dydt.x()-=(B1.z()*iny.y()-B1.y()*iny.z());
				dydt.x()-=(B1.x()*iny.z()-B1.z()*iny.x());
				dydt.x()-=(B1.y()*iny.x()-B1.x()*iny.y());
			}
		}
*/
		inline coord<> realPulse(double t, coord<> &iny, coord<> &dydt, const std::string &sym)
		{
			if(plist)
			{
				//coord<> B1=(plist->CoordPulse(t,sym));
				return cross(iny, plist->CoordPulse(t,sym));
				/*coord<>(
						-(B1.z()*iny.y()-B1.y()*iny.z()),
						-(B1.x()*iny.z()-B1.z()*iny.x()),
						-(B1.y()*iny.x()-B1.x()*iny.y())
					);*/
			}
			return ZeroType<coord<> >::zero();
		}

		template<class Piter>
		inline rmatrix pulseJacobian(double t,Piter *onit)
		{
			static rmatrix J(3,3,0);
			static coord<> B1;
			if(plist)
			{
				B1=(plist->CoordPulse(t,onit->symbol()));
				J(0,0)=0.0; 	J(0,1)=-B1.z(); 	J(0,2)=B1.y();
				J(1,0)=B1.z();	J(1,1)=0.0;			J(1,2)=-B1.x();
				J(2,0)=-B1.y();	J(2,1)=B1.x();		J(2,2)=0.0;
			}
			return J;
		}


		void function(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt){
			//parameters()->reset();
			//int begin=0, end=this->size(), div=0;
			//Range splitR=MPIsplitLoop(begin,end, div);
		//std::cout<<std::endl<<std::endl<<"PULSE: Range: "<<splitR<<" rank: "<<MPIrank<<" size: "<<this->size()<<std::endl<<std::endl;
			//int bell=0;
			iterator myIt(*(this->parameters()));
			//iterator myIt(*(this->parameters()), splitR);
			coord<> tM(0);
			if(otherints_->TotalMag)
			{
				tM=sum(M(Range(0, myIt.size()-1)));
			//	tM=sum(M(splitR));
			//	MPIreduceAndScatter(tM, Reduce::Add);
			}
			//cout<<"MM: "<<M<<endl;
			if(otherints_->PreCalc) otherints_->preCalculate(t,M, dMdt, this->parameters(), tM);
			//int ct=0;
			while(myIt)
			{
				dMdt[myIt.curpos()]=
				   realPulse(t,M[myIt.curpos()], dMdt[myIt.curpos()], myIt.symbol())+
				   otherints_->function(t, M[myIt.curpos()], dMdt[myIt.curpos()], &myIt, tM);

				//cout<<"|P "<<realPulse(t,M[myIt.curpos()], dMdt[myIt.curpos()], myIt.symbol())<<" |O "
				//<<		otherints_->function(t, M[myIt.curpos()], dMdt[myIt.curpos()], &myIt, tM)
				//<<" |M "<<M[myIt.curpos()]
				//		<<endl;
				++myIt;

			}
			//cout<<"MM2: "<<M<<endl;

			//std::cout<<"rank: "<<MPIrank<<" PRE: "<<dMdt<<std::endl;
			//if(div>1) MPIreconstruct(dMdt, begin, end);
			//MPIscatter(dMdt);
			//std::cout<<"rank: "<<MPIrank<<"POST: "<<dMdt<<std::endl;

			//MPI_Barrier(MPI_COMM_WORLD);

			if(CalcVar_==BlochOps::Variational) //calc the variational parts if we want
			{
				variationalFunction( t,M,dMdt, tM);
			}
		}
/*
		void function(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt){
			//parameters()->reset();
			iterator myIt(*(this->parameters()));
			static coord<> tM(0);
			tM(0,0,0);
			if(otherints_->totalMag) tM=sum(M);
			if(otherints_->PreCalc) otherints_->preCalculate(t,M, dMdt, this->parameters(), tM);
			dMdt.fill(ZeroType<coord<> >::zero());
			while(myIt)
			{
				realPulse(t,M[myIt.curpos()], dMdt[myIt.curpos()], myIt.symbol());
				otherints_->function(t, M[myIt.curpos()], dMdt[myIt.curpos()], &myIt, tM);
				++myIt;
			}

			if(CalcVar_==BlochOps::Variational) //calc the variational parts if we want
			{
				variationalFunction( t,M,dMdt, tM);
			}
		}
*/
//function to calculate the variational parts of the diff eq
// this goes from N...N*3 in the 'M' vector
		void variationalFunction(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt, coord<> &tM)
		{
			//int begin=0, end=this->size(), div=1;
			//Range splitR=MPIsplitLoop(begin,end, div);

			int i=size();
			int j=0;
			rmatrix J(3,3,0); //our jaciboian matrix
			iterator myIt(*(this->parameters()));
			while(myIt)
			{
				J=pulseJacobian(t,&myIt)+otherints_->jacobian(t, M[j], &myIt, tM);//+myIt.jacobian()
				dMdt[i]=J*M[i];
				dMdt[i+1]=J*M[i+1];
				dMdt[i+2]=J*M[i+2];
				i+=3;
				++j;
				++myIt;
			}
			//if(div>1) MPIreconstruct(dMdt, begin, end);
			//MPIscatter(dMdt);
		}
//all along the diagonal...
		rmatrix magneticField( Vector<coord<> > &M)
		{
			coord<> tM(0);
			tM(0,0,0);
			if(plist)
			{
				rmatrix tmp(pars->size()*3, pars->size()*3, 0);
				iterator myIt(*pars);
				int ct=0, i=0;
				while(myIt)
				{
					tmp(ct, ct)=B1.x();
					tmp(ct+1, ct+1)=B1.y();
					tmp(ct+2, ct+2)=B1.z();
					++i; ct+=3;
					++myIt;
				}
				if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
			return tmp+
						//BlochBasic<ParamSet>::MagenticField(M)+
						otherints_->MagenticField(0, M, this->parameters(), tM);
			}
			if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
			return //BlochBasic<ParamSet>::MagenticField(M)+
						otherints_->MagenticField(0, M, this->parameters(), tM);
		}

		rmatrix evolutionMatrix(Vector<coord<> > &M)
		{
			coord<> tM(0);
			if(plist)
			{
				rmatrix tmp(pars->size()*3, pars->size()*3, 0);
				iterator myIt(*pars);
				int ct=0, i=0;
				while(myIt)
				{
					coord<> B1=(plist->CoordPulse(t,sym));

					tmp(ct, ct)=(M[i].z()*B1.y()-M[i].y()*B1.z());
					tmp(ct+1, ct+1)=(M[i].x()*B1.z()-M[i].z()*B1.x());
					tmp(ct+2, ct+2)=(M[i].y()*B1.x()-M[i].x()*B1.y());
					++i; ct+=3;
					++myIt;
				}
				if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
				return tmp+
						//BlochBasic<ParamSet>::evolutionMatrix(M)+
						otherints_->evolutionMatrix(0.0, M, this->parameters(), tM);
			}
			if(otherints_->TotalMag) tM=sum(M(Range(Range::Start, size()-1)));
			return //BlochBasic<ParamSet>::evolutionMatrix(M)+
						otherints_->evolutionMatrix(0.0, M, this->parameters(), tM);
		}
};


END_BL_NAMESPACE

#endif
