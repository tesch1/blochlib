

/* bloch_basic.h ********/


 /*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 07-30-01
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
 	bloch_basic.h-->the basic Bloch interactions offsets, T2, T1...serves as the base class for
 	more complex bloch interactions...
 */


#ifndef _bloch_basic_h_
#define _bloch_basic_h_ 1

#include "utils/constants.h"
#include "container/grids/coords.h"

#include "bloch/listblochpars.h"
#include "bloch/pulse.h"
#include <string>

BEGIN_BL_NAMESPACE


//when calculating the Spin tragectories, this flag will also
// alow the calculation of the Variation Equations
// THe variational equations are a set of matrix equations
// that calulate the 'deviation' from the current tragectory path
// dH(x,t)/dt= jacobian(x,t)H(x,t)
// Each interaction will have a jacobian assosicated with it,
// and so for each spin there will be 9 equations (each spin having a (Mx,My,Mz))

//this enum tells the system to calculate them as the tragectory is being calulated
// and thus sets up the nessesary extra data strucutres
// the default will be NOT to calculate the variationals
class BlochOps{
	public:
		enum VariationEqs{Variational, NoVariational};
};


template<class ParamSet, class T=double>
class BlochBasic;

//This is the special class that contains NO pulses...and NO Interactions
// for the 'Double' offset represnetation and the double Mo and Bo rep
template<class GridEngine_t, int BPops, class T>
class BlochBasic<ListBlochParams<GridEngine_t, BPops, double>, T>
{
	private:
		ListBlochParams<GridEngine_t, BPops, double > *pars;			//NOTE::this is a pointer...avoids copying the mass data
		Vector<coord<T> > Mag;  //the first 'N' entries will contain the spin equations, the N..N*3
		Vector<coord<T> > Mag0;
		int size_;
		int tagged_;

		bool todelPars;

		int cursize_;	//will be 'size' or 'size*3' depending on the

	protected:
		BlochOps::VariationEqs CalcVar_; //make this visibile by the sub classes

	public:


		typedef typename ListBlochParams<GridEngine_t, BPops, double >::iterator iterator;

		BlochBasic():
		   size_(1), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=NULL;
			todelPars=false;
		}

		BlochBasic(int N):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, double>(N,"1H");
			todelPars=true;
			setInit();
		}

		BlochBasic(std::string sp):
		  Mag(1,0), size_(1), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, double>(1,sp);
			todelPars=true;
			setInit();
		}

		BlochBasic(int N, std::string sp):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, double>(N,sp);
			todelPars=true;
			setInit();
		}

		BlochBasic(int N, std::string sp, ListBlochParams<GridEngine_t, BPops, double > &in):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=&in;
			todelPars=false;
			setInit();
		}

		BlochBasic(ListBlochParams<GridEngine_t, BPops, double > &in):
		  Mag(in.size(),0), size_(in.size()), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=&in;
			todelPars=false;
			setInit();
		}

		BlochBasic(BlochBasic &in):
		  pars(in.pars), Mag(in.size(),0), size_(in.size()), tagged_(0),
			todelPars(in.todelPars), CalcVar_(BlochOps::NoVariational)
		{
			setInit();
		}

		~BlochBasic()
		{
			if(todelPars){ delete pars;	}
			else{	pars=0;}
		}



		inline void calcVariational()
		{
			CalcVar_=BlochOps::Variational;
			cursize_=size_*3;
			Mag.resizeAndPreserve(size_*4); //(1 size for spins, 3 sizes for variational eqs)
			for(int i=size_;i<size_*4;i+=3)
			{
				Mag[i].x(1); Mag[i].y(0); Mag[i].z(0);
				Mag[i+1].x(0); Mag[i+1].y(1); Mag[i+1].z(0);
				Mag[i+2].x(0); Mag[i+2].y(0); Mag[i+2].z(1);
			}
			Mag0.resizeAndPreserve(cursize_);
		}

		inline BlochOps::VariationEqs VariationalPolicy() const { 	return CalcVar_;	}
		inline BlochOps::VariationEqs variationalPolicy() const { 	return CalcVar_;	}

		inline void noCalcVariational()
		{
			CalcVar_=BlochOps::NoVariational;
			cursize_=size_;
			Mag.resizeAndPreserve(size_, 0);
		}

		Vector<coord<> > curVariational(){
			if(Mag.size()!=4*pars->size()){
				std::cerr<<std::endl<<"Warning: BlochBasic::curVariational()"<<std::endl;
				std::cerr<<" No BlochOps::Variational's initialized...returning null vector"<<std::endl;
				return Vector<coord<> >();
			}else{
				return Mag(Range(size_, 4*size_-1));
			}
		}


		inline ListBlochParams<GridEngine_t, BPops, double > *parameters(){	return pars;	}

		//template<class GridEngine_t, class BPops>
		inline void setParameters(ListBlochParams<GridEngine_t, BPops, double > *in)
		{
			if(todelPars){
				delete pars;
				todelPars=false;
			}
			pars=&in; setInit();
		}

		//template<class GridEngine_t, class BPops>
		inline void setParameters(ListBlochParams<GridEngine_t, BPops, double > &in)
		{
			if(todelPars){
				delete pars;
				todelPars=false;
			}
			pars=&in; setInit();
		}

		inline int size() const{ return size_;	}
		inline int tagged() const { return tagged_;	}
		void tag(int in){
			if(in<size() && in>=0){
				tagged_=in;
			}else{
				std::cerr<<std::endl<<"Error: Bloch.tag()"<<std::endl;
				std::cerr<<" You wish to tag the spin "<<in<<std::endl;
				std::cerr<<" But there are only "<<size()<<" spins in the system"<<std::endl;
				std::cerr<<" setting the tag to the first spin "<<std::endl;
				tagged_=0;
			}
		}

		void setInit(){
			if(pars!=NULL){
				size_=pars->size();
				if(CalcVar_==BlochOps::Variational){
					Mag.resize(size_*4);
				}else{
					Mag.resize(size_);
				}
				if(size_>1)
				{
					iterator myIt(*pars);
					int i=0;
					/*while(myIt && size_){
						Mag[i].x(0);
						Mag[i].y(0);
						Mag[i].z(myIt.Mo());
						++i; ++myIt;
					}*/
					switch(pars->initialCondition()){
						case InitialCondition::HalfUpHalfDown:
							while(myIt && size_){
								Random<UniformRandom<> > angleR(0.0, 1.0);
								myIt.setInitialCondition(pars->initialCondition(), (angleR()>0.5)?0:1);
								Mag[i]=coord<>(0.0,0.0,myIt.Mo());
								++i; ++myIt;
							}
							break;
						case InitialCondition::AllUp: //should be set already
							while(myIt && size_){
								Mag[i]=coord<>(0.0,0.0,myIt.Mo());
								++i; ++myIt;
							}
							break;
						default:
							while(myIt && size_){
								myIt.setInitialCondition(pars->initialCondition());
								Mag[i]=coord<>(0.0,0.0,myIt.Mo());
								++i; ++myIt;
							}
							break;
					}
				//the initial condition for the variation equations is the identity matrix
					if(CalcVar_==BlochOps::Variational)
					{
						for(int i=size_;i<size_*4;i+=3)
						{
							Mag[i].x(1); Mag[i].y(0); Mag[i].z(0);
							Mag[i+1].x(0); Mag[i+1].y(1); Mag[i+1].z(0);
							Mag[i+2].x(0); Mag[i+2].y(0); Mag[i+2].z(1);
						}
					}
				}else{
					Mag[0]=coord<>(0.0,0.0,(*pars)[0].BasicBlochParams<double>::Mo());
				}
				Mag0=Mag;
			}
		}

		void operator=(const BlochBasic &rhs)
		{
			if(this==&rhs) return;
			pars=rhs.pars;
			Mag=rhs.Mag;
			Mag0=rhs.Mag0;
			size_=rhs.size_;
			todelPars=rhs.todelPars;
			tagged_=rhs.tagged_;
			CalcVar_=rhs.CalcVar_;
		}

		typename ListBlochParams<GridEngine_t, BPops, double >::Parameter_T &operator()(int i) { return (*pars)[i];	}

		Vector<coord<T> > initialCondition(){ return Mag0;	}
		Vector<coord<T> > currentMag(){ return Mag;	}

		void setInitialCondition(Vector<coord<> > const &in){ Mag0=in;	Mag=in;}
		void setInitialCondition(coord<>  const &in, int on){ Mag0(on)=in; Mag0(on)=in;	}
		void setCurrentMag(Vector<coord<> > const &in){ Mag=in;	}
		void setCurrentMag(Vector<coord<> > const &in, int on){ Mag(on)=in;	}



		void deltaPulse(double theta, double phi, double psi=0, int on=-1){
			if(on==-1){
				for(int i=0;i<size();i++){
					Mag(i).Rotate3D(theta, phi, psi);
				}
			}else{
				if(on<size())	Mag(on).Rotate3D(theta, phi, psi);
			}
		}

		void deltaPulse(std::string on, double theta, double phi, double psi=0){
			if(on==""){
				for(int i=0;i<size();i++){
					Mag(i).Rotate3D(theta, phi, psi);
				}
			}else{
				for(int i=0;i<size();i++){
					if((*pars)(i).symbol()==on) Mag(i).Rotate3D(theta, phi, psi);
				}
			}
		}

		coord<T> totalM(){	return sum(Mag);		}

		void printData(std::ostream &oo){
			oo<<Mag[tagged_]<<std::endl;
		}

		void printPars(std::ostream &oo){
			oo<<(*pars)<<std::endl;
		}

		bool savePars(std::fstream &out){
			return pars->write(out);
		}

		bool savePars(std::string &out){
			std::fstream oo(out.c_str(), std::ios::out|std::ios::binary);
			return savePars(oo);
		}
		bool readPars(std::fstream &in){
			return pars->read(in);
		}

		bool readPars(std::string &in){
			std::fstream ii(in.c_str(), std::ios::in | std::ios::binary);
			return readPars(ii);
		}

};


//This is the special class that contains NO pulses...and NO Interactions
// for the 'COORD' offset represnetation and the COORD Mo and Bo rep

template<class GridEngine_t, int BPops, class T>
class BlochBasic<ListBlochParams<GridEngine_t, BPops, coord<> >, T>
{
	private:
		ListBlochParams<GridEngine_t, BPops, coord<> > *pars;			//NOTE::this is a pointer...avoids copying the mass data
		Vector<coord<T> > Mag;  //the first 'N' entries will contain the spin equations, the N..N*3
		Vector<coord<T> > Mag0;
		int size_;
		int tagged_;
		bool todelPars;
		int cursize_;	//will be 'size' or 'size*3' depending on the

	protected:
		BlochOps::VariationEqs CalcVar_; //make this visibile by the sub classes

	public:


		typedef typename ListBlochParams<GridEngine_t, BPops, coord<> >::iterator iterator;

		BlochBasic():
		 size_(0), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=NULL;
			todelPars=false;
		}

		BlochBasic(int N):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, coord<> >(N,"1H");
			todelPars=true;
			setInit();
		}

		BlochBasic(std::string sp):
		  Mag(1,0), size_(1), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, coord<> >(1,sp);
			todelPars=true;
			setInit();
		}

		BlochBasic(int N, std::string sp):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=new ListBlochParams<GridEngine_t, BPops, coord<> >(N,sp);
			todelPars=true;
			setInit();
		}

		BlochBasic(int N, std::string sp, ListBlochParams<GridEngine_t, BPops, coord<> >  &in):
		  Mag(N,0), size_(N), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=&in;
			todelPars=false;
			setInit();
		}

		BlochBasic(ListBlochParams<GridEngine_t, BPops, coord<> > &in):
		  Mag(in.size(),0), size_(in.size()), tagged_(0),CalcVar_(BlochOps::NoVariational)
		{
			pars=&in;
			todelPars=false;
			setInit();
		}

		BlochBasic(BlochBasic &in):
		  pars(in.pars), Mag(in.size(),0), size_(in.size()), tagged_(0),
		  todelPars(in.todelPars),	CalcVar_(BlochOps::NoVariational)
		{
			setInit();
		}

		~BlochBasic()
		{
			if(todelPars){	delete pars;	}
			else{ pars=0;	}
		}

		inline void calcVariational()
		{
			CalcVar_=BlochOps::Variational;
			cursize_=size_*4;
			Mag.resizeAndPreserve(size_*4); //(1 size for spins, 3 sizes for variational eqs)
			for(int i=size_;i<size_*4;i+=3)
			{
				Mag[i].x(1); Mag[i].y(0); Mag[i].z(0);
				Mag[i+1].x(0); Mag[i+1].y(1); Mag[i+1].z(0);
				Mag[i+2].x(0); Mag[i+2].y(0); Mag[i+2].z(1);
			}
			Mag0.resizeAndPreserve(cursize_);

		}

		inline BlochOps::VariationEqs VariationalPolicy()const { 	return CalcVar_;	}
		inline BlochOps::VariationEqs variationalPolicy()const { 	return CalcVar_;	}

		inline void  noCalcVariational()
		{
			CalcVar_=BlochOps::NoVariational;
			cursize_=size_;
			Mag.resizeAndPreserve(size_, 0);
		}

		Vector<coord<> > curVariational(){
			if(Mag.size()!=4*pars->size()){
				std::cerr<<std::endl<<"Warning: BlochBasic::curVariational()"<<std::endl;
				std::cerr<<" No BlochOps::Variational's initialized...returning null vector"<<std::endl;
				return Vector<coord<> >();
			}else{
				return Mag(Range(size_, 4*size_-1));
			}
		}


		inline ListBlochParams<GridEngine_t, BPops, coord<> > *parameters(){	return pars;	}
		inline void setParameters(ListBlochParams<GridEngine_t, BPops, coord<> > &in)
		{
			if(todelPars){
				delete pars;
				todelPars=false;
			}
			pars=&in; setInit();
		}

		inline void setParameters(ListBlochParams<GridEngine_t, BPops, coord<> > *in)
		{
			if(todelPars){
				delete pars;
				todelPars=false;
			}
			pars=&in; setInit();
		}

		inline int size() const{ return size_;	}
		inline int tagged() const { return tagged_;	}
		void tag(int in){
			if(in<size() && in>=0){
				tagged_=in;
			}else{
				std::cerr<<std::endl<<"Error: Bloch.tag()"<<std::endl;
				std::cerr<<" You wish to tag the spin "<<in<<std::endl;
				std::cerr<<" But there are only "<<size()<<" spins in the system"<<std::endl;
				std::cerr<<" setting the tag to the first spin "<<std::endl;
				tagged_=0;
			}
		}

		void setInit(){
			if(pars!=NULL){
				size_=pars->size();
				if(CalcVar_==BlochOps::Variational){
					Mag.resize(size_*4);
				}else{
					Mag.resize(size_);
				}
				if(size_>1)
				{
					iterator myIt(*pars);
					int i=0;
					switch(pars->initialCondition()){
						case InitialCondition::HalfUpHalfDown:
							while(myIt && size_){
								Random<UniformRandom<> > angleR(0.0, 1.0);
								myIt.setInitialCondition(pars->initialCondition(), (angleR()>0.5)?0:1);
								Mag[i]=myIt.Mo();
								++i; ++myIt;
							}
							break;
						case InitialCondition::AllUp: //should be set already
							while(myIt && size_){
								Mag[i]=myIt.Mo();
								++i; ++myIt;
							}
							break;
						default:

							while(myIt && size_){
								myIt.setInitialCondition(pars->initialCondition());
								Mag[i]=myIt.Mo();
								++i; ++myIt;
							}
							break;
					}


			//the initial condition for the variation equations is the identity matrix
					if(CalcVar_==BlochOps::Variational)
					{
						for(int i=size_;i<size_*4;i+=3)
						{
							Mag[i].x(1); Mag[i].y(0); Mag[i].z(0);
							Mag[i+1].x(0); Mag[i+1].y(1); Mag[i+1].z(0);
							Mag[i+2].x(0); Mag[i+2].y(0); Mag[i+2].z(1);
						}
						Mag0.resize(cursize_);
					}
				}else{
					Mag[0]=(*pars)[0].BasicBlochParams<coord<> >::Mo();
				}
				Mag0=Mag;
			}
		}

		void operator=(const BlochBasic &rhs)
		{
			if(this==&rhs) return;
			pars=rhs.pars;
			Mag=rhs.Mag;
			Mag0=rhs.Mag0;
			size_=rhs.size_;
			todelPars=rhs.todelPars;
			tagged_=rhs.tagged_;
			CalcVar_=rhs.CalcVar_;
		}

		typename ListBlochParams<GridEngine_t, BPops, coord<> >::Parameter_T &operator()(int i) { return (*pars)[i];	}

		Vector<coord<T> > initialCondition(){ return Mag0;	}
		Vector<coord<T> > currentMag(){ return Mag;	}

		void setInitialCondition(Vector<coord<> > const &in){ Mag0=in;	}
		void setInitialCondition(coord<>  const &in, int on){ Mag0(on)=in;	}


		void setCurrentMag(Vector<coord<> > const &in){ Mag=in;	}
		void setCurrentMag(Vector<coord<> > const &in, int on){ Mag(on)=in;	}


		void deltaPulse(double theta, double phi, double psi=0, int on=-1){
			if(on==-1){
				for(int i=0;i<size();i++){
					Mag(i).Rotate3D(theta, phi, psi);
				}
			}else{
				if(on<size())	Mag(on).Rotate3D(theta, phi, psi);
			}
		}

		void deltaPulse(double theta, double phi, double psi=0, std::string on=""){
			if(on==""){
				for(int i=0;i<size();i++){
					Mag(i).Rotate3D(theta, phi, psi);
				}
			}else{
				for(int i=0;i<size();i++){
					if((*pars)(i).symbol()==on) Mag(i).Rotate3D(theta, phi, psi);
				}
			}
		}

		coord<T> totalM(){	return sum(Mag);		}

/*
		inline void T2bloch(coord<> &iny, coord<> &dydt,const double &T2){
			dydt.x()-=iny.x()/T2;
			dydt.y()-=iny.y()/T2;
		}

		inline void T1bloch(coord<> &iny, coord<> &dydt,const double &T1,const double &Mo){
			dydt.z()+=(Mo-iny.z())/T1;
		}

//the 'master' funtion for the most basic of bloch interactions
		void function(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt){
			int i=0;
			iterator myIt(*pars);
			while(myIt)
			{
				dMdt[i].x(myIt->offset()*M[i].y());
				dMdt[i].y(-myIt->offset()*M[i].x());
				dMdt[i].z(0);
				if(myIt->T2()) T2bloch(M[i], dMdt[i], myIt->T2());
				if(myIt->T1()) T1bloch(M[i], dMdt[i], myIt->T1(), myIt->Mo.z());
				++myIt;
				++i;
			}
			if(CalcVar_==BlochOps::Variational)
			{
				variationalFunction(t, M, dMdt);
			}
		}

//the 'master' variation function for the most basic of bloch interactions
// we get 'lucky' here becuase we have no interactions or nonlinear
//pieces here, the jacobian is NOT dependant on the current posistion
		void variationalFunction(double t, Vector<coord<> > &M, Vector<coord<> > &dMdt)
		{
			int i=size_;
			static rmatrix J(3,3,0); //our jaciboian matrix
			iterator myIt(*pars);
			while(myIt)
			{
				J=myIt.jacobian();
				dMdt[i]=J*M[i];
				dMdt[i+1]=J*M[i+1];
				dMdt[i+2]=J*M[i+2];
				i+=3;
				++myIt;
			}
		}

//all along the diagonal...
		rmatrix magneticField( Vector<coord<> > &M)
		{
			rmatrix tmp(pars->size()*3, pars->size()*3, 0);
			iterator myIt(*pars);
			int ct=0;
			while(myIt)
			{
				tmp(ct, ct)=myIt.offset().x()/myIt.gamma();
				tmp(ct+1, ct+1)=myIt.offset().y()/myIt.gamma();
				tmp(ct+2, ct+2)=myIt.offset().z()/myIt.gamma();
				ct+=3;
				++myIt;
			}
			return tmp;
		}

		rmatrix evolutionMatrix(Vector<coord<> > &M)
		{
			rmatrix tmp(pars->size()*3, pars->size()*3, 0);
			iterator myIt(*pars);
			int ct=0, i=0;
			while(myIt)
			{
				//xhat direction
				tmp(ct, ct)=myIt.offset().z()*M[i].y()-myIt.offset().y()*M[i].z();

				//yhat dir
				tmp(ct+1, ct+1)=myIt.offset().x()*M[i].z()-myIt.offset().z()*M[i].x();

				//zhat
				tmp(ct+2, ct+2)=myIt.offset().y()*M[i].x()-myIt.offset().x()*M[i].y();

				if(myIt.T2())
				{
					tmp(ct,ct)-=M[i].x()/myIt.T2();
					tmp(ct+1,ct+1)-=M[i].y()/myIt.T2();
				}
				tmp(ct+2, ct+2)=myIt.offset();
				if(myIt.T1())
				{
					tmp(ct+2,ct+2)+=(myIt.Mo().z()-M[i].z())/myIt.T1();
				}
				++i; ct+=3;
				++myIt;
			}
			return tmp;
		}
*/

		void printData(std::ostream &oo){
			oo<<Mag[tagged_]<<std::endl;
		}

		void printPars(std::ostream &oo){
			oo<<(*pars)<<std::endl;
		}

		bool savePars(std::fstream &out){
			return pars->write(out);
		}

		bool savePars(std::string &out){
			std::fstream oo(out.c_str(), std::ios::out|std::ios::binary);
			return savePars(oo);
		}
		bool readPars(std::fstream &in){
			return pars->read(in);
		}

		bool readPars(std::string &in){
			std::fstream ii(in.c_str(), std::ios::in | std::ios::binary);
			return readPars(ii);
		}

};


END_BL_NAMESPACE




#endif

