

/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-27-01
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
	oneFID.h--> creates an fid given a slew of info

	Func_T is a class the Must contain these Hamiltonian function

	'H(wr, beta, theta, phi, tbegin, tend)' -->
	'H(t1, t2)' --> hamiltonian AFTER angles have been set
	'setRotorAngles(alpha, beta)' --> set the the current rotor angles
	'setPowderAngles(theta, phi)' --> set the powder angles
	'Fe()' --> identity matrix
	'
*/



#ifndef _oneFID_h_
#define _oneFID_h_ 1


#include "QMspins/driver/compute.h"
#include "QMspins/driver/powder.h"

BEGIN_BL_NAMESPACE


//MatType is either 'hmatrix'
// or dmatrix...

template<class Func_T, class MatrixType_T=hmatrix>
class oneFID
	//: public Func_T
{
	friend class compute<oneFID<Func_T,MatrixType_T>,MatrixType_T >;
	public:

		typedef  Complex<typename MatrixType_T::numtype::numtype> ComplexType_T;
		typedef _matrix<ComplexType_T, FullMatrix> PropMatType_T ;

	private:

		Vector<ComplexType_T> fid_;
		double wr_;
		double rotb_;

		double sw_;

		double startT_;

		const void NoMf()
		{
			std::cerr<<std::endl<<"Error: oneFID::FID(...)"<<std::endl;
			std::cerr<<" NO Function defined...(i.e. a class with "<<std::endl;
			std::cerr<<" hmatrix=Hamiltonian(wr, beta, theta, phi, t1, t2) "<<std::endl;
		}

		const void NptsErr()
		{
			std::cerr<<std::endl<<"Error: oneFID::FID(...)"<<std::endl;
			std::cerr<<" Fid size was '0' please use 'npts(int)' set set it"<<std::endl;
			std::cerr<<" cannot calculate the FID "<<std::endl;
		}

		Func_T *function_;

	public:


		//in order to have this work in
		// parallel, you MUST set this to a valid
		// MPI controller (typically MPIworld)
		MPIcontroller Controller;

		oneFID():
			//Func_T(),
			fid_(0,0), wr_(0), rotb_(0), sw_(20000),startT_(0),
			function_(function()),
			Controller()
		{}

		oneFID(Func_T &in, int npts, double sw=20000, double wr=0, double rotb=0, double startT=0):
			//Func_T(in),
			fid_(npts, 0), wr_(wr), rotb_(rotb), sw_(sw),startT_(startT),
			function_(&in),
			Controller()
		{}

		oneFID(int npts, double sw=20000, double wr=0, double rotb=0, double startT=0):
			fid_(npts, 0), wr_(wr), rotb_(rotb), sw_(sw),startT_(startT),
			function_(NULL),
			Controller()
		{}

		oneFID(const oneFID &rhs):
			//Func_T(rhs),
			fid_(rhs.fid_), wr_(rhs.wr_), rotb_(rhs.rotb_), sw_(rhs.sw_),startT_(rhs.startT_),
			function_(rhs.function_),
			Controller()
		{}

		~oneFID(){function_=NULL;	}

		void operator=(const Func_T &insys)
		{
			function_->operator=(insys);
		}

		void setFunction(Func_T &in)
		{	function_= (&in);	}

		void setFunction(Func_T *in)
		{	function_= (in);	}

		Func_T *function(){	return function_;	}

		inline double &wr(){	return wr_;	}
		inline double &rotorSpeed(){	return wr_;	}
		inline double wr()const{	return wr_;	}
		inline double rotorSpeed()const{	return wr_;	}

		inline void setWr(double in){	 wr_=in;	 }
		inline void setRotorSpeed(double in){	 wr_=in; 	}

		inline double &rotorAngle(){	return rotb_;	}
		inline double &rotorBeta(){	return rotb_;	}
		inline double rotorAngle()const{	return rotb_;	}
		inline double rotorBeta()const{	return rotb_;	}

		inline void setRotorAngle(double in){	 rotb_=in;	}
		inline void setRotorBeta(double in){	 rotb_=in;	}

		inline double &sw(){	return sw_;	}
		inline double &sweep(){	return sw_;	}
		inline double &SweepWidth(){	return sw_;	}
		inline double &sweepWidth(){	return sw_;	}

		inline double sw()const{	return sw_;	}
		inline double sweep()const{	return sw_;	}
		inline double SweepWidth()const{	return sw_;	}
		inline double sweepWidth()const{	return sw_;	}

		inline void setSw(double in){	 sw_=in;		}
		inline void setSweep(double in){	 sw_=in;		}
		inline void setSweepWidth(double in){	 sw_=in;		}

		inline int npts() const {	return fid_.size();	}
		inline int size() const {	return fid_.size();	}

		void npts(int in) {	return fid_.resize(in);	}
		void setNpts(int in) {	return fid_.resize(in);	}
		void setSize(int in) {	return fid_.resize(in);	}

		double startTime()const{	return startT_;	}
		double &startTime(){	return startT_;	}
		void setStartTime(double in) {	startT_=in;	}

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> FID(powder &in, matrix_T &ro, matrix_T2 &det);

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> FID(matrix_T &ro, matrix_T2 &det, double the=0, double phi=0, double gam=0);

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> staticFID(powder &in, matrix_T &ro, matrix_T2 &det);

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> singleStaticFID(matrix_T &ro, matrix_T2 &det, double theta=0, double phi=0, double gamma=0);

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> spinningFID(powder &in, matrix_T &ro, matrix_T2 &det);

		template<class matrix_T, class matrix_T2>
		Vector<ComplexType_T> singleSpinningFID(matrix_T &ro, matrix_T2 &det, double theta=0, double phi=0, double gamma=0);

		MatrixType_T   Hamiltonian( double t1, double t2, double wrIn=0.0);
};

template<class Func_T, class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
 oneFID<Func_T, MatrixType_T>::FID( powder &in, matrix_T &ro, matrix_T2 &det)
{
	if(fid_.size()<=0){	NptsErr();	return fid_;	}
	if(wr_==0.0 || rotb_==0.0){	return staticFID(in, ro, det);	}
	return spinningFID(in,ro,det);
}

template<class Func_T, class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
  oneFID<Func_T, MatrixType_T>::FID( matrix_T &ro, matrix_T2 &det, double the, double phi, double gam)
{
	if(fid_.size()<=0){	NptsErr();	return fid_;	}
	if(wr_==0.0 || rotb_==0.0){	return singleStaticFID(ro, det,the,phi);	}
	return singleSpinningFID(ro,det, the, phi, gam);
}

template<class Func_T, class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
  oneFID<Func_T, MatrixType_T>::staticFID(powder &in, matrix_T &ro, matrix_T2 &det)
{
	powder::iterator myit(in);
	fid_.fill(0);
	double cTh, cPh, cGa, cW;
	int PowTag=30, done=-1, amdone=1;
//first make sure we are Parrellel
	if(Controller.parallel()  && in.size()>1){
	//Perform the "master/slave" loops
		if(Controller.master()){
		//First we see if we send an initial request to each
			for(int rank=1;rank<Controller.size();++rank){
				Controller.put(amdone, rank, PowTag);
				Controller.put(myit.theta(), rank, PowTag+1);
				Controller.put(myit.phi(), rank, PowTag+2);
				Controller.put(myit.gamma(), rank, PowTag+3);
				Controller.put(myit.weight(), rank, PowTag+4);
				++myit;
				if(!myit) break;
			}
			int get, dd;
			while(myit){
				get=Controller.getAny(dd);
				Controller.put(amdone, get, PowTag);
				Controller.put(myit.theta(), get, PowTag+1);
				Controller.put(myit.phi(), get, PowTag+2);
				Controller.put(myit.gamma(), get, PowTag+3);
				Controller.put(myit.weight(), get, PowTag+4);
				++myit;
			}

			for(int rank=1;rank<Controller.size();++rank)
				Controller.put(done, rank, PowTag);

		//Slave
		}else{
			while(1){
				Controller.get(amdone, 0,PowTag);
				if(amdone != done){
					Controller.get(cTh, 0, PowTag+1);
					Controller.get(cPh, 0, PowTag+2);
					Controller.get(cGa, 0, PowTag+3);
					Controller.get(cW, 0, PowTag+4);
					fid_+=cW*singleStaticFID(ro, det,cTh, cPh, cGa);
					Controller.put(amdone, 0);
				}else{
					break;
				}
			}
		}

//Serial Implimentation
	}else{
		while(myit)
		{
			fid_+=myit.weight()*singleStaticFID(ro,det,myit.theta(), myit.phi(), myit.gamma());
			++myit;
		}
	}

	//get the final FID from all the procs
	Controller.reduce(fid_, Reduce::Add);
	//scatter it to all the procs
	Controller.scatter(fid_);

	return fid_;
}


template<class Func_T, class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
 oneFID<Func_T, MatrixType_T>::
  singleStaticFID( matrix_T &ro,  matrix_T2 &de, double th, double ph, double gamma)
{
	double dt=1./sw_;
	function_->setPowderAngles(th, ph, gamma);
	function_->setRotorAngles(0.0, rotb_*DEG2RAD);
	MatrixType_T hamil=function_->Hamiltonian(0.0, dt);
	Vector<ComplexType_T> fid(fid_.size(),ZeroType<ComplexType_T>::zero());
	fid.fill(0);
	ComplexType_T z(0.0, -dt*2.0*Pi);
 	PropMatType_T evect;
	_matrix<ComplexType_T, DiagonalMatrix> eval;
	diag(hamil, eval,evect);
	const double cutoff = 1.0e-10;
	PropMatType_T sig0=adjprop(evect,ro); //adjoint(evect)*ro*(evect);	// Next put sig0 into eigenbase of H
	PropMatType_T Do=adjprop(evect,de); //adjoint(evect)*de*(evect);				// Put detection op. to eigenbase of H

	int hs = hamil.rows();
	int ls = hs*hs;
	eval=exp(z*eval);
 	complex *A=new complex[ls];
	complex *B=new complex[ls];
	int i,j,pos=0;
	for(i=0; i<hs; i++){
		for(j=0; j<hs; j++){
			A[pos] = Do(i,j)*sig0(j,i);
			B[pos] = conj_mul(eval(i),eval(j));
			if(square_norm(A[pos])>cutoff)	pos ++;
    	}
	}
	for(int k=0; k<npts(); k++){
		z = ZeroType<complex>::zero();
		for(int p=0; p<pos; p++){
			z += A[p];
			A[p] *= B[p];
		}
		fid(k)=z ;
	}
	delete [] A;
	delete [] B;
	return fid;
}


template<class Func_T,class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
  oneFID<Func_T, MatrixType_T>::
     spinningFID( powder &in, matrix_T &ro, matrix_T2 &det)
{
	powder::iterator myit(in);
	fid_.fill(0);

	double cTh, cPh, cGa, cW;
	int PowTag=30, done=-1, amdone=1;
//first make sure we are Parrellel
	if(Controller.parallel()  && in.size()>1){
	//Perform the "master/slave" loops
		if(Controller.master()){
		//First we see if we send an initial request to each
			for(int rank=1;rank<Controller.size();++rank){
				Controller.put(amdone, rank, PowTag);
				Controller.put(myit.theta(), rank, PowTag+1);
				Controller.put(myit.phi(), rank, PowTag+2);
				Controller.put(myit.gamma(), rank, PowTag+3);
				Controller.put(myit.weight(), rank, PowTag+4);
				++myit;
				if(!myit) break;
			}
			int get, dd;
			while(myit){
				get=Controller.getAny(dd);
				Controller.put(amdone, get, PowTag);
				Controller.put(myit.theta(), get, PowTag+1);
				Controller.put(myit.phi(), get, PowTag+2);
				Controller.put(myit.gamma(), get, PowTag+3);
				Controller.put(myit.weight(), get, PowTag+4);
				++myit;
			}

			for(int rank=1;rank<Controller.size();++rank)
				Controller.put(done, rank, PowTag);

		//Slave
		}else{
			while(1){
				Controller.get(amdone, 0,PowTag);
				if(amdone != done){
					Controller.get(cTh, 0, PowTag+1);
					Controller.get(cPh, 0, PowTag+2);
					Controller.get(cGa, 0, PowTag+3);
					Controller.get(cW, 0, PowTag+4);
					fid_+=cW*singleSpinningFID(ro, det,cTh, cPh, cGa);
					Controller.put(amdone, 0);
				}else{
					break;
				}
			}
		}

//Serial Implimentation
	}else{
		while(myit)
		{
			fid_+=myit.weight()*singleSpinningFID(ro,det,myit.theta(), myit.phi(), myit.gamma());
			++myit;
		}
	}

	//get the final FID from all the procs
	Controller.reduce(fid_, Reduce::Add);
	//scatter it to all the procs
	Controller.scatter(fid_);

	return fid_;
}

template<class Func_T, class MatrixType_T>
MatrixType_T oneFID<Func_T, MatrixType_T>::Hamiltonian( double t1, double t2,double wrIn)
{
	function_->setRotorAngles(PI2*wrIn*(t1+(t2-t1)/2.0), rotb_*DEG2RAD);
	return function_->Hamiltonian(t1,t2,wrIn);
}

template<class Func_T, class MatrixType_T>
template<class matrix_T, class matrix_T2>
Vector<typename oneFID<Func_T, MatrixType_T>::ComplexType_T>
  oneFID<Func_T, MatrixType_T>::
    singleSpinningFID(matrix_T &ro, matrix_T2 &det,double th, double ph, double gamma)
{
	function_->setRotorAngles(0.0, rotb_*DEG2RAD);
	function_->setPowderAngles(th, ph, gamma);
	//Vector<complex> fid(npts(), 0);
	compute<Func_T, MatrixType_T> loo(*(function_), wr_, sw_, startT_, startT_+1.0/wr_);
	sw_=loo.sw();
	return loo.FID(ro,det, npts());
}


END_BL_NAMESPACE



#endif



