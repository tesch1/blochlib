
/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 09-25-01
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
	compute.h-->
this little class develops stroposcopically observed
 *spectra using the 'COMPUTE' method given in
 *JMR 120, 56 (1996) Eden, M, Lee, Y, and Levitt, Malcolm
 *it calcualtes a single propogatr for some modulation
 *period, T.  It uses all the little compute_step used to calculate
 *propogator to reconstruct the entire frequecy range
 *and thus a fid from the frequencies
 *
 *it also calculates propogators in a piece by peice method
 *i.e. U(t)=Prod(exp(-i dt H(t)))

 the 'function_t' class MUST have a function called
 'hmatrix Hamiltonian(double TIME1, double TIME2, double WR)'
  where 'TIME1=the begining of a delta T step
        'TIME2=the END of a delta T step
        'WR'= the Rotor Speed
 */

#ifndef _compute_h_
#define _compute_h_ 1



#include "container/matrix/matrix.h"
#include "container/Vector/Vector.h"

BEGIN_BL_NAMESPACE


/***** MatrixType_T MUST BE a COMPLEX matrix of some kind.. ***/

template<class function_t,class MatrixType_T=hmatrix>
class compute {
public:

	typedef  Complex<typename MatrixType_T::numtype::numtype> ComplexType_T;
	typedef _matrix<ComplexType_T, FullMatrix> PropMatType_T ;


private:
//	Vector<double> time;
//	static Vector<matrix> Ss;



	static Vector<PropMatType_T > As;
	function_t *mf;
	int rosym;		//1 if ro==1/2(det+adjoint(det)), 0=false, 2=not calculated YET
	template<class matrix_T, class matrix_T2>
	int isroSYM(const matrix_T &ro,const matrix_T2 &det)
	{
		if(rosym==2){
			if(ro==0.5*(det+adjoint(det))) return 1;
			else return 0;
		}else{
			return rosym;
		}
	}


	void CalcComputeStep()
	{
		if(wr_==0.0) return;
		compute_step=int(floor(sw_/wr_+0.5));
		if(gamma_step>=compute_step){
			gammaloop=gamma_step/compute_step;
		}
		if(gammaloop<1) gammaloop=1;
		gamma_step=gammaloop*compute_step;
	//	compute_time=1./(double(compute_step)*wr_);
		sw_=double(compute_step*wr_);
	}

public:

	PropMatType_T Uf;
	int pmax;
	int compute_step;
	int gamma_step;
	int gammaloop;
	double sw_,wr_;
	double tmin;
	double tmax;
	double tau;

	compute();
	compute(function_t &);
	compute(function_t &, int pe);
	compute(function_t &, double wr,double sw, double tmin, double tmax);
	compute(function_t &, int pe, double tmin, double tmax);

	~compute()
	{
		mf=NULL;
	}

	inline double wr()const{	return wr_;	}
	inline double RotorSpeed()const{	return wr_;	}

	inline void SetWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void SetRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}
	inline void setWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void setRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}

	inline double sw()const{	return sw_;	}
	inline double sweep()const{	return sw_;	}
	inline double SweepWidth()const{	return sw_;	}
	inline double sweepWidth()const{	return sw_;	}

	inline void SetSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}

	inline int GammaStep() const {	return gamma_step;	}
	inline int gammaStep() const {	return gamma_step;	}
	inline void setGammaStep(int in){	SetGammaStep(in);	}
	inline void SetGammaStep(int in)
	{
		RunTimeAssert(in>=1);
		gamma_step=in;
		CalcComputeStep();
	}

	inline void calcU(){	calcU(Uf);	}

	template<class matrix_T2>
	void calcU(matrix_T2 &);
	void calcUFID(){	calcUFID(1);	}
	void calcUFID(int gammaon);

	template<class matrix_T, class matrix_T2>
	Vector<ComplexType_T> FID(matrix_T &ro, matrix_T2 &det, int npts);

	template<class matrix_T2, class matrix_T3>
	Vector<ComplexType_T> calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts);

	template<class matrix_T, class matrix_T2>
	inline Vector<ComplexType_T> calcFID( matrix_T &ro, matrix_T2 &det, int npts)
	{		return calcFID(Uf, ro,det,npts);		}


};

template<class function_t,class MatrixType_T>
Vector<typename compute<function_t,MatrixType_T>::PropMatType_T> compute<function_t, MatrixType_T>::As(1,PropMatType_T());

template<class function_t,class MatrixType_T>
compute<function_t, MatrixType_T>::compute()
{
	mf=NULL;
	compute_step=0;pmax=0;
	tmin=0.; tmax=0.; tau=0.;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t,class MatrixType_T>
compute<function_t, MatrixType_T>::compute(function_t &in)
{
	mf=&in;
	tmin=0.;compute_step=0;
	tmax=0.; tau=0.;pmax=0;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t,class MatrixType_T>
compute<function_t, MatrixType_T>::compute(function_t &in, int piin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1, mf->Fe());
	Uf=mf->Fe();
	pmax=compute_step+1;
	tmin=0.;tmax=0.; tau=0.;
	rosym=2;
}

template<class function_t,class MatrixType_T>
compute<function_t, MatrixType_T>::compute(function_t &in, double wr,double sw, double tminin, double tmaxin)
{
	mf=&in;
	wr_=wr;
	sw_=sw;
	gammaloop=1;
	gamma_step=10;
	CalcComputeStep();
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	//time.resize(compute_step+1, 0);
	//wmod=2.*Pi/(tmax-tmin);
	rosym=2;

}

template<class function_t,class MatrixType_T>
compute<function_t, MatrixType_T>::compute(function_t &in, int piin, double tminin, double tmaxin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	gammaloop=1;
	gamma_step=10;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	rosym=2;
}


template<class function_t,class MatrixType_T>
void compute<function_t, MatrixType_T>::calcUFID(int gammaon)
{
	double tadd=PI2*double(gammaon)/double(gammaloop*compute_step)/wr_;
	double t1=tmin+tadd;
	double t2=tmin+tadd+tau;
	static PropMatType_T hh;
	for(int i=0;i<compute_step;i++){
		hh=Mexp(mf->Hamiltonian(t1, t2, wr_), -ComplexType_T(0,1)*tau*PI2);
		if(i==0){ As[0].identity(hh.rows());}
		As[i+1]=hh*As[i];
		//std::cout<<As[i+1]<<std::endl;
		t1+=tau;
		t2+=tau;

	}
	Uf=(As[compute_step]);
}


template<class function_t,class MatrixType_T>
template<class matrix_T2>
void compute<function_t, MatrixType_T>::calcU(matrix_T2 &in)
{
	double t1=tmin;
	double t2=tmin+tau;
	static PropMatType_T hh;
	for(int i=0;i<compute_step;i++){
		hh=Mexp(mf->Hamiltonian(t1, t2,wr_), ComplexType_T(0,1)*tau*PI2);
		//hh=Mexp(hh);
		in*=hh;
		t1+=tau;
		t2+=tau;
	}
	//std::cout<<"U"<<std::endl<<in<<std::endl;
}

template<class function_t,class MatrixType_T>
template<class matrix_T,class matrx_T2>
Vector<typename compute<function_t,MatrixType_T>::ComplexType_T>
    compute<function_t, MatrixType_T>::FID(matrix_T &ro, matrx_T2 &det, int npts)
{
	//std::cout<<gammaloop<<" "<<gamma_step<<" "<<compute_step<<" "<<wr_<<" "<<sw_<<std::endl;
	Vector<ComplexType_T> fid(npts, 0);
	for(int q=0;q<gammaloop;q++){
		calcUFID(q);
		fid+=calcFID( ro, det, npts);
	}
	return fid;
}

template<class function_t,class MatrixType_T>
template<class matrix_T2, class matrix_T3>
Vector<typename compute<function_t,MatrixType_T>::ComplexType_T>
   compute<function_t, MatrixType_T>::calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts)
{
	Vector<ComplexType_T> fid(npts,ZeroType<complex>::zero());
	rosym=isroSYM(ro, det);
	PropMatType_T evect;
	_matrix<ComplexType_T, DiagonalMatrix> eval;
	diag(Uin, eval, evect);
	int N=Uin.rows();
	int i=0, j=0,   p=0, r=0, s=0;
	static Vector<ComplexType_T> ev;
	ev.resize(N, ZeroType<ComplexType_T>::zero());
	static matrix wrs;
	wrs.resize(N, N);
	double tau2PI=tau;
	//the transition matrix
	ComplexType_T tott=ComplexType_T(0.,double(compute_step)*tau2PI);
	for(i=0;i<N;i++){
		ev[i]=log(eval(i,i));
	}

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			wrs(i,j)=chop((ev[i]-ev[j])/tott, 1e-10);
		}
	}

	//the gamma-compute algorithm use the symetry relation between the gamma powder
	//angles DIVIDED into 2Pi/compute_step sections...the detection operator
	//so for each gamma angle we would think that we would need (compute_step) propogators
	//for each rotor cycle division (which we set to be equal to (compute_step) also) to
	//corispond to each different gamma angle...well not so, becuase we have some
	//nice time symetry between the gamma angle and teh time eveolved we only
	//need to calculate the propogators for gamma=0...from this one in select
	//combinations we can generate all the (compute_step) propogators from the gamma=0 ones
	//we still need to divide our rotor cycle up however, also into (compute_step) propogators
	//for the remaining notes i will use this labaling convention
	//
	//  pQs_i --> the transformed detection operator for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)
	//
	//	pRo --> the transformed density matrix for the for a gamma angle of 'p'*2Pi/(compute_step)
	//			NOTE:: the rotor division info is contained in the Qs
	//
	//	pAs_i --> the unitary trasformation for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)...the operators 0As_i were calculated
	//			  in the function 'calcUFID(int)'
	//
	//  the 'pth' one of all of these is realated back to the '0th' operator by some series
	//  of internal multiplications
	//

	ComplexType_T tmp1;
	static Vector<PropMatType_T> Qs;
	Qs.resize(compute_step);
	static Vector<PropMatType_T> rotmp;
	rotmp.resize(compute_step);


	//calculating
	// the ith density matrix (0rotmp_i)^d=event_inv*(0As_i)^d*ro*(0As_i)*evect
	// the ith detection op for (0Qs_i)^d=vent_inv*(0As_i)^d*ro*(0As_i)*evect
	//
	//IF ro=1/2(det+adjoint(det)) then
	// the ith density matrix is (0rotmp_i)^d=0Qs_i+(0QS_i)^d
	//
	//   the '^d' is a adjoint operation

	for(i=0;i<compute_step;i++){
	//	matrix loo=adjoint(As[i+1]);
		Qs[i]=adjprop(evect,adjprop(As[i+1], det));// evect_inv*loo*det*As[i+1]*evect; //Qs[i].chopper(1.e-10);
		if(rosym==0) rotmp[i]=adjprop(evect, adjprop(As[i+1], ro));//evect_inv*loo*ro*As[i+1]*evect; //rotmp[i].chopper(1.e-10);
		else	rotmp[i]=Qs[i]+adjoint(Qs[i]);
	}

	//The signal is then a nice sum over the transition matrix and the pQs's and pRos
	// Of course this is where we manipulate the '0th' operators to create the 'pth'
	//and combine them all into a 'f' matrix which contains the amplitudes
	//
	// pF_i(r,s)= means the 'pth' gamma angle for the 'ith' rotor division element (r,s)
	//
	// pF_i(r,s)=exp[i m wrs(r,s) tott] * (p%compute_step)Ro(s,r) * 0Qs_[i+p%compute_step](r,s) exp[-i wrs(r,s) j tott]
	//
	// here m=int((i+p)/compute_step)) -int(p/compute_step)

	//of course we have many 'p' sections the best way to calculate the F's is to average them
	//
	// Fave_i(r,s)=1/compute_step * Sum_(p=0)^(p=n-1) { pF_i(r,s) }
	// 		= 1/(compute_step)*Sum(...) { [(p%compute_step)R](s,r) exp(i [p-int(p/compute_step)*compute_step] wrs(r,s) tau) }*
	//								[0Qs_(i+p%n)(r,s)] exp(-i [j+p-int((j+p)/compute_step))*compute_step] wrs(r,s) tau)


	static Vector<PropMatType_T> Fave;
	Fave.resize(compute_step, mf->Fz());

	int ind1, ind2;
	for(i=0;i<compute_step;i++){
		for(p=0;p<compute_step;p++){
			ind1=(i+p)%(compute_step);		//for Q selection
			if((i+p)>=compute_step){
				ind2=p-compute_step;
			}else{
				ind2=p;
			}
			tmp1=ComplexType_T(0.,-double(ind2)*tau2PI);
			PropMatType_T tmm(N,N);
			for(r=0;r<N;r++){
				for(s=0;s<N;s++){
					//Fave[p]+=Qs[ind1]*transpose(rotmp[i])*exp(tmp1*wrs);
					Fave[p](r,s)+=Qs[ind1](r,s)*rotmp[i](s,r)*exp(tmp1*wrs(r,s));
					//tmm(r,s)=Qs[ind1](r,s)*rotmp[i](s,r);
				}
			}
			//std::cout<<"transpose test:"<<std::endl<<Qs[ind1]*transpose(rotmp[i])<<std::endl<<"tmm:" <<std::endl<<tmm<<std::endl;

			//std::cout<<"F: "<<p<<std::endl<<Fave[p]<<std::endl;
		}
	}
	ComplexType_T tmp=ComplexType_T(0.,(tau2PI));  //i*tau

	//a little computation time save...calculate the exp(i*tau*wrs) once
	//tmpmm;
	wrs=exp(tmp*wrs);
	PropMatType_T tmpmm(wrs);//(wrs);
	//the FID at intervals of i*tau is then given by
	//
	//  s(i*tau)=Sum_(r,s) { Fave_i(r,s)*exp[i wrs(r,s) i* tau] }

	//here is the j=0 point...saves us a exp calculation

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fid[0]+=Fave[0](i,j);
		}
	}

	for(i=1;i<npts;i++){
		p=i%compute_step;		//to select the Fave_p
		for(r=0;r<N;++r){
			for(s=0;s<N;++s){
				fid[i]+=Fave[p](r,s)*tmpmm(r,s);
				tmpmm(r,s)*=wrs(r,s); //advance 'time' exp(dt*wij)
			}
		}
		//cout<<i<<" "<<1e8*max(Re(tmpmm))<<" "<<1e8*max(Re(wrs))<<" "<<fid(i)<<std::endl;

		//std::cout<<"I: "<<i<<std::endl<<tmpmm<<std::endl<<"fid(i):"<<fid(i)<<std::endl;
	}
	int ff= rosym==1?2:1;
	fid*=double(1./double(compute_step/ff));
	return fid;
}


/***************** DIAGONAL MATRIX VERSION **************/



template<class function_t>
class compute<function_t, dmatrix> {

public:

	typedef Complex<double> ComplexType_T;
	typedef dmatrix PropMatType_T ;



private:
//	Vector<double> time;
//	static Vector<matrix> Ss;
	static Vector<PropMatType_T> As;
	function_t *mf;
	int rosym;		//1 if ro==1/2(det+adjoint(det)), 0=false, 2=not calculated YET
	template<class matrix_T, class matrix_T2>
	int isroSYM(const matrix_T &ro,const matrix_T2 &det)
	{
		if(rosym==2){
			if(ro==0.5*(det+adjoint(det))) return 1;
			else return 0;
		}else{
			return rosym;
		}
	}


	void CalcComputeStep()
	{
		if(wr_==0.0) return;
		compute_step=int(floor(sw_/wr_+0.5));
		if(gamma_step>=compute_step){
			gammaloop=gamma_step/compute_step;
		}
		if(gammaloop<1) gammaloop=1;
		gamma_step=gammaloop*compute_step;
	//	compute_time=1./(double(compute_step)*wr_);
		sw_=double(compute_step*wr_);
	}

public:
	PropMatType_T Uf;
	int pmax;
	int compute_step;
	int gamma_step;
	int gammaloop;
	double sw_,wr_;
	double tmin;
	double tmax;
	double tau;

	compute();
	compute(function_t &);
	compute(function_t &, int pe);
	compute(function_t &, double wr,double sw, double tmin, double tmax);
	compute(function_t &, int pe, double tmin, double tmax);

	~compute()
	{
		mf=NULL;
	}

	inline double wr()const{	return wr_;	}
	inline double RotorSpeed()const{	return wr_;	}

	inline void SetWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void SetRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}
	inline void setWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void setRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}

	inline double sw()const{	return sw_;	}
	inline double sweep()const{	return sw_;	}
	inline double SweepWidth()const{	return sw_;	}
	inline double sweepWidth()const{	return sw_;	}

	inline void SetSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}

	inline int GammaStep() const {	return gamma_step;	}
	inline int gammaStep() const {	return gamma_step;	}
	inline void setGammaStep(int in){	SetGammaStep(in);	}
	inline void SetGammaStep(int in)
	{
		RunTimeAssert(in>=1);
		gamma_step=in;
		CalcComputeStep();
	}

	inline void calcU(){	calcU(Uf);	}
	void calcU(PropMatType_T &);
	void calcUFID(){	calcUFID(1);	}
	void calcUFID(int gammaon);

	template<class matrix_T, class matrix_T2>
	Vector<complex> FID(matrix_T &ro, matrix_T2 &det, int npts);

	template< class matrix_T2, class matrix_T3>
	Vector<complex> calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts);

	template<class matrix_T, class matrix_T2>
	inline Vector<complex> calcFID( matrix_T &ro, matrix_T2 &det, int npts)
	{		return calcFID(Uf, ro,det,npts);		}


};

template<class function_t>
Vector<dmatrix> compute<function_t, dmatrix>::As(1,dmatrix());

template<class function_t>
compute<function_t, dmatrix>::compute()
{
	mf=NULL;
	compute_step=0;pmax=0;
	tmin=0.; tmax=0.; tau=0.;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t>
compute<function_t, dmatrix>::compute(function_t &in)
{
	mf=&in;
	tmin=0.;compute_step=0;
	tmax=0.; tau=0.;pmax=0;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t>
compute<function_t, dmatrix>::compute(function_t &in, int piin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1);
	As.fill(mf->Fe());
	Uf=mf->Fe();
	pmax=compute_step+1;
	tmin=0.;tmax=0.; tau=0.;
	rosym=2;
}

template<class function_t>
compute<function_t, dmatrix>::compute(function_t &in, double wr,double sw, double tminin, double tmaxin)
{
	mf=&in;
	wr_=wr;
	sw_=sw;
	gammaloop=1;
	gamma_step=10;
	CalcComputeStep();
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	//time.resize(compute_step+1, 0);
	//wmod=2.*Pi/(tmax-tmin);
	rosym=2;

}

template<class function_t>
compute<function_t, dmatrix>::compute(function_t &in, int piin, double tminin, double tmaxin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	gammaloop=1;
	gamma_step=10;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	rosym=2;
}


template<class function_t>
void compute<function_t, dmatrix>::calcUFID(int gammaon)
{
	double tadd=PI2*double(gammaon)/double(gammaloop*compute_step)/wr_;
	double t1=tmin+tadd;
	double t2=tmin+tadd+tau;
	static PropMatType_T hh;
	for(int i=0;i<compute_step;i++){
		hh=exp(mf->Hamiltonian(t1, t2, wr_)*(-complexi*tau*PI2));
		if(i==0){ As[0].identity(hh.rows());}
		As[i+1]=hh*As[i];
		//std::cout<<As[i+1]<<std::endl;
		t1+=tau;
		t2+=tau;

	}
	Uf=(As[compute_step]);
}


template<class function_t>
void compute<function_t, dmatrix>::calcU(PropMatType_T &in)
{
	double t1=tmin;
	double t2=tmin+tau;
	PropMatType_T hh;
	for(int i=0;i<compute_step;i++){
		hh=exp(mf->Hamiltonian(t1, t2,wr_)*complexi*tau*PI2);
		//hh=Mexp(hh);
		in*=hh;
		t1+=tau;
		t2+=tau;
	}
	//std::cout<<"U"<<std::endl<<in<<std::endl;
}

template<class function_t>
template<class matrix_T, class matrix_T2>
Vector<complex> compute<function_t, dmatrix>::FID(matrix_T &ro, matrix_T2 &det, int npts)
{
	//std::cout<<gammaloop<<" "<<gamma_step<<" "<<compute_step<<" "<<wr_<<" "<<sw_<<std::endl;
	Vector<complex> fid(npts, 0);
	for(int q=0;q<gammaloop;q++){
		calcUFID(q);
		fid+=calcFID( ro, det, npts);
	}
	return fid;
}

template<class function_t>
template<class matrix_T2, class matrix_T3>
Vector<complex> compute<function_t, dmatrix>::calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts)
{
	Vector<complex> fid(npts,ZeroType<complex>::zero());
	rosym=isroSYM(ro, det);
	dmatrix evect=mf->Fe();
	//dmatrix eval=exp(Uin);
	//diag(Uin, eval, evect);
	int N=Uin.rows();
	int i=0, j=0,   p=0, r=0, s=0;
	static Vector<complex> ev;
	ev.resize(N, ZeroType<complex>::zero());
	static matrix wrs;
	wrs.resize(N, N);
	double tau2PI=tau;
	//the transition matrix
	complex tott=complex(0.,double(compute_step)*tau2PI);
	for(i=0;i<N;i++){
		ev[i]=log(Uin(i,i));
	}

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			wrs(i,j)=chop((ev[i]-ev[j])/tott, 1e-10);
		}
	}

	//the gamma-compute algorithm use the symetry relation between the gamma powder
	//angles DIVIDED into 2Pi/compute_step sections...the detection operator
	//so for each gamma angle we would think that we would need (compute_step) propogators
	//for each rotor cycle division (which we set to be equal to (compute_step) also) to
	//corispond to each different gamma angle...well not so, becuase we have some
	//nice time symetry between the gamma angle and teh time eveolved we only
	//need to calculate the propogators for gamma=0...from this one in select
	//combinations we can generate all the (compute_step) propogators from the gamma=0 ones
	//we still need to divide our rotor cycle up however, also into (compute_step) propogators
	//for the remaining notes i will use this labaling convention
	//
	//  pQs_i --> the transformed detection operator for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)
	//
	//	pRo --> the transformed density matrix for the for a gamma angle of 'p'*2Pi/(compute_step)
	//			NOTE:: the rotor division info is contained in the Qs
	//
	//	pAs_i --> the unitary trasformation for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)...the operators 0As_i were calculated
	//			  in the function 'calcUFID(int)'
	//
	//  the 'pth' one of all of these is realated back to the '0th' operator by some series
	//  of internal multiplications
	//

	complex tmp1;
	static Vector<matrix> Qs(compute_step);
	Qs.resize(compute_step);
	static Vector<matrix> rotmp(compute_step);
	rotmp.resize(compute_step);


	//calculating
	// the ith density matrix (0rotmp_i)^d=event_inv*(0As_i)^d*ro*(0As_i)*evect
	// the ith detection op for (0Qs_i)^d=vent_inv*(0As_i)^d*ro*(0As_i)*evect
	//
	//IF ro=1/2(det+adjoint(det)) then
	// the ith density matrix is (0rotmp_i)^d=0Qs_i+(0QS_i)^d
	//
	//   the '^d' is a adjoint operation

	for(i=0;i<compute_step;i++){
	//	matrix loo=adjoint(As[i+1]);
		Qs[i]=adjprop(evect,adjprop(As[i+1], det));// evect_inv*loo*det*As[i+1]*evect; //Qs[i].chopper(1.e-10);
		if(rosym==0) rotmp[i]=adjprop(evect, adjprop(As[i+1], ro));//evect_inv*loo*ro*As[i+1]*evect; //rotmp[i].chopper(1.e-10);
		else	rotmp[i]=Qs[i]+adjoint(Qs[i]);
	}

	//The signal is then a nice sum over the transition matrix and the pQs's and pRos
	// Of course this is where we manipulate the '0th' operators to create the 'pth'
	//and combine them all into a 'f' matrix which contains the amplitudes
	//
	// pF_i(r,s)= means the 'pth' gamma angle for the 'ith' rotor division element (r,s)
	//
	// pF_i(r,s)=exp[i m wrs(r,s) tott] * (p%compute_step)Ro(s,r) * 0Qs_[i+p%compute_step](r,s) exp[-i wrs(r,s) j tott]
	//
	// here m=int((i+p)/compute_step)) -int(p/compute_step)

	//of course we have many 'p' sections the best way to calculate the F's is to average them
	//
	// Fave_i(r,s)=1/compute_step * Sum_(p=0)^(p=n-1) { pF_i(r,s) }
	// 		= 1/(compute_step)*Sum(...) { [(p%compute_step)R](s,r) exp(i [p-int(p/compute_step)*compute_step] wrs(r,s) tau) }*
	//								[0Qs_(i+p%n)(r,s)] exp(-i [j+p-int((j+p)/compute_step))*compute_step] wrs(r,s) tau)


	static Vector<matrix> Fave;
	Fave.resize(compute_step, mf->Fz());

	int ind1, ind2;
	for(i=0;i<compute_step;i++){
		for(p=0;p<compute_step;p++){
			ind1=(i+p)%(compute_step);		//for Q selection
			if((i+p)>=compute_step){
				ind2=p-compute_step;
			}else{
				ind2=p;
			}
			tmp1=complex(0.,-double(ind2)*tau2PI);
			matrix tmm(N,N);
			for(r=0;r<N;r++){
				for(s=0;s<N;s++){
					//Fave[p]=Qs[ind1]*transpose(rotmp[i])*exp(tmp1*wrs);
					Fave[p](r,s)+=Qs[ind1](r,s)*rotmp[i](s,r)*exp(tmp1*wrs(r,s));
					//tmm(r,s)=Qs[ind1](r,s)*rotmp[i](s,r);
				}
			}
			//std::cout<<"transpose test:"<<std::endl<<Qs[ind1]*transpose(rotmp[i])<<std::endl<<"tmm:" <<std::endl<<tmm<<std::endl;

			//std::cout<<"F: "<<p<<std::endl<<Fave[p]<<std::endl;
		}
	}
	complex tmp=complex(0.,(tau2PI));  //i*tau

	//a little computation time save...calculate the exp(i*tau*wrs) once
	//tmpmm;
	wrs=exp(tmp*wrs);
	matrix tmpmm(wrs);//(wrs);
	//the FID at intervals of i*tau is then given by
	//
	//  s(i*tau)=Sum_(r,s) { Fave_i(r,s)*exp[i wrs(r,s) i* tau] }

	//here is the j=0 point...saves us a exp calculation

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fid[0]+=Fave[0](i,j);
		}
	}

	for(i=1;i<npts;i++){
		p=i%compute_step;		//to select the Fave_p
		for(r=0;r<N;++r){
			for(s=0;s<N;++s){
				fid[i]+=Fave[p](r,s)*tmpmm(r,s);
				tmpmm(r,s)*=wrs(r,s); //advance 'time' exp(dt*wij)
			}
		}
		//cout<<i<<" "<<1e8*max(Re(tmpmm))<<" "<<1e8*max(Re(wrs))<<" "<<fid(i)<<std::endl;

		//std::cout<<"I: "<<i<<std::endl<<tmpmm<<std::endl<<"fid(i):"<<fid(i)<<std::endl;
	}
	int ff= rosym==1?2:1;
	fid*=double(1./double(compute_step/ff));
	return fid;
}


/*** SINGLE PRECISION VERIOSAN OF DIAGONAL ***/

template<class function_t>
class compute<function_t, dmatrixs> {

public:

	typedef scomplex ComplexType_T;
	typedef dmatrixs PropMatType_T ;
private:
//	Vector<double> time;
//	static Vector<matrix> Ss;
	static Vector<dmatrixs> As;
	function_t *mf;
	int rosym;		//1 if ro==1/2(det+adjoint(det)), 0=false, 2=not calculated YET
	template<class matrix_T, class matrix_T2>
	int isroSYM(const matrix_T &ro,const matrix_T2 &det)
	{
		if(rosym==2){
			if(ro==0.5*(det+adjoint(det))) return 1;
			else return 0;
		}else{
			return rosym;
		}
	}


	void CalcComputeStep()
	{
		if(wr_==0.0) return;
		compute_step=int(floor(sw_/wr_+0.5));
		if(gamma_step>=compute_step){
			gammaloop=gamma_step/compute_step;
		}
		if(gammaloop<1) gammaloop=1;
		gamma_step=gammaloop*compute_step;
	//	compute_time=1./(double(compute_step)*wr_);
		sw_=double(compute_step*wr_);
	}

public:
	dmatrixs Uf;
	int pmax;
	int compute_step;
	int gamma_step;
	int gammaloop;
	double sw_,wr_;
	double tmin;
	double tmax;
	double tau;

	compute();
	compute(function_t &);
	compute(function_t &, int pe);
	compute(function_t &, double wr,double sw, double tmin, double tmax);
	compute(function_t &, int pe, double tmin, double tmax);

	~compute()
	{
		mf=NULL;
	}

	inline double wr()const{	return wr_;	}
	inline double RotorSpeed()const{	return wr_;	}

	inline void SetWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void SetRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}
	inline void setWr(double in){	 wr_=in;	CalcComputeStep(); }
	inline void setRotorSpeed(double in){	 wr_=in; CalcComputeStep();	}

	inline double sw()const{	return sw_;	}
	inline double sweep()const{	return sw_;	}
	inline double SweepWidth()const{	return sw_;	}
	inline double sweepWidth()const{	return sw_;	}

	inline void SetSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void SetSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSw(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweep(double in){	 sw_=in;	CalcComputeStep();	}
	inline void setSweepWidth(double in){	 sw_=in;	CalcComputeStep();	}

	inline int GammaStep() const {	return gamma_step;	}
	inline int gammaStep() const {	return gamma_step;	}
	inline void setGammaStep(int in){	SetGammaStep(in);	}
	inline void SetGammaStep(int in)
	{
		RunTimeAssert(in>=1);
		gamma_step=in;
		CalcComputeStep();
	}

	inline void calcU(){	calcU(Uf);	}
	void calcU(dmatrixs &);
	void calcUFID(){	calcUFID(1);	}
	void calcUFID(int gammaon);

	template<class matrix_T, class matrix_T2>
	Vector<scomplex> FID(matrix_T &ro, matrix_T2 &det, int npts);

	template<class matrix_T2, class matrix_T3>
	Vector<scomplex> calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts);

	template<class matrix_T, class matrix_T2>
	inline Vector<scomplex> calcFID( matrix_T &ro, matrix_T2 &det, int npts)
	{		return calcFID(Uf, ro,det,npts);		}


};

template<class function_t>
Vector<dmatrixs> compute<function_t, dmatrixs>::As(1,dmatrixs());

template<class function_t>
compute<function_t, dmatrixs>::compute()
{
	mf=NULL;
	compute_step=0;pmax=0;
	tmin=0.; tmax=0.; tau=0.;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t>
compute<function_t, dmatrixs>::compute(function_t &in)
{
	mf=&in;
	tmin=0.;compute_step=0;
	tmax=0.; tau=0.;pmax=0;
	rosym=2;
	gammaloop=1;
	gamma_step=10;
	wr_=0;
	sw_=0;
}

template<class function_t>
compute<function_t, dmatrixs>::compute(function_t &in, int piin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1);
	As.fill(mf->Fe());
	Uf=mf->Fe();
	pmax=compute_step+1;
	tmin=0.;tmax=0.; tau=0.;
	rosym=2;
}

template<class function_t>
compute<function_t, dmatrixs>::compute(function_t &in, double wr,double sw, double tminin, double tmaxin)
{
	mf=&in;
	wr_=wr;
	sw_=sw;
	gammaloop=1;
	gamma_step=10;
	CalcComputeStep();
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	//time.resize(compute_step+1, 0);
	//wmod=2.*Pi/(tmax-tmin);
	rosym=2;

}

template<class function_t>
compute<function_t, dmatrixs>::compute(function_t &in, int piin, double tminin, double tmaxin)
{
	mf=&in;
	compute_step=piin;
	As.resize(compute_step+1, mf->Fe());
	pmax=compute_step+1;
	gammaloop=1;
	gamma_step=10;
	Uf=mf->Fe();
	tmin=tminin;
	tmax=tmaxin;
	tau=(tmax-tmin)/compute_step;
	if(tau<=0){
		BLEXCEPTION(" your time for the propgator is negative")
	}
	rosym=2;
}


template<class function_t>
void compute<function_t, dmatrixs>::calcUFID(int gammaon)
{
	double tadd=PI2*double(gammaon)/double(gammaloop*compute_step)/wr_;
	double t1=tmin+tadd;
	double t2=tmin+tadd+tau;
	static dmatrixs hh;
	for(int i=0;i<compute_step;i++){
		hh=exp(mf->Hamiltonian(t1, t2, wr_)*(-scomplexi*tau*PI2));
		if(i==0){ As[0].identity(hh.rows());}
		As[i+1]=hh*As[i];
		//std::cout<<As[i+1]<<std::endl;
		t1+=tau;
		t2+=tau;

	}
	Uf=(As[compute_step]);
}


template<class function_t>
void compute<function_t, dmatrixs>::calcU(dmatrixs &in)
{
	double t1=tmin;
	double t2=tmin+tau;
	dmatrixs hh;
	for(int i=0;i<compute_step;i++){
		hh=exp(mf->Hamiltonian(t1, t2,wr_)*scomplexi*tau*PI2);
		//hh=Mexp(hh);
		in*=hh;
		t1+=tau;
		t2+=tau;
	}
	//std::cout<<"U"<<std::endl<<in<<std::endl;
}

template<class function_t>
template<class matrix_T, class matrix_T2>
Vector<scomplex>
 compute<function_t, dmatrixs>::
   FID(matrix_T &ro, matrix_T2 &det, int npts)
{
	//std::cout<<gammaloop<<" "<<gamma_step<<" "<<compute_step<<" "<<wr_<<" "<<sw_<<std::endl;
	Vector<ComplexType_T> fid(npts, 0);
	for(int q=0;q<gammaloop;q++){
		calcUFID(q);
		fid+=calcFID( ro, det, npts);
	}
	return fid;
}

template<class function_t>
template<class matrix_T2, class matrix_T3>
Vector<scomplex>
  compute<function_t, dmatrixs>::
    calcFID( PropMatType_T &Uin,  matrix_T2 &ro, matrix_T3 &det, int npts)
{
	Vector<ComplexType_T> fid(npts,ZeroType<ComplexType_T>::zero());
	rosym=isroSYM(ro, det);
	dmatrixs evect=mf->Fe();
	//dmatrix eval=exp(Uin);
	//diag(Uin, eval, evect);
	int N=Uin.rows();
	int i=0, j=0,   p=0, r=0, s=0;
	static Vector<ComplexType_T> ev;
	ev.resize(N, ZeroType<ComplexType_T>::zero());
	static matrixs wrs;
	wrs.resize(N, N);
	double tau2PI=tau;
	//the transition matrix
	ComplexType_T tott=ComplexType_T(0.,float(compute_step)*tau2PI);
	for(i=0;i<N;i++){
		ev[i]=log(Uin(i,i));
	}

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			wrs(i,j)=chop((ev[i]-ev[j])/tott, 1e-10);
		}
	}

	//the gamma-compute algorithm use the symetry relation between the gamma powder
	//angles DIVIDED into 2Pi/compute_step sections...the detection operator
	//so for each gamma angle we would think that we would need (compute_step) propogators
	//for each rotor cycle division (which we set to be equal to (compute_step) also) to
	//corispond to each different gamma angle...well not so, becuase we have some
	//nice time symetry between the gamma angle and teh time eveolved we only
	//need to calculate the propogators for gamma=0...from this one in select
	//combinations we can generate all the (compute_step) propogators from the gamma=0 ones
	//we still need to divide our rotor cycle up however, also into (compute_step) propogators
	//for the remaining notes i will use this labaling convention
	//
	//  pQs_i --> the transformed detection operator for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)
	//
	//	pRo --> the transformed density matrix for the for a gamma angle of 'p'*2Pi/(compute_step)
	//			NOTE:: the rotor division info is contained in the Qs
	//
	//	pAs_i --> the unitary trasformation for the ith rotor division for a
	//			  gamma angle of 'p'*2Pi/(compute_step)...the operators 0As_i were calculated
	//			  in the function 'calcUFID(int)'
	//
	//  the 'pth' one of all of these is realated back to the '0th' operator by some series
	//  of internal multiplications
	//

	ComplexType_T tmp1;
	static Vector<matrixs> Qs(compute_step);
	Qs.resize(compute_step);
	static Vector<matrixs> rotmp(compute_step);
	rotmp.resize(compute_step);


	//calculating
	// the ith density matrix (0rotmp_i)^d=event_inv*(0As_i)^d*ro*(0As_i)*evect
	// the ith detection op for (0Qs_i)^d=vent_inv*(0As_i)^d*ro*(0As_i)*evect
	//
	//IF ro=1/2(det+adjoint(det)) then
	// the ith density matrix is (0rotmp_i)^d=0Qs_i+(0QS_i)^d
	//
	//   the '^d' is a adjoint operation

	for(i=0;i<compute_step;i++){
	//	matrix loo=adjoint(As[i+1]);
		Qs[i]=adjprop(evect,adjprop(As[i+1], det));// evect_inv*loo*det*As[i+1]*evect; //Qs[i].chopper(1.e-10);
		if(rosym==0) rotmp[i]=adjprop(evect, adjprop(As[i+1], ro));//evect_inv*loo*ro*As[i+1]*evect; //rotmp[i].chopper(1.e-10);
		else	rotmp[i]=Qs[i]+adjoint(Qs[i]);
	}

	//The signal is then a nice sum over the transition matrix and the pQs's and pRos
	// Of course this is where we manipulate the '0th' operators to create the 'pth'
	//and combine them all into a 'f' matrix which contains the amplitudes
	//
	// pF_i(r,s)= means the 'pth' gamma angle for the 'ith' rotor division element (r,s)
	//
	// pF_i(r,s)=exp[i m wrs(r,s) tott] * (p%compute_step)Ro(s,r) * 0Qs_[i+p%compute_step](r,s) exp[-i wrs(r,s) j tott]
	//
	// here m=int((i+p)/compute_step)) -int(p/compute_step)

	//of course we have many 'p' sections the best way to calculate the F's is to average them
	//
	// Fave_i(r,s)=1/compute_step * Sum_(p=0)^(p=n-1) { pF_i(r,s) }
	// 		= 1/(compute_step)*Sum(...) { [(p%compute_step)R](s,r) exp(i [p-int(p/compute_step)*compute_step] wrs(r,s) tau) }*
	//								[0Qs_(i+p%n)(r,s)] exp(-i [j+p-int((j+p)/compute_step))*compute_step] wrs(r,s) tau)


	static Vector<matrixs> Fave;
	Fave.resize(compute_step, mf->Fz());

	int ind1, ind2;
	for(i=0;i<compute_step;i++){
		for(p=0;p<compute_step;p++){
			ind1=(i+p)%(compute_step);		//for Q selection
			if((i+p)>=compute_step){
				ind2=p-compute_step;
			}else{
				ind2=p;
			}
			tmp1=ComplexType_T(0.,-double(ind2)*tau2PI);
			matrixs tmm(N,N);
			for(r=0;r<N;r++){
				for(s=0;s<N;s++){
					//Fave[p]=Qs[ind1]*transpose(rotmp[i])*exp(tmp1*wrs);
					Fave[p](r,s)+=Qs[ind1](r,s)*rotmp[i](s,r)*exp(tmp1*wrs(r,s));
					//tmm(r,s)=Qs[ind1](r,s)*rotmp[i](s,r);
				}
			}
			//std::cout<<"transpose test:"<<std::endl<<Qs[ind1]*transpose(rotmp[i])<<std::endl<<"tmm:" <<std::endl<<tmm<<std::endl;

			//std::cout<<"F: "<<p<<std::endl<<Fave[p]<<std::endl;
		}
	}
	ComplexType_T tmp=ComplexType_T(0.,(tau2PI));  //i*tau

	//a little computation time save...calculate the exp(i*tau*wrs) once
	//tmpmm;
	wrs=exp(tmp*wrs);
	matrixs tmpmm(wrs);//(wrs);
	//the FID at intervals of i*tau is then given by
	//
	//  s(i*tau)=Sum_(r,s) { Fave_i(r,s)*exp[i wrs(r,s) i* tau] }

	//here is the j=0 point...saves us a exp calculation

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			fid[0]+=Fave[0](i,j);
		}
	}

	for(i=1;i<npts;i++){
		p=i%compute_step;		//to select the Fave_p
		for(r=0;r<N;++r){
			for(s=0;s<N;++s){
				fid[i]+=Fave[p](r,s)*tmpmm(r,s);
				tmpmm(r,s)*=wrs(r,s); //advance 'time' exp(dt*wij)
			}
		}
		//cout<<i<<" "<<1e8*max(Re(tmpmm))<<" "<<1e8*max(Re(wrs))<<" "<<fid(i)<<std::endl;

		//std::cout<<"I: "<<i<<std::endl<<tmpmm<<std::endl<<"fid(i):"<<fid(i)<<std::endl;
	}
	int ff= rosym==1?2:1;
	fid*=float(1./float(compute_step/ff));
	return fid;
}


END_BL_NAMESPACE


#endif

