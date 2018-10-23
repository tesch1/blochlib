


/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
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
	bs_meth.h-->Bulirsch-Stoer/Richard Extrapolation method diff eq solver methods

*/

#ifndef _bs_meth_h_
#define _bs_meth_h_ 1

#include "container/Vector/Vector.h"
#include "utils/utils.h"


BEGIN_BL_NAMESPACE

//this can change!!!  to suit the ease of the system
//easy integration can use a smaller array of steps
// Since this integrator uses adaptive step sizes, if you choose
// an array that is too small it should just up the step size...resulting in a similar
// time for the integration...however, it may save several seconds here and
// there if you use the smaller array sizes...However, if you choose an array too small
// the number of steps it increases will increase you total time...a bit of
// profiling is nessesary (look to the function "mmid" to get a good idea of the time)
// in the profile report....


//tougher equations (i.e. highly nonlinear & stiff type) should use
/*
static const int __nseqNotStiff[10]={0,2,4,6,8,10,12,14, 16, 18};

template<class func, class T, class Container >
const int  bs<func, T,  Container >::kmaxx=8;

template<class func, class T, class Container>
const int  bs<func, T,  Container >::imaxx=9;
*/
//'medium' (highly coupled) equations should use this

static const int __nseqNotStiff[8]={0,2,4,6,8,10,12,14};
template<class func, class T, class Container >
const int  bs<func, T,Container >::kmaxx=6;
template<class func, class T, class Container >
const int  bs<func, T,Container >::imaxx=7;

//and linear simple equations can go even farther down
/*
static const int __nseqNotStiff[6]={0,2,4,6,8,10};
template<class func, class T, class Container >
const int  bss<func, T,Container >::kmaxx=4;
template<class func, class T, class Container >
const int  bs<func, T,Container >::imaxx=5;
*/

template<class func, class T, class Container >
const double  bs<func, T,Container >::SC_max=1;

template<class func, class T, class Container >
Container bs<func, T,Container >::yn(1,0);  //to make compiler happy...they truely set in the set_rest_i() function...

template<class func, class T, class Container >
Container bs<func, T,Container >::ym(1,0); //to make compiler happy...they truely set in the set_rest_i() function...


template<class func, class T, class Container >
void bs<func, T,Container >::reset(){
	basicODE<T, Container>::reset();
	yn=ysav;
	ym=ysav;

	x.fill(0);
	d.resize(num_func, kmaxx+1, 0);
	found_order=false;
	findorder();
}

template<class func, class T, class Container >
void bs<func,T,Container >::set_yscal(){
	//yn=(y);		//just a dummy vector..temporary (avoid initialization step for it over and over again)
//	my_sys->function(cur_t, y, yn);
//	func_calls++;
//	yscal=tiny +abs(cur_step)*abs(yn)+abs(y);
}

/***** FINDING THE ORDER ****/
template<class func, class T, class Container  >
Vector<double> bs<func, T,Container >::a(imaxx+1); //used to detemin the max order

template<class func, class T, class Container  >
rmatrix bs<func, T,Container >::alf(imaxx+1, imaxx+1); //the work matrix

template<class func, class T, class Container  >
Vector<double>  bs<func, T,Container >::err(imaxx+1); //local error

template<class func, class T, class Container  >
void bs<func, T,Container >::findorder()
{
	kmax=0;
	kopt=0;
	stepnext=-1.0e29;
	double eps1=safty1*eps;
	a[1]=__nseqNotStiff[1]+1;
	for(int k=1;k<=kmaxx;k++){
		a[k+1]=a[k]+__nseqNotStiff[k+1];
	}
	for(int iq=2;iq<=kmaxx;iq++){
		for(int k=1;k<iq;k++){
			alf(k,iq)=pow(eps1, (a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.)*(2.*k+1)));
		}
	}
//figure out the maximum # of stiffbs sections to divide each step into
	for(kopt=2;kopt<kmaxx;kopt++){
		if(a[kopt+1]>a[kopt]*alf(kopt-1,kopt)) {	break; }
	}
	kmax=kopt;

	kmax=kopt;
	order=2*(kmax-1);
	found_order=true;
}

template<class func, class T, class Container >
void bs<func, T,Container >::mmid(Container &ysavO, double &step, double &nstep, Container &yseqO){

//varous temporary numbers i'll need
	double temp_2step=0, temp_step=0, temp_t=0;
	T swap;

//our two ends of the midpoint calc are "yn" and "ym"

//divide our tot steps into the little ones needed
	temp_step=step/nstep;

//calculate the f(x,t) at the first y (values held in yseqO)
	my_sys->function(cur_t, ysavO, yn);
	//cout<<"BEGIN: ["<<yn<<"]"<<endl;

//now step up the two ends
	ym=ysavO;
	yn=ysavO+temp_step*yn;


//advance our time by a smidge
	temp_t=cur_t+temp_step;

//recalcuate the far end point
	my_sys->function(temp_t, yn, yseqO);
	func_calls++;

//move two temp steps
	temp_2step=2.*temp_step;

//now we do a midpoint calc taking smaller and smaller divisions
//until we hit the max divisions, nstep
	for(int n=2;n<=nstep;n++){
		for(int i=0;i<num_func;i++){
			swap=ym[i]+temp_2step*yseqO[i];
			ym[i]=yn[i];
			yn[i]=swap;
		}
		temp_t+=temp_step;
		my_sys->function(temp_t, yn, yseqO);
		func_calls++;
	}

//assign the out value
	yseqO=0.5*(ym+yn+temp_step*yseqO);
	//cout<<"END: ["<<yseqO<<"]"<<endl;
}

//this is the polynomial extropoaltion rioutine
//it evaluates num_func functions at x=0 by fitting
//to a series of estimates with progressively smaller values x=xest
//and the function vectors 'yseq.' ('yseq' was given by the
//midpoint routine above). the output is put into
//'y' and the estimated error is put into 'yerr'

template<class func, class T, class Container  >
void bs<func, T,Container >::polyex(int &xiest, T &xest, Container &estim, Container &extrap)
{

	T q,f2,f1,delta;
	static Container c(num_func);

	x[xiest]=xest;

//set everythin equal at first
	yerr=estim;
	extrap=estim;

	if(xiest==1){
		for(int i=0;i<num_func;i++){
			d(i,1)=estim[i];
		}
	}else{
		c=estim;
		for(int k=1;k<xiest;k++){
			delta=1.0/(x[xiest-k]-xest);
			f1=xest*delta;
			f2=x[xiest-k]*delta;
			for(int j=0;j<num_func;j++){
				q=d(j,k);
				d(j,k)=yerr[j];
				delta=c[j]-q;
				yerr[j]=f1*delta;
				c[j]=f2*delta;
				extrap[j]+=yerr[j];
			}
		}
		for(int i=0;i<num_func;i++){
			d(i,xiest)=yerr[i];
		}
	}
}

template<class func, class T, class Container  >
void bs<func, T,Container >::bsstep(Container &yO, Container &yscalO)
{
	static int  xiest=0;
	static double epsold=-1.0, temp_t=0.;
	static T xest;
	static double temp_step=0;

//	const int init=imaxx+1;
//	static Vector<double> a(init);
//	static rmatrix alf(init,init);
//	static Vector<double> err(init);

//	static int reduct, exitflag=0;

	Container yseq; yseq=yO;
	//yseq.fill(0);
//first we set up the initial error testing structures a[] and alf[]
	if(!found_order || epsold!=eps){ findorder(); epsold=eps;	}


//set the step size to the one given initially
	temp_step=cur_step;

//	ysav=yO;

//	if(cur_t!=temp_t || temp_step !=stepnext){
//		first=1;
//		kopt=kmax;
//	}

//	reduct=0;
//	for(;;){
		//std::cerr<<"cur_t: "<<cur_t<<" cur_step: "<<cur_step<<" temp stp: "<<temp_step<<" ngood: "<<ngood<<" nbad: "<<nbad<<std::endl;
		for(int k=1;k<=kmax;k++){
			xiest=k;
			temp_t=cur_t+temp_step;
			
			if(temp_t==cur_t && temp_t!=t2 && temp_t!=t1 && temp_step==0.0){
				basicODE<T, Container>::error(std::cerr, "ERROR: bs::bsstep(...)::STEP SIZE UNDERFLOW");
				BLEXCEPTION("ERROR: bs::bsstep(...)::STEP SIZE UNDERFLOW")
				break;
			}
			double nseq_on=double(__nseqNotStiff[k]);
			mmid(ysav, temp_step, nseq_on, yseq);
			xest=T(pow(temp_step/nseq_on, 2.0));
			polyex(xiest, xest, yseq, yO);
		}
			/*
			if(k!=1){
				errmax=tiny;
				double ratio_tmp=0;
				for(int i=0;i<num_func;i++){
					ratio_tmp=max(abs(yerr[i])/abs(yscalO[i]));
					errmax=std::max(errmax, ratio_tmp);
				}
				errmax=errmax/eps;
				km=k-1;
				err[km]=pow(errmax/safty1, 1./(2.*(km)+1));
			}
			if(k!=1&&(k>=kopt-1 || first)){
				if(errmax<1.0){
					exitflag=1;
					break;
				}
				if(k==kmax||k==kopt+1){
					red=safty2/err[km];
					break;
				}else if(k==kopt && alf(kopt-1,kopt)<err[km]){
					red=1./err[km];
					break;
				}else if(kopt==kmax && alf(km,kmax-1)<err[km]){
					red=alf(km,kmax-1)*safty2/err[km];
					break;
				}else if(alf(km,kopt)<err[km]){
					red=min(alf(km,kopt-1))/err[km];
					break;
				}
			}
		}
		if(exitflag){
			break;
		}
		red=std::min(red, red_min);
		red=std::max(red, red_max);

		temp_step *= red;
		reduct=1;
	}
	cur_t=temp_t;
	stepdid=temp_step;
	first=0;
	wrkmin=1.0e35;
	for(int k=1;k<=km;k++){
		fact=std::max(err[k], scal_max);
		work=min(fact*a[k+1]);
		if(work<wrkmin){
			scale=fact;
			wrkmin=work;
			kopt=k+1;
		}
	}
	stepnext=temp_step/scale;

	if(kopt>=xiest && kopt!=kmax && !reduct){
		fact=std::max(scale/alf(kopt-1,kopt), scal_max);
		if(a[kopt+1]*fact<=wrkmin){
			stepnext=temp_step/fact;
			kopt++;
		}
	}
	*/
}

template<class func, class T, class Container >
void bs<func,T,Container >::odeint(){
//set the correect sign of the step (i.e. moving backwards or forwards in time)
	double diff_tmp=t2-t1;
	cur_step=(double(sign(diff_tmp))*stepi);

	int numstep=0;

//we do not want to go on forever so max sure we stop after max_step itterations
//ofstream logf("doo", ios::app);
	for(numstep=0;numstep<max_step;numstep++){

//this makes sure we do not aver shoot our final time 't2'
		if((cur_t+cur_step-t2)*(cur_t+cur_step-t1)>0.0){
			cur_step=(t2-cur_t);
		}
//calcualte the yscal factor for this step
		//set_yscal();
		ysav=y;
		bsstep(y, yscal);							//do a single step

		double D0=absToler+relToler*(max(max(abs(ysav)))+cur_step*max(max(abs(y))));
		double D1=max(max(abs(yerr)));
		double r= D1/D0;
		double ShrinkForce=1.0;
	//if something crazy happens (like a bad func eval)
	//make sure we catch it
#ifdef HAVE_ISNAN
		if(hasnan(r)|| hasnan(D1)|| D1>1.0  || hasnan(D0)){ r=2.; ShrinkForce=0.1;}
#endif
		//cout<<"Dt: "<<cur_step<<" r: "<<r<<" D0: "<<D0<<" D1: "<<D1<<endl;
		if(r>1.0/safty){ //stesp decrease...must repest our section
			stepnext=safty*cur_step*pow(ShrinkForce, 1.0/double(order+1))*ShrinkForce;
			cur_step=stepnext;
			stepdid=cur_step;
			nbad++;
			y=ysav;
		}else if(r<0.5){
			double tmf=safty*pow(1.0/r, 1.0/double(order));
			stepnext=cur_step*(tmf>1.0)?tmf:1.0;
			//make sure we do not got out of control in the step size increase
			if(hold_step)	stepnext=min(5.0*cur_step, stepnext);
			ngood++;
			//did a step nicely, so get out of the main loop
			stepdid=cur_step;
			cur_t+=stepdid;
			//logf<<cur_t<<" "<<y<<endl;

		}else{
			ngood++;
			stepdid=cur_step;
			cur_t+=stepdid;
			//logf<<cur_t<<" "<<y<<endl;
		}
		//if we have finished to t2 break out
		if((cur_t-t2)*(t2-t1)>=0){
			//cout<<"cur_t: "<<cur_t<<" t2: "<<t2<<endl;
			break;
		}

		if(abs(stepnext)<=stepmin){
			basicODE<T, Container>::error(std::cerr,"Error: bs.solve()..STEP SIZE IS TOO SMALL");
			BLEXCEPTION("Error: bs.solve()..STEP SIZE IS TOO SMALL")
		}
		cur_step= stepnext;
	}
	if(numstep==max_step){
		basicODE<T, Container> ::error(std::cerr, "Warning::bs::solve()--TOO MANY STEPS "+itost(numstep));
		BLEXCEPTION("Warning::bs::solve()--TOO MANY STEPS "+itost(numstep))		
	}
}


//the master 'vector' solver



END_BL_NAMESPACE

#endif
