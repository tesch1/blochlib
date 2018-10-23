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
	gear_meth.h-->a 5th order Gear Diff eq solver methods

*/

#ifndef _gear_meth_h_
#define _gear_meth_h_ 1

BEGIN_BL_NAMESPACE


template<class func, class T>
const double  gear<func, T>::AA=1.0;
template<class func, class T>
const double gear<func, T>::a0=AA*1.0;
template<class func, class T>
const double gear<func, T>::a1=AA*1.0;
template<class func, class T>
const double gear<func, T>::BB=1.0/3.0;
template<class func, class T>
const double gear<func, T>::b0=BB*4.0;
template<class func, class T>
const double gear<func, T>::b1=BB*(-1.0);
template<class func, class T>
const double gear<func, T>::b2=BB*2.0;
template<class func, class T>
const double  gear<func, T>::CC=1./11.0;
template<class func, class T>
const double gear<func, T>::c0=CC*18.0;
template<class func, class T>
const double gear<func, T>::c1=CC*(-9.0);
template<class func, class T>
const double gear<func, T>::c2=CC*2.0;
template<class func, class T>
const double gear<func, T>::c3=CC*6.0;
template<class func, class T>
const double  gear<func, T>::DD=1.0/25.0;
template<class func, class T>
const double gear<func, T>::d0=DD*48.0;
template<class func, class T>
const double gear<func, T>::d1=DD*(-36.0);
template<class func, class T>
const double gear<func, T>::d2=DD*(16.0);
template<class func, class T>
const double  gear<func, T>::d3=DD*(-3.0);
template<class func, class T>
const double gear<func, T>::d4=DD*(12.0);
template<class func, class T>
const double gear<func, T>::EE=1.0/137.0;
template<class func, class T>
const double  gear<func, T>::e0=EE*300.0;
template<class func, class T>
const double gear<func, T>::e1=EE*(-300.0);
template<class func, class T>
const double gear<func, T>::e2=EE*(200.0);
template<class func, class T>
 const double  gear<func, T>::e3=EE*(-75.0);
template<class func, class T>
const double gear<func, T>::e4=EE*12.0;
template<class func, class T>
const double gear<func, T>::e5=EE*+60.0;
template<class func, class T>
const double  gear<func, T>::FF=1.0/147.0;
template<class func, class T>
const double gear<func, T>::f0=360.0;
template<class func, class T>
const double  gear<func, T>::f1=FF*(-450.0);
template<class func, class T>
const double  gear<func, T>::f2=FF*400.0;
template<class func, class T>
const double gear<func, T>::f3=FF*(-225.0);
template<class func, class T>
const double gear<func, T>::f4=FF*(72.0);
template<class func, class T>
const double gear<func, T>::f5=FF*(-10.0);
template<class func, class T>
const double gear<func, T>::f6=FF*(60.0);


template<class func, class T>
double gear<func, T>::grow=-0.2;

template<class func, class T>
double gear<func, T>::shrink=-0.25;

template<class func, class T>
double gear<func, T>::safty=0.9;

template<class func, class T>
double gear<func, T>::err_max=1.89e-7;

template<class func, class T>
int gear<func, T>::max_step=1000000;

template<class func, class T>
int gear<func, T>::max_step_within=10000;

template<class func, class T>
double gear<func, T>::tiny=1.0e-30;

template<class func, class T>
double gear<func, T>::eps=1.0e-8;



template<class func, class T>
void gear<func, T>::set_rest_i(){
	yerr.resize(num_func,T(0));
	yscal.resize(num_func,T(0));
	hold_step=0;

}



/*Vector<T> gear::set_y_mat(matrix &iny){
	Vector<T> tmp;
	tmp=my_sys->array_from_mat(iny);
	return tmp;
}*/

template<class func, class T>
void gear<func, T>::set_yscal(){
	//Vector<T> tmp(num_func, T(0));
	//float floo=0;
	//my_sys->function(cur_t, y, tmp);
	/*for(int i=0;i<num_func;i++){
		if(abs(y[i])>=0.1)	{
			floo=abs(cur_step)*abs(tmp[i])+tiny+abs(y[i]);
			yscal[i]=floo;
		}else{
			floo=abs(cur_step)*abs(tmp[i])+tiny+abs(y[i]);
			yscal[i]=floo;
		}
	}*/
	static T  CC(1);
	for(int i=0;i<num_func;++i)
	{
		if(max(CC)>max(y[i]))
		{
			yscal[i]=CC;
		}else{
			yscal[i]=y[i];
		}
	}
}

//The MAIN function stepper....

template<class func, class T>
void gear<func, T>::gr(int stepct){

//here are some temp bins for the 'for' loop
	static Vector<T> ak1(num_func,0.);
	static Vector<T> ak2(num_func,0.);
	static Vector<T> ak3(num_func,0.);
	static Vector<T> ak4(num_func,0.);
	static Vector<T> ak5(num_func,0.);
	static Vector<T> ak6(num_func,0.);
	static Vector<T> ytemp(num_func,0.);
	static Vector<T> ytemp1(num_func,0.);
	static double tmt=cur_t;
	ytemp=y;

/*
	the intial points are typically the 'hard' ones and most inacuarte becuase we
	have NO history of the system eveolution, so we have to 'wait' to buld up the
	history...
	so it IMPORTANT that the inital step size entered by the user is small enough
	not to cause horrible problems...as we cannot really check the local error
	until we hit the 5th order step...

	the first step is first order and uses the 'A' coief
	the results is stored in 'ak1' and 'y' (our master variable)

	the second step uses the 'B' and is stored in 'ak2'
	....
	the 5th step is stored in ak5...

	after all 5 steps are obtained we can start the 'normal' progression
	and use the preivous pts for integration, moving the last pt to 'ak_-1'


	****THe construcutor of this class MAKES SURE that 'cur_step' < (tend-tbegin)/5
	so we can at least get the first first 5 back predict steps


*/

	//ak1
		if(stepct==0)
		{
			my_sys->function(tmt, ytemp, ak1);
			ak1=ytemp+cur_step*ak1;

			tmt+=cur_step;
			my_sys->function(tmt, ytemp, ak2);
			ak2=b0*ytemp+b1*ak1+b2*cur_step*ak2;

			tmt+=cur_step;
			my_sys->function(tmt, ytemp, ak3);
			ak3=c0*ytemp+c1*ak1+c2*ak2+c3*cur_step*ak3;

			tmt+=cur_step;
			my_sys->function(tmt, ytemp, ak4);
			ak4=d0*ytemp+d1*ak1+d2*ak2+d3*ak3+d4*cur_step*ak4;

			tmt+=cur_step;
			my_sys->function(tmt, ytemp, ak5);
			y=e0*ytemp+e1*ak1+e2*ak2+e3*ak3+e4*ak4+e5*cur_step*ak5;
			//our error step checker

		//	my_sys->function(tmt, ytemp, ak6);
		//	ak6=f0*ytemp+f1*ak1+f2*ak2+f3*ak3+f4*ak4+f5*y+f6*cur_step*ak6;
		//	yerr=ak6;

		}
		else
		{
			ak1=ak2;
			ak3=ak2;
			ak3=ak4;
			ak4=ak5;
			for(int i=0;i<5;i++){
				my_sys->function(cur_t+cur_step, ytemp, ak5);
				ak5=e0*ytemp+e1*ak1+e2*ak2+e3*ak3+e4*ak4+e5*cur_step*ak5;
			}
			y=ak5;

			//our error step checker
		//	my_sys->function(cur_t+2.0*cur_step, ytemp, ak6);
		//	ak6=f0*ytemp+f1*ak1+f2*ak2+f3*ak3+f4*ak4+f5*y+f6*cur_step*ak6;
		//	yerr=ak6;


		}





}

template<class func, class T>
void gear<func, T>::round1(){
	double temp_err=0;
	double temp_step=0;
	double temp_t=0;
	int count=0;
	int stepct=0;
	double tmpcurstep=cur_step;
	while(count<max_step_within){
		gr(stepct);
		tmpcurstep=cur_step;
		temp_err=0;
		/*for(int i=0;i<num_func;i++){
			float loof=yerr[i]/yscal[i];
			loof=abs(loof);
			temp_err=max(temp_err, loof);
		}*/

		temp_err=max(temp_err,max(abs(yerr/yscal)) );
		cout<<"ERR: "<<temp_err<<" "<<max(yerr)<<" "<<max(yscal)<<endl;
		temp_err=temp_err/eps;
		if(temp_err<=1.0)	break;

		temp_step=safty*cur_step*pow(temp_err, shrink);
		T temp_compar=0.1*cur_step;
		if(cur_step>=0){
			cur_step=max(temp_step, temp_compar);
		}else{
			cur_step=min(temp_step, temp_compar);
		}
		if(cur_step!=tmpcurstep) stepct=0;
		else stepct=1;
		temp_t=cur_t+cur_step;
		if(temp_t==cur_t&&cur_step==0.){
			std::cout<<"STEP SIZE UNDERFLOW!!"<<std::endl;
			std::cout<<"moving the step size so that i may continue"<<std::endl;
			std::cout<<"cur_t: "<<cur_t<<" cur_step: "<<cur_step<<" temp_step: "<<temp_step<<std::endl;
			temp_t*=2+tiny;
			break;
		}
		cur_t=temp_t;
		count++;
	}
	if(temp_err>err_max){
		stepnext=safty*cur_step*pow(temp_err, grow);
	}else{
		stepnext=5.0*cur_step;
	}
	stepdid=cur_step;
	cur_t=cur_t+stepdid;
}

template<class func, class T>
void gear<func, T>::odeint(){
//set the correect sign of the step (i.e. moving backwards or forwards in time)
	float diff=t2-t1;
	cur_step=sign(diff)*cur_step;
	int numstep=0;
//we do not want to go on forever so max sure we stop after max_step itterations
	for(numstep=0;numstep<max_step;numstep++){

//this makes sure we do not aver shoot our final time 't2'
		if((cur_t+cur_step-t2)*(cur_t+cur_step-t1)>0.0){
			cur_step=t2-cur_t;
		}

//calcualte the yscal factor for this step
		set_yscal();
		//round1();							//do a single step
		if(numstep==0) gr(1);
		else gr(0);
		stepnext=cur_step;
		cur_t+=cur_step;

		if(stepdid==cur_step){
			ngood++;
		}else{
			nbad++;
		}
//if we have finished to t2 break out
		if((cur_t-t2)*(t2-t1)>=0){
			break;
		}

		if(abs(stepnext)<=stepmin){
			std::cout<<"STEP SIZE IS TOO SMALL in odeint...!!"<<std::endl;
			std::cout<<"your step size has exceeded your min step value"<<std::endl;
		}
		if(hold_step){
			if(stepnext>5.*stepi){
				stepnext=5.*stepi;
			}
		}
		cur_step=stepnext;
	}
	if(numstep==max_step){
		std::cout<<"TOO MANY STEPS!!"<<std::endl;
		std::cout<<"i have to stop here, but you can continue from"<<std::endl;
		std::cout<<"the last step if you want...so make sure your"<<std::endl;
		std::cout<<"'front end' can handle it...that is if you want to"<<std::endl;
	}
}

END_BL_NAMESPACE



#endif

