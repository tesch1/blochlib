

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
	ckrk_meth.h-->a 5th order Runge-Kutta Diff eq solver methods

*/

#ifndef _ckrk_meth_h_
#define _ckrk_meth_h_ 1


BEGIN_BL_NAMESPACE

template<class func, class T, class Container>
const float  ckrk<func, T, Container>::a2=0.2;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::a3=0.3;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::a4=0.6;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::a5=1.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::a6=0.875;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b21=0.2;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::b31=3.0/40.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b32=9.0/40.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b41=0.3;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b42 = -0.9;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b43=1.2;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::b51 = -11.0/54.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b52=2.5;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b53 = -70.0/27.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b54=35.0/27.0;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::b61=1631.0/55296.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b62=175.0/512.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b63=575.0/13824.0;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::b64=44275.0/110592.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::b65=253.0/4096.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::c1=37.0/378.0;
template<class func, class T, class Container>
 const float  ckrk<func, T, Container>::c3=250.0/621.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::c4=125.0/594.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::c6=512.0/1771.0;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::dc1= ckrk<func, T, Container>::c1-2825.0/27648.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::dc3= ckrk<func, T, Container>::c3-18575.0/48384.0;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::dc5 = -277.00/14336.0;
template<class func, class T, class Container>
const float  ckrk<func, T, Container>::dc4= ckrk<func, T, Container>::c4-13525.0/55296.0;
template<class func, class T, class Container>
const float ckrk<func, T, Container>::dc6= ckrk<func, T, Container>::c6-0.25;

template<class func, class T, class Container>
double ckrk<func, T, Container>::grow=-0.2;

template<class func, class T, class Container>
double ckrk<func, T, Container>::shrink=-0.25;

template<class func, class T, class Container>
int ckrk<func, T, Container>::max_step_within=10000;



/*Container ckrk::set_y_mat(matrix &iny){
	Container tmp;
	tmp=my_sys->array_from_mat(iny);
	return tmp;
}*/

template<class func, class T, class Container>
void ckrk<func, T, Container>::set_yscal(){
	Container tmp(num_func, T(0));
	//float floo=0;
	my_sys->function(cur_t, y, tmp);
	func_calls++;
	/*for(int i=0;i<num_func;i++){
		if(abs(y[i])>=0.1)	{
			floo=abs(cur_step)*abs(tmp[i])+tiny+abs(y[i]);
			yscal[i]=floo;
		}else{
			floo=abs(cur_step)*abs(tmp[i])+tiny+abs(y[i]);
			yscal[i]=floo;
		}
	}*/
	yscal=abs(cur_step)*abs(tmp)+tiny+abs(y);

}

template<class func, class T, class Container>
void ckrk<func, T, Container>::rk(){

//here are some temp bins for the 'for' loop
	static Container ak1(num_func,0.);
	static Container ak2(num_func,0.);
	static Container ak3(num_func,0.);
	static Container ak4(num_func,0.);
	static Container ak5(num_func,0.);
	static Container ak6(num_func,0.);
	static Container ytemp(num_func,0.);
	static Container ytemp1(num_func,0.);
	static double timestep=0;
	ytemp=y;
//the first one in the 6 step part is getting k1=step*f(tn, yn)
//it is done this way (differeent from the others, becuase
//we have a bondray value problem that we hopefully entered
//at the first my_sys->function call...i.e at t=t1, y1=y(t1) and f(t1, y1)
//it gets advanced at the end of the sequence

//the next 5 sets relay on the one previous
//here we move t a bit, use our old values * some constant
//
//the one commented out below is for n number of NON coupled ODE's
//not really usefull for most things, but it is here in case you need
//to use it...it is faster then the next one

	/*for(int i=0;i<num_func;i++){
		cury=y;
		ak1[i]=cur_step*my_sys->function::func(cur_t, cury->num,i);
		ak2[i]=cur_step*my_sys->function::func(cur_t+a2*cur_step, cury->num+b21*ak1[i],i);
		ak3[i]=cur_step*my_sys->function::func(cur_t+a3*cur_step, cury->num+b31*ak1[i]+b32*ak2[i], i);
		ak4[i]=cur_step*my_sys->function::func(cur_t+a4*cur_step, cury->num+b41*ak1[i]+b42*ak2[i]+b43*ak3[i], i);
		ak5[i]=cur_step*my_sys->function::func(cur_t+a5*cur_step, cury->num+b51*ak1[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i], i);
		ak6[i]=cur_step*my_sys->function::func(cur_t+a6*cur_step, cury->num+b61*ak1[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i],i);
		//std::cout<<"f("<<cur_t<<","<<cury->num<<") = "<<ak1[i]/cur_step<<std::endl;
		cury=cury->nextaddr;
	}*/

	//the coupling is done in the 'my_sys->function'
	//class soooo, make sure it is correct...it goes through
	//the sam rotine as above except that it is a bit more difficult to
	//call a coupled equation if you are looping through one variable at a time
	//so we need to send the entire variable package to the my_sys->function
	//and afterwards creat the next variable set (i.e. ytemp1)
	//send it off...and so on
	//
		my_sys->function(cur_t, ytemp, ak1);


		/*for(int i=0;i<num_func;i++){
			ak1[i]*=cur_step;
			ytemp1[i]=ytemp[i]+b21*ak1[i];
		}*/

		ak1*=cur_step;
		ytemp1=ytemp+b21*ak1;

		timestep=	cur_t+a2*cur_step;
		my_sys->function(timestep, ytemp1, ak2);

		/*for(int i=0;i<num_func;i++){
			ak2[i]=ak2[i]*cur_step;
			ytemp1[i]=ytemp[i]+b31*ak1[i]+b32*ak2[i];
		}*/
		ak2*=cur_step;
		ytemp1=ytemp+b31*ak1+b32*ak2;

		timestep=cur_t+a3*cur_step;
		my_sys->function(timestep, ytemp1, ak3);


		/*for(int i=0;i<num_func;i++){
			ak3[i]=ak3[i]*cur_step;
			ytemp1[i]=ytemp[i]+b41*ak1[i]+b42*ak2[i]+b43*ak3[i];
		}*/
		ak3*=cur_step;
		ytemp1=ytemp+b41*ak1+b42*ak2+b43*ak3;

		timestep=cur_t+a4*cur_step;
		my_sys->function(timestep, ytemp1, ak4);


		/*for(int i=0;i<num_func;i++){
			ak4[i]=ak4[i]*cur_step;
			ytemp1[i]=ytemp[i]+b51*ak1[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i];
		}*/
		ak4*=cur_step;
		ytemp1=ytemp+b51*ak1+b52*ak2+b53*ak3+b54*ak4;

		timestep=cur_t+a5*cur_step;
		my_sys->function(timestep, ytemp1, ak5);
		/*for(int i=0;i<num_func;i++){
			ak5[i]=ak5[i]*cur_step;
			ytemp1[i]=ytemp[i]+b61*ak1[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i];
		}*/
		ak5*=cur_step;
		ytemp1=ytemp+b61*ak1+b62*ak2+b63*ak3+b64*ak4+b65*ak5;

		timestep=cur_t+a6*cur_step;
		my_sys->function(timestep, ytemp1, ak6);

		/*for(int i=0;i<num_func;i++){
			y[i]=y[i]+(c1*ak1[i]+c3*ak3[i]+c4*ak4[i]+cur_step*c6*ak6[i]);
			yerr[i]=(dc1*ak1[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+cur_step*dc6*ak6[i]);
		}*/

		y+=(c1*ak1+c3*ak3+c4*ak4+cur_step*c6*ak6);
		yerr=(dc1*ak1+dc3*ak3+dc4*ak4+dc5*ak5+cur_step*dc6*ak6);
		func_calls+=6;

}

template<class func, class T, class Container>
void ckrk<func, T, Container>::round1(){
	double temp_err=0;
	double temp_step=0;
	double temp_t=0;
	int count=0;
	static double err_max=eps/2.0;

	while(count<max_step_within){
		rk();

		temp_err=0;
		/*for(int i=0;i<num_func;i++){
			float loof=yerr[i]/yscal[i];
			loof=abs(loof);
			temp_err=max(temp_err, loof);
		}*/

		temp_err=max(temp_err,max(abs(yerr/yscal)) );

		temp_err=temp_err/eps;
		if(temp_err<=1.0)	break;

		temp_step=safty1*cur_step*pow(temp_err, shrink);
		double temp_compar=0.1*cur_step;
		if(cur_step>=0){
			cur_step=max(temp_step, temp_compar);
		}else{
			cur_step=min(temp_step, temp_compar);
		}
		temp_t=cur_t+cur_step;
		if(temp_t==cur_t&&cur_step==0.){
			basicODE<T, Container>::error(std::cerr, "Error::ckrk::round1()...STEP SIZE UNDERFLOW (upping step size)");
			temp_t*=2;
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

template<class func, class T, class Container>
void ckrk<func, T, Container>::odeint(){
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
		round1();							//do a single step


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
			basicODE<T, Container>::error(std::cerr,"Error: ckrk::odeint()...STEP SIZE IS TOO SMALL");
			BLEXCEPTION("Error: ckrk::odeint()...STEP SIZE IS TOO SMALL")
		}
		if(hold_step){
			if(stepnext>5.*stepi){
				stepnext=5.*stepi;
			}
		}
		cur_step=stepnext;
	}
	if(numstep==max_step){
		basicODE<T, Container> ::error(std::cerr, "Warning::bs::solve()--TOO MANY STEPS "+itost(numstep));
	}
}


END_BL_NAMESPACE
#endif

