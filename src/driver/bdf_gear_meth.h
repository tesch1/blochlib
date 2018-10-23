


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

#ifndef _bdf_gear_meth_h_
#define _bdf_gear_meth_h_ 1

#include "utils/utils.h"

BEGIN_BL_NAMESPACE


//These are the constants for the integrator
template<class func, class T, class Container >
const int  bdf<func, T,Container >::MaxOrder=5;

template<class func, class T, class Container >
const int  bdf<func, T,Container >::MaxIter=4;

template<class func, class T, class Container >
Vector<double>   bdf<func, T,Container >::G(MaxOrder);

template<class func, class T, class Container >
Vector<double>   bdf<func, T,Container >::invG(MaxOrder);

template<class func, class T, class Container >
Vector<double>  bdf<func, T,Container >::alpha(MaxOrder);

template<class func, class T, class Container >
Vector<double> bdf<func, T,Container >::errConsts(MaxOrder);

template<class func, class T, class Container >
rmatrix bdf<func, T,Container >::difU(MaxOrder,MaxOrder);

template<class func, class T, class Container >
rmatrix  bdf<func, T,Container >::kJ(MaxOrder,MaxOrder);

template<class func, class T, class Container >
rmatrix  bdf<func, T,Container >::kI(MaxOrder,MaxOrder);

template<class func, class T, class Container >
void bdf<func, T,Container >::initConsts(int curorder)
{
	order=curorder;
	double Gs[]={1, 3.0/2.0, 11.0/6.0, 25.0/12.0, 137.0/60.0};
	G=Vector<double>(Gs,MaxOrder);
	double als[]={-37.0/200.0, -1.0/9.0, -0.0823, -0.0415, 0.0};
	alpha=Vector<double>(als,MaxOrder);
	invG=1.0/(G * (1.0 - alpha));
	errConsts=alpha *G + (1.0 /Vector<double>(Spread<double>(2.0,6.0)));
	difU.resize(MaxOrder, MaxOrder, 0.0);
	difU(0,0)=-1; 	difU(0,1)=-2; 	difU(0,2)=-3; 	difU(0,3)=-4; 	difU(0,4)=-5;
	difU(1,0)=0; 	difU(1,1)=1;	difU(1,2)=3; 	difU(1,3)=6; 	difU(1,4)=10;
	difU(2,0)=0; 	difU(2,1)=0;	difU(2,2)=-1; 	difU(2,3)=-4; 	difU(2,4)=-10;
	difU(3,0)=0; 	difU(3,1)=0;	difU(3,2)=0; 	difU(3,3)=1; 	difU(3,4)=5;
	difU(4,0)=0; 	difU(4,1)=0;	difU(4,2)=0; 	difU(4,3)=0; 	difU(4,4)=-1;
	difU.resizeAndPreserve(curorder, curorder);

	kJ.resize(MaxOrder, MaxOrder, 0.0);
	kJ(0,0)=1; 	kJ(0,1)=2; 	kJ(0,2)=3; 	kJ(0,3)=4; 	kJ(0,4)=5;
	kJ(1,0)=1; 	kJ(1,1)=2;	kJ(1,2)=3; 	kJ(1,3)=4; 	kJ(1,4)=5;
	kJ(2,0)=1; 	kJ(2,1)=2;	kJ(2,2)=3; 	kJ(2,3)=4; 	kJ(2,4)=5;
	kJ(3,0)=1; 	kJ(3,1)=2;	kJ(3,2)=3; 	kJ(3,3)=4; 	kJ(3,4)=5;
	kJ(4,0)=1; 	kJ(4,1)=2;	kJ(4,2)=3; 	kJ(4,3)=4; 	kJ(4,4)=5;
	kJ.resizeAndPreserve(curorder, curorder);

	kI.resize(MaxOrder, MaxOrder, 0.0);
	kI(0,0)=1; 	kI(0,1)=1; 	kI(0,2)=1; 	kI(0,3)=1; 	kI(0,4)=1;
	kI(1,0)=2; 	kI(1,1)=2;	kI(1,2)=2; 	kI(1,3)=2; 	kI(1,4)=2;
	kI(2,0)=3; 	kI(2,1)=3;	kI(2,2)=3; 	kI(2,3)=3; 	kI(2,4)=3;
	kI(3,0)=4; 	kI(3,1)=4;	kI(3,2)=4; 	kI(3,3)=4; 	kI(3,4)=4;
	kI(4,0)=5; 	kI(4,1)=5;	kI(4,2)=5; 	kI(4,3)=5; 	kI(4,4)=5;
	kI.resizeAndPreserve(curorder, curorder);

	k=1;
	lastk=1;
	K=Range(0,k-1);
	dif.resize(num_func);
}


template<class func, class T, class Container >
void bdf<func, T,Container >::reset(){
	basicODE<T, Container>::reset();
	got_jacobi=false;
	initConsts(order);
	have_started=false;
	have_rate=false;
}

//this is the 'starter' function for the solver...
// sets the steps sizes, and initial values for
// all of our vaious data containers

template<class func, class T, class Container  >
void bdf<func, T,Container >::start()
{
	//first set a minimum step if not set
	if(stepmin==0.0) stepmin=16.0*eps*abs(t1+tiny);

	//calc the inital jacobain
	my_sys->jacobian(cur_t, y, dfdy);
	got_jacobi=true; //set this so we do not need to calc jacobi at cur_t again
	jac_calls++;

	//need to make a function call to get some inial guesses on the
	// steps sizes needed (the 'yscale')
	my_sys->function(cur_t, y, yscal);
	func_calls++;

	double curThres=absToler/relToler;
	double wt = max(max(abs(y)),curThres);
	double rh = 1.25 * max(max(abs(yscal / wt))) / sqrt(relToler);

	//set the current step size
	cur_step = (stepmax!=0)?min(stepmax, abs(t2 - cur_t)):abs(t2 - cur_t);

	if(cur_step * rh > 1){
		cur_step = 1.0 / rh;
		cur_step = max(cur_step, stepmin);
	}

//set the coreect sign of the step (i.e. moving backwards or forwards in time)
	tdirection=double(sign(t2-cur_t));

// The error of BDF1 is 0.5*h^2*y''(t), so we can determine the optimal h.
	abs_step=cur_step;
	cur_step=tdirection*cur_step;
//out first 'difference' point
	double dt = (cur_t + tdirection*min(sqrt(eps)*max(abs(cur_t),abs(cur_t+cur_step)),abs_step) - cur_t);

//the ysav is simply used as a continer ...pay no attention to its name
	my_sys->function(cur_t+dt,y,ysav);
	func_calls++;

//the yerr is simply used as a continer ...pay no attention to its name
	yerr = (ysav - yscal) / dt;

	rh = 1.25 * sqrt(0.5 * max(max(abs(yerr + dfdy*yscal))) / wt / relToler);
	abs_step = min(stepmax, abs(t2- cur_t));
	if(abs_step * rh > 1)  abs_step = 1.0 / rh;

	abs_step = max(abs_step, stepmin);
	cur_step=tdirection*abs_step;

	// Initialize.
	k = 1;                                  // start at order 1 with BDF1
	K=Range(0,k-1);                                // K = 1:k
	lastk = k;
	stepdid = cur_step;

	dif.resize(order+2, Container(num_func, 0));
	dif[k-1]=cur_step * yscal;

	JacTemp = cur_step * invG(k-1)*dfdy;
	JacTemp.LU(L,U);
	LU_calls++;

	cout<<"Current step: "<<cur_step<<endl;
	cout<<JacTemp<<L<<U<<endl;
	ysav.print(cout,"\n");cout<<endl;
	invG.print(cout,"\n"); cout<<endl;

	have_started=true;
}



template<class func, class T, class Container >
void bdf<func,T,Container >::odeint()
{
	int numstep=0;
	if(!have_started) start();
	double d;
  	mattype difRU, tmpRhs;
  	int i,j;
  	bool done=false;
//		difRU=kI-1.0-kJ*abs(cur_step)/abs(stepdid);
//		cout<<"difRU: "<<endl<<difRU;
//Stretch the step if within 10% of t2-t.
	if(1.1*abs_step >= abs(t2 - cur_t))
	{
		cur_step = t2 - cur_t;
		abs_step = abs(cur_step);
		done = true;
	}

//update our differences

 	if((cur_step != stepdid) || (k != lastk))
  	{
//		cout<<"difRU: "<<endl<<kI-kJ;
		difRU=kI-1.0-kJ*abs(cur_step)/abs(stepdid);
//		cout<<"difRU: "<<endl<<difRU;
		for(i=1;i<difRU.rows();++i)
			for(j=0;j<difRU.cols();++j)
				difRU(i,j)*=difRU(i-1,j);

//		cout<<"cumprod(difRU): "<<endl<<difRU;

		for(i=0;i<difRU.rows();++i)
			for(j=0;j<difRU.cols();++j)
				difRU(i,j)/=kI(i,j);

//		cout<<"difRU/kI: "<<endl<<difRU;

		difRU *=difU;
//		cout<<"difRU*difU: "<<endl<<difRU;
//		cout<<"K: "<<endl<<dif(K)<<" ||"<<K<<" ||"<<difRU(K,K)<<endl;;
		dif(K)=dif(K)* difRU(K,K);

		JacTemp = cur_step * invG(k-1) * dfdy;
		JacTemp.LUdecomp(rowPerm, d);
		LU_calls++;
		have_rate = false;
	}

	bool  notfailed = true;
	Container psi(num_func, 0.0);
	Container difpred(num_func, 0.0);

	double curThres=absToler/relToler;
	double newnrm, oldnrm,errit,minnrm,invwt;

	double tnew;
	while(1)  //main loop to do out total time step....t1-->t2
	{
		bool gotynew = false;                    // is ysav evaluated yet?
		while(!gotynew)
		{
		// Compute the constant terms in the equation for ynew.
			psi.fill(ZeroType<T>::zero());
			for(int o=0;o<k;++o)
				for(int p=0;p<psi.size();++p)
					psi(p) += dif[o][p] * (G(o) * invG(k));

			psi.print(cout, "\n"); cout<<endl;     BLEXCEPTION(__FILE__,__LINE__)

		// Predict a solution at cur_t+cur_step.
			tnew = cur_t + cur_step;
			if(done)	tnew = t2;   //go to the exact end.
			ytmp = y + sum(dif(K),2); //our prediction
			ysav = ytmp;

		// The difference, difpred, between pred and the final accepted
		// ysav is equal to the backward difference of ysaV of order
		// k+1. Initialize to zero for the iteration to compute ysav.
			difpred.fill(0.0);
			invwt = 1.0 / max(max(max(max(abs(y))),max(max(abs(ysav)))),curThres);
			minnrm = 100.0*eps*max(max(abs(ysav * invwt)));

	      // Iterate with simplified Newton method. to find 'ysav' (or next step)
			bool tooslow = false;
			int iter=0;
			for(iter=0; ite<MaxIter;++iter)
			{
				my_sys->function(tnew, ysav,yscal);
				yscal = cur_step*invG(k-1)*yscal -  (psi+difpred);
       			JacTemp.LUbackSub(rowPerm, yscal);

				newnrm = max(max(abs(yscal)))* invwt;

				difpred += yscal;
	        	ysav = ytmp + difpred;

				if(newnrm <= minnrm){
					gotynew = true;
					break;
				}else if(iter == 1){
					if(have_rate)
					{
						errit = newnrm * rate / (1.0 - rate);
						if(errit <= 0.05*relToler) // More stringent when using old rate.
						{
							gotynew = true;
							break;
						}else{
							rate = 0;
						}
					}
				}else if(newnrm > 0.9*oldnrm){
					tooslow = true;
					break;
				}else{
					rate = max(0.9*rate, newnrm / oldnrm);
					have_rate = true;
					errit = newnrm * rate / (1.0 - rate);
					if( errit <= 0.5*relToler){
	            		gotynew = true;
	            		break;
	          		}else if( iter == MaxIter){
	            		tooslow = true;
						break;
					}else if(0.5*relToler < errit*pow(rate,double(MaxIter-iter))){
						tooslow = true;
						break;
					}
				}
				oldnrm = newnrm;
			}       // end of Newton loop
			func_calls+= iter;

			if(tooslow)
			{
				nbad++;
	        // Speed up the iteration by forming new linearization or reducing cur_step.
				if(!got_jacobi){
					my_sys->jacobian(cur_t,y,dfdy);
					jac_calls++;
					got_jacobi = true;

				}else if( abs_step <= stepmin){
				  basicODE<T, Container>::error(std::cerr, "Error: bdf.odeint()...STEP SIZE TOO SMALL");
				  return;
				}else{
					stepdid = cur_step;
					abs_step = max(0.3 * abs_step, stepmin); //reduce dramtically the stepsize
					cur_step = tdirection * abs_step;
					done = false;
					difRU=kI-1.0-kJ*abs(cur_step)/abs(stepdid);
					for(i=1;i<difRU.rows();++i)
						for(j=0;j<difRU.cols();++j)
							difRU(i,j)*=difRU(i-1,j);

					for(i=0;i<difRU.rows();++i)
						for(j=0;j<difRU.cols();++j)
							difRU(i,j)/=kI(i,j);

					difRU *=difU;
					dif(K)=dif(K)* difRU(K,K);
				}
				JacTemp =  - cur_step*invG(k-1) * dfdy;
				JacTemp.LUdecomp(rowPerm,d);
				LU_calls++;
				have_rate = false;
			} 				//End "Too Slow"
		}    // end of while loop for getting the New Y (ysav)

	    //difpred is now the backward difference of ysav of order k+1.
	    double err = max(max(abs(difpred)) )* invwt * errConsts(k-1);

		if( err > relToler )	// Failed step error is to big
		{
			nbad++;
			if(abs_step <= stepmin) //step size too small
			{
      			basicODE<T, Container>::error(std::cerr, "Error: bdf.odeint()...Cannot meet relToler requirements");
				return;
			}

			stepdid = cur_step;
			if(notfailed)
			{
				notfailed = false;
				double step_opt = abs_step * max(0.1, 0.833*pow(relToler/err, 1.0/(double(k)+1.0)) ); //shrink by 17%
				if( k > 1) //our backward 'error' control
				{
					double errkm1 = max(max(abs((dif(k-1) + difpred)))) * invwt * errConsts(k-1);
					double step_km1 = abs_step * max(0.1, 0.769*pow(relToler/errkm1,1.0/double(k) )); // the back step
					if(step_km1 > step_opt)
					{
						step_opt = min(abs_step,step_km1);      // don't allow step size increase
						k = k - 1;
						K = Range(0, k-1);
					}
				}
				abs_step = max(stepmin, step_opt);
			}else{
				abs_step = max(stepmin, 0.5 * abs_step);
			}

			cur_step = tdirection * abs_step;

			if(abs(abs_step) < abs(stepdid))	done = false;

			difRU=kI-1.0-kJ*abs(cur_step)/abs(stepdid);

			for(i=1;i<difRU.rows();++i)
				for(j=0;j<difRU.cols();++j)		difRU(i,j)*=difRU(i-1,j);

			for(i=0;i<difRU.rows();++i)
				for(j=0;j<difRU.cols();++j)		difRU(i,j)/=kI(i,j);

			difRU *=difU;
			dif(K)=dif(K)* difRU(K,K);
			JacTemp =  - cur_step * invG(k-1) * dfdy;
			JacTemp.LUdecomp(rowPerm,d);
			LU_calls++;
			have_rate = false;

		}else{                               // Successful step
			break;
		}
	}	//End main loop to do out total time step....t1-->t2 Successful step
	y=ysav; //place the 'good point' as the new y value

/*
//we do not want to go on forever so max sure we stop after max_step itterations
	for(numstep=0;numstep<max_step;numstep++){

//this makes sure we do not aver shoot our final time 't2'
		if((cur_t+cur_step-t2)*(cur_t+cur_step-t1)>0.0){
			cur_step=(t2-cur_t);
		}
//calcualte the yscal factor for this step
		set_yscal();
		bsstep(y, yscal);							//do a single step
		if(stepdid==cur_step){
			ngood++;
		}else{
			nbad++;
		}
//if we have finished to t2 break out
		if((cur_t-t2)*(t2-t1)>=0){
			//cout<<"cur_t: "<<cur_t<<" t2: "<<t2<<endl;
			break;
		}

		if(abs(stepnext)<=stepmin){
			basicODE<T, Container> ::error(std::cerr,"Error: bdf.solve()..STEP SIZE IS TOO SMALL");
			    BLEXCEPTION(__FILE__,__LINE__)
		}
		if(hold_step){
			if(stepnext>5.*stepi){
				stepnext=5.*stepi;
			}
		}

		cur_step= stepnext;
	}
	if(numstep==max_step){
		basicODE<T, Container> ::error(std::cerr, "Warning::bdf::solve()--TOO MANY STEPS "+itost(numstep));
	}
*/
}

END_BL_NAMESPACE

#endif


