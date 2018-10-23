
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
	stiffbs_meth.h-->Bulirsch-Stoer/Richard Extrapolation method diff eq solver methods

*/

#ifndef _stiffbs_meth_h_
#define _stiffbs_meth_h_ 1

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

static const int __nseqstiff[9]={0,2,6,10,14,22,34,50,70};

template<class func, class T, class Container >
const int  stiffbs<func, T,  Container >::kmaxx=7;

template<class func, class T, class Container>
const int  stiffbs<func, T,  Container >::imaxx=8;


template<class func, class T, class Container >
const double  stiffbs<func, T,Container >::SC_max=1;

//to make compiler happy...they truely set in the set_rest_i() function...
template<class func, class T, class Container >
Container stiffbs<func, T,Container >::yn(1,0);

//to make compiler happy...they truely set in the set_rest_i() function...
template<class func, class T, class Container >
Container stiffbs<func, T,Container >::ym(1,0);

template<class func, class T, class Container >
void stiffbs<func, T,Container >::reset()
{
	basicODE<T, Container>::reset();
	yn.resize(num_func,0);
	ym.resize(num_func,0);
	dfdt.resize(num_func,0);
	x.resize(imaxx,0);
	d.resize(num_func, kmaxx+1, 0);
	found_order=false;
	kmax=0;
	kopt=0;
	found_order=false;
	findorder();
}


template<class func, class T, class Container >
void stiffbs<func,T,Container >::set_yscal()
{
//	yn=(y);		//just a dummy vector..temporary (avoid initialization step for it over and over again)
	my_sys->function(cur_t, y, yscal);
	func_calls++;
//	yscal=tiny +abs(cur_step)*abs(yn)+abs(y);
//	yscal=abs(yn);
//	for(int i=0;i<yscal.size();++i) yscal[i]=max(tiny+abs(y[i])+abs(yscal[i])*abs(cur_step), 1.0);
	for(int i=0;i<yscal.size();++i) yscal[i]=max(max(abs(y[i]-yscal[i])), 1.0);
//	yscal=max(abs(y), 1.0);
}




template<class func, class T, class Container >
int stiffbs<func, T,Container >::mmid(Container &ysavO, double &step, double &nstep, Container &yseqO)
{

//varous temporary numbers i'll need
	double temp_step=0, temp_t=cur_t;
	float d; //used for lu decomp
	int i=0;

//divide our tot steps into the little ones needed
	temp_step=step/nstep;

//scale the jacobian by the current steps size
	JacTemp=-temp_step*jacobi;
	for(i=0;i<JacTemp.rows();++i)  JacTemp(i,i)+=1.0;

	LUwarn=false;
	if(JacTemp.LUdecomp(rowPerm, d)==-1) return 0;
	LU_calls++;
	//JacTemp.inv();

//calculate the f(x,t) at the first y (values held in yseqO=dydt)
	my_sys->function(temp_t, ysavO, yseqO);
	my_sys->function(temp_t+step, ysavO, dfdt); //to get dfdt
	dfdt=(dfdt-yseqO)/step;
	func_calls++;

//now step up the initialstep
	yseqO=temp_step*yseqO+temp_step*temp_step*dfdt;

//LU backsub
	JacTemp.LUbackSub(rowPerm, yseqO);
	yn=yseqO;
	ym=ysavO+yn; //move to the next y

	temp_t+=temp_step; //advance time

	my_sys->function(temp_t, ym, yseqO); //get the next point
	func_calls++;
	for(int nn=1;nn<nstep;++nn){
		yseqO=temp_step*yseqO-yn;
		JacTemp.LUbackSub(rowPerm, yseqO); //get the next solutions
		yn+=2.0*yseqO;
		ym+=yn;
		temp_t+=temp_step;
		my_sys->function(temp_t, ym, yseqO); //get the next point
		func_calls++;

	}
	//the last step
	yseqO=temp_step*yseqO-yn;
	JacTemp.LUbackSub(rowPerm, yseqO);
	yseqO=yseqO+ym;
	LUwarn=true;
	return 1;

}

//this is the polynomial extropoaltion rioutine
//it evaluates num_func functions at x=0 by fitting
//to a series of estimates with progressively smaller values x=xest
//and the function vectors 'yseq.' ('yseq' was given by the
//midpoint routine above). the output is put into
//'y' and the estimated error is put into 'yerr'

template<class func, class T, class Container  >
void stiffbs<func, T,Container >::polyex(int &xiest, T &xest, Container &estim, Container &extrap)
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

/***** FINDING THE ORDER ****/
template<class func, class T, class Container  >
Vector<double> stiffbs<func, T,Container >::a(imaxx+1); //used to detemin the max order

template<class func, class T, class Container  >
rmatrix stiffbs<func, T,Container >::alf(imaxx+1, imaxx+1); //the work matrix

template<class func, class T, class Container  >
Vector<double>  stiffbs<func, T,Container >::err(imaxx+1); //local error

template<class func, class T, class Container  >
void stiffbs<func, T,Container >::findorder()
{
	kmax=0;
	kopt=0;
	stepnext=-1.0e29;
	double eps1=safty1*eps;
	a[1]=__nseqstiff[1]+1;
	for(int k=1;k<=kmaxx;k++){
		a[k+1]=a[k]+__nseqstiff[k+1];
	}
	for(int iq=2;iq<=kmaxx;iq++){
		for(int k=1;k<iq;k++){
			alf(k,iq)=pow(eps1, (a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.)*(2.*k+1)));
		}
	}

//update the work cost for the jacobian evaluation
	a[1]+=num_func;
	for(int k=1;k<=kmaxx;++k) a[k+1]=a[k]+__nseqstiff[k+1];

//figure out the maximum # of stiffbs sections to divide each step into
	for(kopt=2;kopt<kmaxx;kopt++){
		if(a[kopt+1]>a[kopt]*alf(kopt-1,kopt)) {	break; }
	}
	kmax=kopt;
	order=2*(kmax-1);
	found_order=true;
}

template<class func, class T, class Container  >
void stiffbs<func, T,Container >::stiffbsstep(Container &yO, Container &yscalO)
{
	static int  xiest=0;
	static double epsold=-1.0, temp_t=0.;
	static T xest;
	static double temp_step=0;


//	static int reduct, exitflag=0;
	Container yseq=yO;

//set the work matrices if we need to
	if(!found_order || epsold!=eps){ findorder(); epsold=eps;	}

//set the step size to the one given initially
	temp_step=cur_step;

	ysav=yO;
	my_sys->jacobian(cur_t, yO, jacobi);
	jac_calls++;

//	if(cur_t!=temp_t || temp_step !=stepnext){
//		first=1;
//		kopt=kmax;
//	}

//	reduct=0;
	//for(;;){
		for(int k=1;k<=kmax;k++)
		{
			xiest=k;
			temp_t=cur_t+temp_step;
			//my_sys->jacobian(temp_t, yseq, jacobi);
			//jac_calls++;
			//cout<<"Curetn T: " <<cur_t<<" nextstep: "<<stepnext<<endl;;
			if(temp_t==cur_t && temp_t!=t2 && temp_t!=t1 && temp_step==0.0){
				basicODE<T, Container> ::error(std::cerr, "ERROR: stiffbs::bsstep(...)::STEP SIZE UNDERFLOW");
				BLEXCEPTION("ERROR: stiffbs::bsstep(...)::STEP SIZE UNDERFLOW")
				break;
			}
			double nseq_on=double(__nseqstiff[k]);
			//different the the bs one...takes into account the Jacobian

			//if the matrix signular...reduce step size and start again
			while(!mmid(ysav, temp_step, nseq_on, yseq))
			{
				temp_step/=5.0;
				temp_t=cur_t+temp_step;
				k=1;
				nseq_on=double(__nseqstiff[k]);
				xiest=k;
				if(temp_t==cur_t && temp_t!=t2 && temp_t!=t1 && temp_step==0.0){
					basicODE<T, Container> ::error(std::cerr, "ERROR: stiffbs::bsstep(...)::STEP SIZE UNDERFLOW");
					BLEXCEPTION("ERROR: stiffbs::bsstep(...)::STEP SIZE UNDERFLOW")
					break;
				}
			}
		//	error(cout, "\n");
			xest=T(pow(temp_step/nseq_on, 2.0));
			polyex(xiest, xest, yseq, yO);
		}
/*
			if(k!=1){
				errmax=tiny;
				ratio_tmp=max(max(abs(yerr)/abs(yscalO)));
			#ifdef HAVE_ISNAN
				if(hasnan(ratio_tmp)) ratio_tmp=1e5;
			#endif
				errmax=max(errmax, ratio_tmp)/eps;
				km=k-1;
				err[km]=pow(errmax/safty1, 1./(2.*(km)+1));
			//	cout<<"error: "<<errmax<<" nextstep: "<<stepnext<<" temp_step: "<<temp_step<<" scale: "<<scale<<"\n";
			//	cout<<yerr<<endl<<yscalO<<endl;
				//	cout<<"Err Ratio: "<<ratio_tmp<<endl;
			//	cout<<"Err Max: "<<errmax<<endl;
			//	cout<<"Error: "<<err[km]<<endl;
			//	cout<<"k: "<<km<<endl;
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

		if(exitflag)	break;

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
		fact=max(err[k], SC_max);
		work=min(fact*a[k+1]);
		if(work<wrkmin){
			scale=fact;
			wrkmin=work;
			kopt=k+1;
		}
	}
	stepnext=temp_step/scale;

	if(kopt>=xiest && kopt!=kmax && !reduct){
		fact=max(scale/alf(kopt-1,kopt), SC_max);
		if(a[kopt+1]*fact<=wrkmin){
			stepnext=temp_step/fact;
			kopt++;
		}
	}
	//cout<<"Recuction: "<<red<<endl;
		//cout<<"scale: "<<scale<<endl;

*/



//	}
}

template<class func, class T, class Container >
void stiffbs<func,T,Container >::odeint()
{
//set the correect sign of the step (i.e. moving backwards or forwards in time)
	double diff_tmp=t2-t1;
	cur_step=(double(sign(diff_tmp))*stepi);

	int numstep=0;

	//first set a minimum step if not set
	if(stepmin==0.0) stepmin=16.0*eps*abs(t1+tiny);



//we do not want to go on forever so max sure we stop after max_step itterations
//ofstream doodoo("doo", ios::app);
	for(numstep=0;numstep<max_step;numstep++){

//this makes sure we do not aver shoot our final time 't2'
		if((cur_t+cur_step-t2)*(cur_t+cur_step-t1)>0.0){
			cur_step=(t2-cur_t);
		}
//calcualte the yscal factor for this step
		ysav=y;
		//set_yscal();
		stiffbsstep(y, yscal);							//do a single step
/* The standard control object is a four parameter heuristic
 * defined as follows:
 *    D0 = eps_abs + eps_rel * (a_y |y| + a_dydt h |y'|)
 *    D1 = |yerr|
 *    q  = consistency order of method (q=4 for 4(5) embedded RK)
 *    S  = safety factor (0.9 say)
 *
 *                      /  (D0/D1)^(1/(q+1))  D0 >= D1
 *    h_NEW = S h_OLD * |
 *                      \  (D0/D1)^(1/q)      D0 < D1
 */

		double D0=absToler+relToler*(max(max(abs(ysav)))+cur_step*max(max(abs(y))));
		double D1=max(max(abs(yerr)));
		double r=D1/D0;
		double ShrinkForce=1.;
		if(hasnan(r) || hasnan(D1) || D1>1.0 || hasnan(D0)){ r=2.; ShrinkForce=0.1;}
	//	cout<<"r: "<<r<<" "<<D0<<" "<<D1<<endl;
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
			//doodoo<<cur_t<<" "<<y<<endl;

		}else{
			ngood++;
			stepdid=cur_step;
			cur_t+=stepdid;
			//doodoo<<cur_t<<" "<<y<<endl;
		}

		//if(stepdid==cur_step){
		//	ngood++;
		//}else{
		//	nbad++;
		//}
//if we have finished to t2 break out
		if((cur_t-t2)*(t2-t1)>=0){
			//cout<<"cur_t: "<<cur_t<<" t2: "<<t2<<endl;
			break;
		}

		if(abs(stepnext)<=stepmin){
			basicODE<T, Container> ::error(std::cerr,"Error: bs.solve()..STEP SIZE IS TOO SMALL");
			BLEXCEPTION("Error: bs.solve()..STEP SIZE IS TOO SMALL")
		}
		cur_step= stepnext;
	}
	if(numstep==max_step){
		basicODE<T, Container> ::error(std::cerr, "Warning::bs::solve()--TOO MANY STEPS "+itost(numstep));
	}
}


END_BL_NAMESPACE


#endif


